import os
import sys
from pathlib import Path

inversefold_path = Path(__file__).resolve().parents[4] / "examples" / "inverse_folding"
os.chdir(inversefold_path)

sys.path.insert(0, str(inversefold_path))

import torch
from esm import FastaBatchedDataset

import socket
import time
from pydantic import BaseModel
from fastapi import HTTPException

from comps import (
    CustomLogger,
    register_microservice,
    register_statistics,
    opea_microservices,
    statistics_dict,
)
import esm
from esm import pretrained, MSATransformer
import base64
import json
import uuid
from typing import List, Dict, Any,Optional
from comps import OpeaComponent, OpeaComponentRegistry

from biotite.sequence.io.fasta import FastaFile, get_sequences
import numpy as np
import torch.nn.functional as F
from tqdm import tqdm
import esm.inverse_folding

import score_log_likelihoods
from types import SimpleNamespace



def store_client_input(input_pdb_str, input_fasta_str, server_id: int = 0):
    try:
        work_dir = f"/tmp/inverse_fold_score/server_{server_id}"
        os.makedirs(work_dir, exist_ok=True)
        # Handle PDB input
        if input_pdb_str:
            if not input_pdb_str.strip():
                raise HTTPException(status_code=400, detail="inference.input_pdb is empty")
            input_pdb_file = os.path.join(work_dir, "input.pdb")
            with open(input_pdb_file, "w") as f:
                f.write(input_pdb_str + "\n")
            logger.info("Saved input_pdb to input.pdb")

        # Handle FASTA input (base64-encoded)
        if input_fasta_str:
            decoded = base64.b64decode(input_fasta_str)
            input_fasta_file = os.path.join(work_dir, "input.fasta")
            with open(input_fasta_file, "wb") as f:
                f.write(decoded)
            logger.info(f"Decoded base64 FASTA and saved to {input_fasta_file}")
        else:
            raise HTTPException(status_code=400, detail="FASTA base64 string is missing.")

    except HTTPException as http_exc:
        raise http_exc

    except Exception as e:
        logger.error(f"Error in ESM - Embedding microservice: {e}")
        raise HTTPException(status_code=500, detail=str(e))

    return input_pdb_file,input_fasta_file



def validate_client_inference_input(input_pdb_str):
    try:
        
        if input_pdb_str is None:
            logger.error("Received null PDB string in request.")
            raise HTTPException(
                status_code=400,
                detail="Missing required input: input_pdb_str"
            )

    except HTTPException as http_exc:
        raise http_exc

    except Exception as e:
        logger.error(f"Error in ESM - Embedding microservice: {e}") 
        raise HTTPException(status_code=500, detail=str(e))

def get_machine_ip():
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(("8.8.8.8", 80))
        ip = s.getsockname()[0]
        s.close()
        return ip
    except Exception:
        return "127.0.0.1"

WORKLOAD_NAME = "esm_inversefold"
logger = CustomLogger(f"{WORKLOAD_NAME}_microservice")

class ESM_inversefoldscoreInput(BaseModel):
    pdb_str: Optional[str] = None
    fasta_str: Optional[str] = None
    chain: Optional[str] = None  # new
    multichain_backbone: Optional[bool] = False
    singlechain_backbone: Optional[bool] = False
    port: Optional[int] = None

class ESM_inversefoldscoreOutput(BaseModel):
    status: str
    result: Optional[str]=None
    message: Optional[str] = None

@OpeaComponentRegistry.register("OPEA_OMICS_esminversefold_score")
class Opea_ESM_inversefold_score(OpeaComponent):
    
    def __init__(self, name: str, description: str, config: dict = None):
        super().__init__(name, description, config)

        self.use_bf16 = False
        if config.get("bf16"):
            self.use_bf16 = config["bf16"]
            logger.info(f"Info: Running in bfloat16 mode.")

        model_load_time = time.time()
        self.model, self.alphabet = self._load_model()
        print("Model load time",time.time() - model_load_time)
        
    
    def _load_model(self):
        
        model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
        model = model.eval()
        
        if self.use_bf16:
            model=model.bfloat16()
            
        return model,alphabet
    
    def _run_inference(self,model,alphabet,pdbfile,seqfile,chain,multichain_backbone):
        args = SimpleNamespace(
            pdbfile=pdbfile,
            seqfile=seqfile,
            chain=chain,
            multichain_backbone=multichain_backbone,
            outpath=None
        )
        
        device_type ="cpu"
        with torch.no_grad():
            with torch.amp.autocast(device_type=device_type , enabled=self.use_bf16):
                
                if multichain_backbone:
                    sampled_score = score_log_likelihoods.score_multichain_backbone(model, alphabet, args)
                else:
                    sampled_score =score_log_likelihoods.score_singlechain_backbone(model, alphabet,args)

        return sampled_score
    
        
    
    async def invoke(self, input: ESM_inversefoldscoreInput) -> ESM_inversefoldscoreOutput:

        try:
            start_time = time.time()
            
            validate_client_inference_input(input.pdb_str)
            
            store_input_pdb,store_input_fasta  = store_client_input(input.pdb_str,input.fasta_str,input.port)

            results = self._run_inference(self.model, self.alphabet, store_input_pdb,store_input_fasta,input.chain,input.multichain_backbone)

            result_encoded = base64.b64encode(json.dumps(results).encode("utf-8")).decode("utf-8")

            print("total run time",time.time() - start_time)
            return ESM_inversefoldscoreOutput(
                status="success",
                result=result_encoded,
                message="Inference completed successfully."
            )

        except HTTPException as http_exc:
            raise http_exc

    def check_health(self) -> bool:
        """Checks the health of the ESM2 InverseFold Score service."""
        if self.model is None:
            logger.error("ESM2 InverseFold Score model is not initialized.")
            return False

        return True