# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import os
import sys
from pathlib import Path

fold_path = Path(__file__).resolve().parents[3] / "scripts"
os.chdir(fold_path)
sys.path.insert(0, str(fold_path))

from esm import FastaBatchedDataset
import os
import socket
import argparse
import time
from typing import List, Optional
from pydantic import BaseModel
from fastapi import HTTPException

from comps import CustomLogger

import torch
from esm import pretrained, MSATransformer
import base64
import uuid
from typing import Dict, Any
from comps import CustomLogger, OpeaComponent, OpeaComponentRegistry

from pathlib import Path
import sys
import argparse
import logging
import typing as T
from pathlib import Path
from timeit import default_timer as timer

import torch

import esm
from esm.data import read_fasta
from fold import enable_cpu_offloading,init_model_on_gpu_with_cpu_offloading,create_batched_sequence_datasest



logger = logging.getLogger()
logger.setLevel(logging.INFO)

formatter = logging.Formatter(
    "%(asctime)s | %(levelname)s | %(name)s | %(message)s",
    datefmt="%y/%m/%d %H:%M:%S",
)

console_handler = logging.StreamHandler(sys.stdout)
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(formatter)
logger.addHandler(console_handler)


PathLike = T.Union[str, Path]


def store_client_input(fasta_base64,server_id: int = 0):
    try:
        work_dir = f"/tmp/esmfold_input/server_{server_id}"
        os.makedirs(work_dir, exist_ok=True)
        # Support base64 input
        if fasta_base64:
            decoded = base64.b64decode(fasta_base64)
            input_file = os.path.join(work_dir, "input.fasta")
            with open(input_file, "wb") as f:
                f.write(decoded)
            logger.info(f"Decoded base64 FASTA and saved to {input_file}")
        else:
            raise HTTPException(status_code=400, detail="FASTA base64 string is missing.")
        
    except HTTPException as http_exc:
        raise http_exc

    except Exception as e:
        logger.error(f"Error in ESMFold microservice: {e}")
        raise HTTPException(status_code=500, detail=str(e))
    return input_file


def validate_client_inference_input(fasta_base64):
    try:
        
        if fasta_base64 is None:
            logger.error("Received null fasta_base64 in request.")
            raise HTTPException(
                status_code=400,
                detail="Missing required input: fasta_base64"
            )

    except HTTPException as http_exc:
        raise http_exc

    except Exception as e:
        logger.error(f"Error in ESMFold microservice: {e}") 
        raise HTTPException(status_code=500, detail=str(e))



WORKLOAD_NAME = "esm_fold"
logger = CustomLogger(f"{WORKLOAD_NAME}_microservice")

model_path = None
model = None
alphabet = None

class ESM_foldInput(BaseModel):
    fasta_base64: str
    max_tokens_per_batch: Optional[int] = 4096
    num_recycles: Optional[List[int]] = [-1]
    output_dir: Optional[str] = "output"
    port: Optional[int] = None

class ESM_foldOutput(BaseModel):
    status: str
    result: Optional[Dict[str, Any]] = None
    message: Optional[str] = None

@OpeaComponentRegistry.register("OPEA_OMICS_ESMFOLD")
class Opea_ESM_fold(OpeaComponent):

    def __init__(self, name: str, description: str, config: dict = None):
        super().__init__(name, description, config)
        
        logger.info("Loading model")
        self.use_bf16 = False
        if config.get("bf16"):
            self.use_bf16 = config["bf16"]
            logger.info(f"Info: Running in bfloat16 mode.")
        model_load_time = time.time()
        self.model = self._load_model(config["model_dir"],config["chunk_size"],config["cpu_offload"])
        print("Model load time",time.time() - model_load_time)
        
    def _load_model(self,model_dir,chunk_size,cpu_offload):
        if model_dir is not None:
            # if pretrained model path is available
            torch.hub.set_dir(model_dir)

        model = esm.pretrained.esmfold_v1()


        model = model.eval()
        model.set_chunk_size(chunk_size)
        noipex = True
        cpu_only =True
        if cpu_only:
            model.esm.float()  # convert to fp32 as ESM-2 in fp16 is not supported on CPU
            model.cpu()
        elif cpu_offload:
            model = init_model_on_gpu_with_cpu_offloading(model)
        else:
            model.cuda()
        if not noipex:
            dtype = torch.bfloat16 if self.use_bf16 else torch.float32
            import intel_extension_for_pytorch as ipex
            model = ipex.optimize(model, dtype=dtype)
        if noipex and self.use_bf16:
            model=model.bfloat16()
        return model
        
    def _run_inference(self,fasta_path,max_tokens_per_batch,num_recycles):
        all_sequences = sorted(read_fasta(fasta_path), key=lambda header_seq: len(header_seq[1]))
        logger.info("Starting Predictions")
        batched_sequences = create_batched_sequence_datasest(all_sequences, max_tokens_per_batch)

        num_completed = 0
        num_sequences = len(all_sequences)
        for headers, sequences in batched_sequences:
            start = timer()
            try:
                output = self.model.infer(sequences, num_recycles=num_recycles)
            except RuntimeError as e:
                if e.args[0].startswith("CUDA out of memory"):
                    if len(sequences) > 1:
                        logger.info(
                            f"Failed (CUDA out of memory) to predict batch of size {len(sequences)}. "
                            "Try lowering `--max-tokens-per-batch`."
                        )
                    else:
                        logger.info(
                            f"Failed (CUDA out of memory) on sequence {headers[0]} of length {len(sequences[0])}."
                        )

                    continue
                raise

            folding_result = {}  # Dictionary to hold {fastaname: pdb_string}

            output = {key: value.cpu() for key, value in output.items()}
            pdbs =  self.model.output_to_pdb(output)
            tottime = timer() - start
            time_string = f"{tottime / len(headers):0.1f}s"
            if len(sequences) > 1:
                time_string = time_string + f" (amortized, batch size {len(sequences)})"

            for header, seq, pdb_string, mean_plddt, ptm in zip(
                headers, sequences, pdbs, output["mean_plddt"], output["ptm"]
                    ):
                # Store result in dictionary instead of writing to file
                folding_result[header] = pdb_string

                num_completed += 1
                logger.info(
                    f"Predicted structure for {header} with length {len(seq)}, pLDDT {mean_plddt:0.1f}, "
                    f"pTM {ptm:0.3f} in {time_string}. "
                    f"{num_completed} / {num_sequences} completed."
                )

            return folding_result  # Returns dict: {fastaname: pdb_string}

        
    async def invoke(self, input: ESM_foldInput) -> ESM_foldOutput:

        try:
            start_time = time.time()
            validate_client_inference_input(input.fasta_base64)
            
            store_input_fasta  = store_client_input(input.fasta_base64,input.port)

            results = self._run_inference( store_input_fasta, input.max_tokens_per_batch,input.num_recycles)
            
            print("total run time",time.time() - start_time)
            return ESM_foldOutput(
                status="success",
                result=results,
                message="Inference completed successfully."
            )

        except HTTPException as http_exc:
            raise http_exc



    def check_health(self) -> bool:
        """Checks the health of the ESM2 fold service."""
        if self.model is None:
            logger.error("ESM2 fold model is not initialized.")
            return False

        return True
