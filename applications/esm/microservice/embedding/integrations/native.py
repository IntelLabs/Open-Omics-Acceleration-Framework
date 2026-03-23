# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

from pathlib import Path
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


def store_client_input(fasta_base64, server_id: int = 0):
    try:
        work_dir = f"/tmp/embedding_inputs/server_{server_id}"
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
        logger.error(f"Error in ESM - Embedding microservice: {e}")
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
        logger.error(f"Error in ESM - Embedding microservice: {e}") 
        raise HTTPException(status_code=500, detail=str(e))



WORKLOAD_NAME = "esm_embbeding"
logger = CustomLogger(f"{WORKLOAD_NAME}_microservice")

model_path = None
model = None
alphabet = None

class ESM_embeddingInput(BaseModel):
    model_name: Optional[str] = "esm2_t6_8M_UR50D"
    fasta_base64: str
    toks_per_batch: Optional[int] = 4096
    repr_layers: Optional[List[int]] = [-1]
    truncation_seq_length: Optional[int] = 1022
    include: Optional[List[str]] = ["mean"]
    bf16: Optional[bool] = False
    output_dir: Optional[str] = "output"
    port: Optional[int] = None

class ESM_embeddingOutput(BaseModel):
    status: str
    result: Optional[List[Dict[str, Any]]] = None
    message: Optional[str] = None

@OpeaComponentRegistry.register("OPEA_OMICS_esmembedding")
class Opea_ESM_embedding(OpeaComponent):

    def __init__(self, name: str, description: str, config: dict = None):
        super().__init__(name, description, config)

        model_name_or_path = os.getenv("MODEL") or config.get("model_name") if config else None

        if not model_name_or_path:
            model_name_or_path = "esm2_t33_650M_UR50D"

        self.use_bf16 = False
        if config.get("bf16"):
            self.use_bf16 = config["bf16"]
            logger.info(f"Info: Running in bfloat16 mode.")

        print(f"Using model: {model_name_or_path}")
        model_load_time = time.time()
        self.model, self.alphabet = self._load_model(model_name_or_path)
        print("Model load time",time.time() - model_load_time)
    
    def _load_model(self,model_path):
        model, alphabet = pretrained.load_model_and_alphabet(model_path)
        model.eval()
        if self.use_bf16:
            model.bfloat16()

        if isinstance(model, MSATransformer):
            raise ValueError(
                "This script currently does not handle models with MSA input (MSA Transformer)."
            )
        if torch.cuda.is_available():
            model = model.cuda()
            print("Transferred model to GPU")
            
        return model, alphabet

    def _run_inference(self,model, alphabet, fasta_path, toks_per_batch,truncation_seq_length,include,repr_layers):
        dataset = FastaBatchedDataset.from_file(fasta_path)
        batches = dataset.get_batch_indices(toks_per_batch, extra_toks_per_seq=1)
        data_loader = torch.utils.data.DataLoader(
            dataset,
            collate_fn=alphabet.get_batch_converter(truncation_seq_length),
            batch_sampler=batches,
        )
        print(f"Read {fasta_path} with {len(dataset)} sequences")

        return_contacts = "contacts" in include
        enable_autocast = self.use_bf16
        device_type = "cpu" 

        assert all(-(model.num_layers + 1) <= i <= model.num_layers for i in repr_layers)
        repr_layers = [(i + model.num_layers + 1) % (model.num_layers + 1) for i in repr_layers]

        embedding_results = []

        with torch.no_grad():
            for batch_idx, (labels, strs, toks) in enumerate(data_loader):
                print(f"Processing {batch_idx + 1} of {len(batches)} batches ({toks.size(0)} sequences)")

                if torch.cuda.is_available():
                    toks = toks.to(device="cuda", non_blocking=True)

                with torch.amp.autocast(device_type=device_type, enabled=enable_autocast):
                    out = model(toks, repr_layers=repr_layers, return_contacts=return_contacts)

                representations = {
                    layer: t.to("cpu") for layer, t in out["representations"].items()
                }
                if return_contacts:
                    contacts = out["contacts"].to("cpu")

                for i, label in enumerate(labels):
                    truncate_len = min(truncation_seq_length, len(strs[i]))
                    result = {"label": label}

                    if "per_tok" in include:
                        result["representations"] = {
                            str(layer): t[i, 1 : truncate_len + 1].clone().tolist()
                            for layer, t in representations.items()
                        }

                    if "mean" in include:
                        result["mean_representations"] = {
                            str(layer): t[i, 1 : truncate_len + 1].mean(0).clone().tolist()
                            for layer, t in representations.items()
                        }

                    if "bos" in include:
                        result["bos_representations"] = {
                            str(layer): t[i, 0].clone().tolist()
                            for layer, t in representations.items()
                        }

                    if return_contacts:
                        result["contacts"] = contacts[i, :truncate_len, :truncate_len].clone().tolist()

                    embedding_results.append(result)

        return embedding_results

    async def invoke(self, input: ESM_embeddingInput) -> ESM_embeddingOutput:

        try:
            start_time = time.time()
            
            validate_client_inference_input(input.fasta_base64)
            
            store_input_fasta  = store_client_input(input.fasta_base64,input.port)

            results = self._run_inference(self.model, self.alphabet, store_input_fasta, input.toks_per_batch,input.truncation_seq_length,input.include,input.repr_layers)
            
            print("total run time",time.time() - start_time)
            return ESM_embeddingOutput(
                status="success",
                result=results,
                message="Inference completed successfully."
            )

        except HTTPException as http_exc:
            raise http_exc

    def check_health(self) -> bool:
        """Checks the health of the ESM2 embedding service."""
        if self.model is None:
            logger.error("ESM2 embedding model is not initialized.")
            return False

        return True