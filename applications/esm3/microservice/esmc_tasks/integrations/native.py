# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License
import os
import sys
from pathlib import Path

inversefold_path = Path(__file__).resolve().parents[3]  / "scripts"
os.chdir(inversefold_path)
sys.path.insert(0, str(inversefold_path))

import ESMC_logits_embedding_task

import time
from typing import Optional
from pydantic import BaseModel
from fastapi import HTTPException

from comps import CustomLogger

import base64
from comps import CustomLogger, OpeaComponent, OpeaComponentRegistry
from esm.models.esmc import ESMC
from esm.sdk.api import (
    ESMCInferenceClient,
)

from types import SimpleNamespace
import json


def store_client_input(fasta_base64, server_id: int = 0):
    try:
        work_dir = f"/tmp/esmc_inputs/server_{server_id}"
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
        logger.error(f"Error in ESMC - Embedding microservice: {e}")
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
        logger.error(f"Error in ESMC - Embedding microservice: {e}")
        raise HTTPException(status_code=500, detail=str(e))



WORKLOAD_NAME = "esmc_embbeding"
logger = CustomLogger(f"{WORKLOAD_NAME}_microservice")


class ESMC_embeddingInput(BaseModel):
    fasta_base64: str
    sequence: Optional[bool] = False
    structure: Optional[bool] = False
    secondary_structure: Optional[bool] = False
    sasa: Optional[bool] = False
    function: Optional[bool] = False

    residue_annotations: Optional[bool] = False
    return_embeddings: Optional[bool] = False
    return_hidden_states: Optional[bool] = False
    ith_hidden_layer: Optional[float] = -1
    protein_complex: Optional[bool] = False
    port: Optional[int] = None

class ESMC_embeddingOutput(BaseModel):
    status: str
    result: Optional[str] = None
    message: Optional[str] = None

@OpeaComponentRegistry.register("OPEA_OMICS_esmcembedding")
class Opea_ESMC_embedding(OpeaComponent):

    def __init__(self, name: str, description: str, config: dict = None):
        super().__init__(name, description, config)

        model_name_or_path = os.getenv("MODEL") or config.get("model_name") if config else None

        if not model_name_or_path:
            model_name_or_path = "esmc_600m"

        self.use_bf16 = False
        if config.get("bf16"):
            self.use_bf16 = config["bf16"]
            logger.info(f"Info: Running in bfloat16 mode.")

        print(f"Using model: {model_name_or_path}")
        model_load_time = time.time()

        self.model = ESMC.from_pretrained(model_name_or_path, bf16=self.use_bf16)
        print("Model load time",time.time() - model_load_time)





    def _run_inference(self,  client: ESMCInferenceClient,fasta_file,sequence,structure,secondary_structure,sasa,function, residue_annotations,return_embeddings,ith_hidden_layer,
                       return_hidden_states,protein_complex,bf16):
        args = SimpleNamespace(
            sequence=sequence,
            structure=structure,
            secondary_structure=secondary_structure,
            sasa=sasa,
            function=function,
            residue_annotations=residue_annotations,
            return_embeddings=return_embeddings,
            ith_hidden_layer=ith_hidden_layer,
            return_hidden_states=return_hidden_states,
            protein_complex=protein_complex,
            bf16=bf16
        )

        result = ESMC_logits_embedding_task.processing_fasta(client, fasta_file, args)
        return result

    async def invoke(self, input: ESMC_embeddingInput) -> ESMC_embeddingOutput:

        try:
            start_time = time.time()

            validate_client_inference_input(input.fasta_base64)

            store_input_fasta  = store_client_input(input.fasta_base64,input.port)

            results = self._run_inference( self.model,store_input_fasta,input.sequence,input.structure,input.secondary_structure,input.sasa,input.function,input.residue_annotations,input.return_embeddings,input.ith_hidden_layer,input.return_hidden_states,input.protein_complex,bf16=self.use_bf16)

            def make_serializable(obj):
                if isinstance(obj, list):
                    return [make_serializable(x) for x in obj]
                elif hasattr(obj, "model_dump"):  # Pydantic v2 (used in ESM3)
                    return obj.model_dump()
                elif hasattr(obj, "__dict__"):
                    return {k: make_serializable(v) for k, v in obj.__dict__.items()}
                elif isinstance(obj, (dict, tuple, set)):
                    return json.loads(json.dumps(obj, default=str))
                else:
                    return obj

            serializable_data = make_serializable(results)
            encoded_result = base64.b64encode(json.dumps(serializable_data).encode("utf-8")).decode("utf-8")

            return ESMC_embeddingOutput(
                status="success",
                result=encoded_result,
                message="Embedding generated successfully"
            )
        except HTTPException as http_exc:
            raise http_exc

    def check_health(self) -> bool:
        """Checks the health of the ESMC embedding service."""
        if self.model is None:
            logger.error("ESMC embedding model is not initialized.")
            return False

        return True
