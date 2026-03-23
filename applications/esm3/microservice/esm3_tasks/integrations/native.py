# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License
import os
import sys
from pathlib import Path

inversefold_path = Path(__file__).resolve().parents[3] / "scripts"
os.chdir(inversefold_path)
sys.path.insert(0, str(inversefold_path))


import ESM3_logits_embedding_task
import ESM3_folding_task
import ESM3_prompt_sequence
import ESM3_chain_of_thought
import ESM3_function_prediction_task
import ESM3_inversefold_task

from typing import List, Optional
from pydantic import BaseModel
from fastapi import HTTPException

from comps import CustomLogger

import base64
from comps import CustomLogger, OpeaComponent, OpeaComponentRegistry
import os
import time

from esm.sdk import client

from Bio import SeqIO
import numpy as np

from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient

from types import SimpleNamespace


def store_client_input(input_base64,task, server_id: int = 0):
    try:
        work_dir = f"/tmp/esm3_inputs/server_{server_id}"
        os.makedirs(work_dir, exist_ok=True)
        # Support base64 input
        if task=="logits_embedding" or task=="fold" or task=="prompt_task":
            if input_base64:
                decoded = base64.b64decode(input_base64)
                input_file = os.path.join(work_dir, "input.fasta")
                with open(input_file, "wb") as f:
                    f.write(decoded)
                logger.info(f"Decoded base64 FASTA and saved to {input_file}")
            else:
                raise HTTPException(status_code=400, detail="FASTA base64 string is missing.")

        if task=="inversefold" or task=="function_prediction":
            decoded = base64.b64decode(input_base64)
            input_file = os.path.join(work_dir, "input.pdb")
            with open(input_file, "wb") as f:
                f.write(decoded)
            logger.info(f"Decoded base64 PDB and saved to {input_file}")

        if task=="chain_of_thought":
            decoded = base64.b64decode(input_base64)
            input_file = os.path.join(work_dir, "input.csv")
            with open(input_file, "wb") as f:
                f.write(decoded)
            logger.info(f"Decoded base64 PDB and saved to {input_file}")

    except HTTPException as http_exc:
        raise http_exc

    except Exception as e:
        logger.error(f"Error in ESM3 - Task microservice: {e}")
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
        logger.error(f"Error in ESM3 - Task microservice: {e}")
        raise HTTPException(status_code=500, detail=str(e))



WORKLOAD_NAME = "esm3_task"
logger = CustomLogger(f"{WORKLOAD_NAME}_microservice")

class ESM3_taskInput(BaseModel):
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
    schedule: Optional[str] = "cosine"
    strategy: Optional[str] = "entropy"
    num_steps: Optional[int] = 1
    temperature: Optional[float] = -1
    temperature_annealing: Optional[bool] = False
    top_p: Optional[float] = 1.0
    condition_on_coordinates_only: Optional[bool] = False
    task: Optional[str] = None
    num_sequences: Optional[int] = 8
    batch_run: Optional[bool] = False
    sequence_length: Optional[int] = None
    port: Optional[int] = None



class ESM3_taskOutput(BaseModel):
    status: str
    result: Optional[str] = None
    message: Optional[str] = None

@OpeaComponentRegistry.register("OPEA_OMICS_esm3task")
class Opea_ESM3_task(OpeaComponent):

    def __init__(self, name: str, description: str, config: dict = None):
        super().__init__(name, description, config)

        self.use_bf16 = False
        if config.get("bf16"):
            self.use_bf16 = config["bf16"]
            logger.info(f"Info: Running in bfloat16 mode.")

        # print(f"Using model: {model_name_or_path}")
        model_load_time = time.time()



        if os.environ.get("ESM_API_KEY", ""):
            print("ESM_API_KEY found. Using model from Forge API...")
            self.client = ESM3InferenceClient()
        else:
            print("No ESM_API_KEY found. Loading model locally...")
            self.client = ESM3.from_pretrained("esm3_sm_open_v1", bf16=self.use_bf16)

        print("Model load time",time.time() - model_load_time)



    def _run_inference(
        self,
        client: ESM3InferenceClient,
        input_file,
        sequence,
        structure,
        secondary_structure,
        sasa,
        function,
        residue_annotations,
        return_embeddings,
        ith_hidden_layer,
        return_hidden_states,
        protein_complex,bf16,
        schedule,strategy,num_steps,temperature,temperature_annealing,top_p,condition_on_coordinates_only,task,
        num_sequences,batch_run,sequence_length


    ):


        # Wrap all individual parameters into an args object
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
            bf16=bf16,
            schedule=schedule,
            strategy=strategy,
            num_steps=num_steps,
            temperature=temperature,
            temperature_annealing=temperature_annealing,
            top_p=top_p,
            condition_on_coordinates_only=condition_on_coordinates_only,
            num_sequences=num_sequences,
            batch_run=batch_run,
            sequence_length=sequence_length
        )

        # output_dir = None

        if task == "logits_embedding":
            result = ESM3_logits_embedding_task.processing_fasta(client, input_file,  args)
        elif  task == "fold":
            result = ESM3_folding_task.processing_fasta(client, input_file,  args)
        elif task =="inversefold":
            result=ESM3_inversefold_task.processing_pdb(client, input_file,  args)
        elif task=="chain_of_thought":
            result=ESM3_chain_of_thought.chain_of_thought(client, input_file, args)
        elif task=="function_prediction":
            result=ESM3_function_prediction_task.processing_pdb(client, input_file,  args)
        elif task=="prompt_task":
            result=ESM3_prompt_sequence.processing_fasta(client, input_file,  args)

        return result

    async def invoke(self, input: ESM3_taskInput) -> ESM3_taskOutput:

        try:
            validate_client_inference_input(input.fasta_base64)

            store_input_fasta  = store_client_input(input.fasta_base64,input.task,input.port)
            results = self._run_inference( self.client,store_input_fasta,input.sequence,input.structure,input.secondary_structure,input.sasa,input.function,input.residue_annotations,input.return_embeddings,input.ith_hidden_layer,input.return_hidden_states,input.protein_complex,self.use_bf16,
                                          input.schedule,input.strategy,input.num_steps,input.temperature,input.temperature_annealing,input.top_p,input.condition_on_coordinates_only,input.task,input.num_sequences,input.batch_run,
                                          input.sequence_length)

            import json, base64

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

            return ESM3_taskOutput(
                status="success",
                result=encoded_result,
                message="EMS3 generated successfully"
            )
        except HTTPException as http_exc:
            raise http_exc

    def check_health(self) -> bool:
        """Checks the health of the ESM3 Task service."""
        if self.model is None:
            logger.error("ESM3 task model is not initialized.")
            return False

        return True
