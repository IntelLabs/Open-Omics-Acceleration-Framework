# Copyright (C) 2024 Intel Corporation
# SPDX-License-Identifier: Apache-2.0

import time
import socket
import torch
import netifaces
from fastapi import HTTPException
from typing import List, Optional
from pydantic import BaseModel
from transformers import pipeline
from comps import CustomLogger, OpeaComponent, OpeaComponentRegistry


try:
    import intel_extension_for_pytorch as ipex
except ImportError:
    ipex = None

logger = CustomLogger("ProtGPT2")

def validate_client_inference_input(max_length, num_return_sequences):
    try:
        if max_length is None or num_return_sequences is None:
            logger.error("Received null max_length or num_return_sequences in request.")
            raise HTTPException(
                status_code=400,
                detail="Missing required input: max_length and num_return_sequences"
            )
        return True  # valid case

    except HTTPException as http_exc:
        raise http_exc

    except Exception as e:
        logger.error(f"Error in ProtGPT2 microservice: {e}")
        raise HTTPException(status_code=500, detail=str(e))

# -------------------------
# Input / Output schemas
# -------------------------
class ProtGPT2Input(BaseModel):
    max_length: int
    num_return_sequences: int
    do_sample: Optional[bool] = True
    top_k: Optional[int] = 950
    repetition_penalty: Optional[float] = 1.2
    eos_token_id: Optional[int] = 0
    dtype: Optional[str] = False
    model_dir: Optional[str] = None     # path to custom model
    iterations: Optional[int] = 1
    seed: Optional[int] = None
class ProtGPT2Output(BaseModel):
    sequences: List[str]


def load_model(dtype_str: str, model_dir: Optional[str] = None):
    if dtype_str == "bf16":
        dtype = torch.bfloat16
    else:
        dtype_str = "fp32"
        dtype = torch.float32

    model_name = model_dir if model_dir else "nferruz/ProtGPT2"
    logger.info(f"Loading ProtGPT2 model: {model_name} with dtype={dtype_str}")

    model = pipeline("text-generation", model=model_name, torch_dtype=dtype)

    if ipex:
        logger.info("Applying IPEX optimization")
        model.model = ipex.optimize(model.model, dtype=dtype)

    return model


# -------------------------
# OPEA ProtGPT2 Component
# -------------------------
@OpeaComponentRegistry.register("OPEA_OMICS_PROTGPT2")
class Opea_ProtGPT2(OpeaComponent):
    """Wrapper for ProtGPT2 model used in OPEA microservice."""

    def __init__(self, name: str, description: str, config: dict = None):
        super().__init__(name, description, config)

        logger.info("Loading model")

        # default bfloat16
        self.use_bf16 = True
        if config and config.get("fp32"):   # expects True/False flag
            self.use_bf16 = False
            logger.info("Info: Running in float32 mode.")
        else:
            logger.info("Info: Running in bfloat16 mode.")

        model_dir = config.get("model_dir") if config else None
        model_load_time = time.time()

        # pass dtype string to loader
        dtype_str = "bf16" if self.use_bf16 else "fp32"
        logger.info(f"Info: Running in {dtype_str} mode.")
        self.model = load_model(dtype_str, model_dir)

        logger.info(f"Model loaded in {time.time() - model_load_time:.2f} seconds")

    async def invoke(self, input: ProtGPT2Input) -> ProtGPT2Output:
        try:
            # Validate mandatory fields
            validate_client_inference_input(input.max_length, input.num_return_sequences)

            tic = time.time()
            if input.seed is not None:
                torch.manual_seed(input.seed)
                logger.info(f"Using seed {input.seed}")
            for i in range(input.iterations):
                logger.info(f"Iteration {i+1}/{input.iterations}")
                t0 = time.time()
                results = self.model(
                    "<|endoftext|>",
                    max_length=input.max_length,
                    do_sample=input.do_sample,
                    top_k=input.top_k,
                    repetition_penalty=input.repetition_penalty,
                    num_return_sequences=input.num_return_sequences,
                    eos_token_id=input.eos_token_id
                )
                t1 = time.time()
                print('Time taken for', i, 'iteration:', t1 - t0, 'seconds')

            seqs = [r["generated_text"] for r in results]
            toc = time.time()
            print(f"Total time taken: {toc - tic:.2f} seconds")
            print(f"Average time per iteration: {(toc - tic) / input.iterations:.2f} seconds")
            logger.info(f"Request served in {toc - tic:.2f} seconds")

            return ProtGPT2Output(
                status="success",
                sequences=seqs,
                message="Inference completed successfully."
            )

        except HTTPException as http_exc:
            raise http_exc
        except Exception as e:
            logger.error(f"Unexpected error in ProtGPT2 inference: {e}")
            raise HTTPException(status_code=500, detail=str(e))

    def check_health(self) -> bool:
        if self.model is None:
            logger.error("ProtGPT2 model is not initialized.")
            return False
        return True
