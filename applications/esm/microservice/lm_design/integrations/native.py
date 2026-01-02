# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import os
import sys
from pathlib import Path

lmdesign_path = Path(__file__).resolve().parents[3] / "examples" / "lm-design"
os.chdir(lmdesign_path)
sys.path.insert(0, str(lmdesign_path))

from pathlib import Path

import lm_design

import time
from typing import Optional, List
from typing import Dict, Any, Optional
import socket

from comps import (
    CustomLogger,
    ServiceType,
    opea_microservices,
    register_microservice,
    register_statistics,
    statistics_dict,
)
logger = CustomLogger("lmdesign_microservice")
import re
import time, pickle
import torch
from omegaconf import OmegaConf
import hydra
import logging
from hydra.core.hydra_config import HydraConfig
import numpy as np
import random
import glob
from omegaconf import DictConfig
from hydra import initialize_config_dir, compose
from hydra.core.global_hydra import GlobalHydra
import argparse
from fastapi import HTTPException
import base64
import json

from comps import CustomLogger, OpeaComponent, OpeaComponentRegistry
from pydantic import BaseModel

from lm_design import Designer



logger = CustomLogger("opea")



def store_client_input(cfg_dict, server_id: int = 0):
    try:
        work_dir = f"/tmp/lm_design/server_{server_id}"
        os.makedirs(work_dir, exist_ok=True)
        
        pdb_data = cfg_dict.get("pdb_fn")
        if not pdb_data:
            logger.info("No pdb_fn provided — skipping file write.")
            return cfg_dict  # nothing to save

        # Detect if pdb_fn is a PDB text string or path
        if isinstance(pdb_data, str) and pdb_data.strip().startswith("ATOM"):
            # Save as file
            output_path = os.path.join(work_dir, "input.pdb")
            with open(output_path, "w") as f:
                f.write(pdb_data.strip() + "\n")
            cfg_dict["pdb_fn"] = output_path
            logger.info(f"Saved inline PDB string to {output_path}")

        elif isinstance(pdb_data, str):
            # Assume it's a file path — verify it exists
            if not os.path.exists(pdb_data):
                raise HTTPException(
                    status_code=400,
                    detail=f"PDB file not found at: {pdb_data}"
                )
            logger.info(f"Using existing PDB file: {pdb_data}")

        elif isinstance(pdb_data, list) and len(pdb_data) > 0:
            # Handle case of list of paths
            cfg_dict["pdb_fn"] = pdb_data

        else:
            raise HTTPException(
                status_code=400,
                detail="Invalid 'pdb_fn' format: must be a path, list of paths, or inline PDB string."
            )

        return cfg_dict

    except HTTPException as http_exc:
        raise http_exc

    except Exception as e:
        logger.error(f"Error in RFdiffusion microservice: {e}")
        raise HTTPException(status_code=500, detail=str(e))

def validate_client_inference_input(cfg_dict):
    try:
        if cfg_dict is None:
            logger.error("Received null cfg_dict in request.")
            raise HTTPException(
                status_code=400,
                detail="Missing required input: cfg_dict"
            )

        pdb_fn = cfg_dict.get("pdb_fn")
        task = cfg_dict.get("task")

        if task not in ["fixedbb", "free_generation"]:
            raise HTTPException(
                status_code=400,
                detail="Missing or invalid 'task'. Must be either 'fixedbb' or 'free_generation'."
            )

        if task == "fixedbb":
            if not pdb_fn or not str(pdb_fn).strip():
                raise HTTPException(
                    status_code=400,
                    detail="'pdb_fn' is required for task='fixedbb'."
                )

        elif task == "free_generation":
            if not pdb_fn:
                logger.info("No pdb_fn provided for free_generation — continuing without pdb_fn.")

    except HTTPException as http_exc:
        raise http_exc

    except Exception as e:
        logger.error(f"Error in RFdiffusion microservice: {e}")
        raise HTTPException(status_code=500, detail=str(e))
    
def encode_base64_json(obj):
    """
    Encode a JSON-serializable Python object into base64 string.
    Useful for storing complex data like metadata.
    """
    try:
        json_bytes = json.dumps(obj, default=str).encode("utf-8")
        return base64.b64encode(json_bytes).decode("utf-8")
    except Exception as e:
        raise ValueError(f"Failed to base64 encode: {e}")


def get_machine_ip():
    try:
        s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        s.connect(("8.8.8.8", 80))
        ip = s.getsockname()[0]
        s.close()
        return ip
    except Exception:
        return "127.0.0.1"

WORKLOAD_NAME = "esm_lmdesign"
logger = CustomLogger(f"{WORKLOAD_NAME}_microservice")

class ESM_lmdesignInput(BaseModel):
    cfg_dict: Optional[Dict[str, Any]] = None
    port: Optional[int] = None

class ESM_lmdesignOutput(BaseModel):
    status: str
    results: Optional[str] = None
    message: Optional[str] = None
    
        
@OpeaComponentRegistry.register("OPEA_OMICS_esmlmdesign")
class Opea_ESM_lmdesign(OpeaComponent):

    def __init__(self, name: str, description: str, config: dict = None):
        super().__init__(name, description, config)
        self.use_bf16 = False
        if config.get("bf16"):
            self.use_bf16 = config["bf16"]
            logger.info(f"Info: Running in bfloat16 mode.")

        model_load_time = time.time()
        self.designer_model_only = self._load_model()
        print("Model load time",time.time() - model_load_time)
        
    def _load_model(self):
        
        cfg = OmegaConf.load("conf/config.yaml")
        designer_model_only = Designer(cfg=cfg)
        return designer_model_only
    
    def _run_infernce(self,cfg):
        
        des = Designer.__new__(Designer)
        des.cfg = cfg
        des.vocab = self.designer_model_only .vocab
        des.vocab_mask_AA = self.designer_model_only .vocab_mask_AA
        des.vocab_mask_AA_idx = self.designer_model_only .vocab_mask_AA_idx
        des.struct_model = self.designer_model_only .struct_model
        des.LM = self.designer_model_only .LM
        des.pdb_loader_params = self.designer_model_only .pdb_loader_params
        des.device = self.designer_model_only .device
        des.allowed_AA = self.designer_model_only .allowed_AA
        
        if self.use_bf16:
            cfg["bf16"] = True
            print("convert to bf16 mode")
        des = Designer(cfg)
            
        if cfg.task == "fixedbb":
            des._init_target(Path(cfg["pdb_fn"]))
        elif cfg.task == "free_generation":
            if not hasattr(des, "L"):
                des.L = cfg.get("length", 100)

        des.init_sequences(cfg.num_seqs)
        start_time = time.time()
        result = des.run_from_cfg()
        elapsed = time.time() - start_time
        return result



    async def invoke(self, input: ESM_lmdesignInput) -> ESM_lmdesignOutput:

        try:
            total_time = time.time()
            cfg_dict = input.cfg_dict
            
            validate_client_inference_input(cfg_dict)
            
            store_input_cfg_dict = store_client_input(cfg_dict,input.port)
            
            pdb_fn = cfg_dict["pdb_fn"]
            cfg = OmegaConf.create(cfg_dict)
            infer_time = time.time()
            result= self._run_infernce(cfg)
            print(f"infer_time = {time.time() - infer_time:.2f} seconds")
            logger.info(f"Inference time: {time.time() - infer_time:.2f} sec")
            result_str = encode_base64_json(result)
            if not result:
                raise HTTPException(status_code=500, detail="Inference completed but returned no output")

            logger.info(f"Total request time: {time.time() - total_time:.2f} sec")
            print(f"total run time = {time.time() - total_time:.2f} seconds")

            return ESM_lmdesignOutput(
                    status="success",
                    results=result_str,
                    message="LMDesign inference complete."
                )
        except HTTPException as http_exc:
            raise http_exc


    def check_health(self) -> bool:
        """Checks the health of the ESM2 InverseFold Sample service."""
        if self.designer_model_only is None:
            logger.error("ESM2 InverseFold Sample Model is not initialized.")
            return False

        return True