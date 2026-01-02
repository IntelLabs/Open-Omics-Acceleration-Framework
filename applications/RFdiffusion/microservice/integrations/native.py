# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import torch
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import time
from typing import Dict, Any, Optional

import re
from omegaconf import OmegaConf
import logging
from rfdiffusion.util import writepdb_str, writepdb_multi_str
import numpy as np
import random
import glob
from omegaconf import DictConfig
from rfdiffusion.inference import model_runners
from fastapi import HTTPException
import base64
import json

from comps import CustomLogger, OpeaComponent, OpeaComponentRegistry
from pydantic import BaseModel


logger = CustomLogger("rfdiffusion_microservice")

class RFDiffusionInput(BaseModel):
    cfg_dict: Optional[Dict[str, Any]] = None
    port: Optional[int] = None

class RFDiffusionOutput(BaseModel):
    status: str
    out_dict: Optional[Dict[str, Any]] = None
    message: Optional[str] = None

def make_deterministic(seed=0):
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)

def sampler_selector(conf: DictConfig, samplers: dict):
    torch.cuda.is_available = lambda : False
    if conf.inference.deterministic:
        make_deterministic()
    if conf.scaffoldguided.scaffoldguided:
        if samplers["scaffolded"] is None:
            samplers["scaffolded"] = model_runners.ScaffoldedSampler(conf)
        else:
            samplers["scaffolded"].initialize(conf)
        return samplers["scaffolded"]

    elif conf.inference.model_runner == "default":
        if samplers["default"] is None:
            samplers["default"] = model_runners.Sampler(conf)
        else:
            samplers["default"].initialize(conf)
        return samplers["default"]

    elif conf.inference.model_runner == "SelfConditioning":
        if samplers["self_conditioning"] is None:
            samplers["self_conditioning"] = model_runners.SelfConditioning(conf)
        else:
            samplers["self_conditioning"].initialize(conf)
        return samplers["self_conditioning"]

    elif conf.inference.model_runner == "ScaffoldedSampler":
        if samplers["scaffolded"] is None:
            samplers["scaffolded"] = model_runners.ScaffoldedSampler(conf)
        else:
            samplers["scaffolded"].initialize(conf)
        return samplers["scaffolded"]

    else:
        raise HTTPException(status_code=400, detail=f"Unrecognized sampler {conf.inference.model_runner}")

def run_inference(sampler):

    log = logging.getLogger(__name__)
    # Loop over number of designs to sample.
    design_startnum = sampler.inf_conf.design_startnum
    if sampler.inf_conf.design_startnum == -1:
        existing = glob.glob(sampler.inf_conf.output_prefix + "*.pdb")
        indices = [-1]
        for e in existing:
            print(e)
            m = re.match(".*_(\d+)\.pdb$", e)
            print(m)
            if not m:
                continue
            m = m.groups()[0]
            indices.append(int(m))
        design_startnum = max(indices) + 1

    for i_des in range(design_startnum, design_startnum + sampler.inf_conf.num_designs):
        if sampler.inf_conf.deterministic:
            make_deterministic(i_des)

        start_time = time.time()

        x_init, seq_init = sampler.sample_init()
        denoised_xyz_stack = []
        px0_xyz_stack = []
        seq_stack = []
        plddt_stack = []

        if sampler.inf_conf.precision == "bfloat16":
            dtype=torch.bfloat16
            x_t = torch.clone(x_init)
            seq_t = torch.clone(seq_init)
            x_t = x_t.to(dtype=dtype)
            seq_t = seq_t.to(dtype=dtype)
        else:
            x_t = torch.clone(x_init)
            seq_t = torch.clone(seq_init)

        # Loop over number of reverse diffusion time steps.
        for t in range(int(sampler.t_step_input), sampler.inf_conf.final_step - 1, -1):
            px0, x_t, seq_t, plddt = sampler.sample_step(
                t=t, x_t=x_t, seq_init=seq_t, final_step=sampler.inf_conf.final_step
            )
            px0_xyz_stack.append(px0)
            denoised_xyz_stack.append(x_t)
            seq_stack.append(seq_t)
            plddt_stack.append(plddt[0])  # remove singleton leading dimension

        # Flip order for better visualization in pymol
        denoised_xyz_stack = torch.stack(denoised_xyz_stack)
        denoised_xyz_stack = torch.flip(
            denoised_xyz_stack,
            [
                0,
            ],
        )
        px0_xyz_stack = torch.stack(px0_xyz_stack)
        px0_xyz_stack = torch.flip(
            px0_xyz_stack,
            [
                0,
            ],
        )

        # For logging -- don't flip
        plddt_stack = torch.stack(plddt_stack)

        final_seq = seq_stack[-1]
        final_seq = torch.where(
            torch.argmax(seq_init, dim=-1) == 21, 7, torch.argmax(seq_init, dim=-1)
        )  # 7 is glycine

        bfacts = torch.ones_like(final_seq.squeeze())
        # make bfact=0 for diffused coordinates
        bfacts[torch.where(torch.argmax(seq_init, dim=-1) == 21, True, False)] = 0

        write_pdb_string = writepdb_str(
            denoised_xyz_stack[0, :, :4],
            final_seq,
            sampler.binderlen,
            chain_idx=sampler.chain_idx,
            bfacts=bfacts,
        )

        # run metadata
        trb = dict(
            config=OmegaConf.to_container(sampler._conf, resolve=True),
            plddt=plddt_stack.cpu().numpy(),
            device='cpu',
            time=time.time() - start_time,
        )
        if hasattr(sampler, "contig_map"):
            for key, value in sampler.contig_map.get_mappings().items():
                trb[key] = value

        xt1_pdb_str = None
        px0_pdb_str = None
        if sampler.inf_conf.write_trajectory:
            xt1_pdb_str = writepdb_multi_str(
                denoised_xyz_stack,
                bfacts,
                final_seq.squeeze(),
                use_hydrogens=False,
                backbone_only=False,
                chain_ids=sampler.chain_idx,
            )

            px0_pdb_str = writepdb_multi_str(
                px0_xyz_stack,
                bfacts,
                final_seq.squeeze(),
                use_hydrogens=False,
                backbone_only=False,
                chain_ids=sampler.chain_idx,
            )

        log.info(f"Design finished in {(time.time() - start_time) / 60:.2f} minutes")

        return {
            "write_pdb_string": write_pdb_string,
            "xt1_traj_pdb_string": xt1_pdb_str,
            "px0_traj_pdb_string": px0_pdb_str,
            "metadata": trb,
        }

def store_client_input(cfg_dict, server_id: int = 0):
    try:
        # Each server gets its own folder: /tmp/rfdiffusion_inputs/server_{id}
        work_dir = f"/tmp/rfdiffusion_inputs/server_{server_id}"
        os.makedirs(work_dir, exist_ok=True)

        # Handle inference.input_pdb
        if "input_pdb" in cfg_dict.get("inference", {}) and cfg_dict["inference"]["input_pdb"]:
            input_pdb_str = cfg_dict["inference"]["input_pdb"]
            if not input_pdb_str.strip():
                raise HTTPException(status_code=400, detail="inference.input_pdb is empty")

            input_pdb_path = os.path.join(work_dir, "input.pdb")
            with open(input_pdb_path, "w") as f:
                f.write(input_pdb_str + "\n")

            cfg_dict["inference"]["input_pdb"] = input_pdb_path
            logger.info(f"[Server {server_id}] Saved input_pdb to {input_pdb_path}")
        else:
            logger.info(f"[Server {server_id}] No inference.input_pdb provided — skipping file write.")

        # Handle scaffoldguided.target_path
        if "scaffoldguided" in cfg_dict and "target_path" in cfg_dict["scaffoldguided"] and cfg_dict["scaffoldguided"]["target_path"]:
            target_pdb_str = cfg_dict["scaffoldguided"]["target_path"]
            if not target_pdb_str.strip():
                raise HTTPException(status_code=400, detail="scaffoldguided.target_path is empty")

            target_pdb_path = os.path.join(work_dir, "target.pdb")
            with open(target_pdb_path, "w") as f:
                f.write(target_pdb_str + "\n")

            cfg_dict["scaffoldguided"]["target_path"] = target_pdb_path
            logger.info(f"[Server {server_id}] Saved scaffoldguided.target_path to {target_pdb_path}")
        else:
            logger.info(f"[Server {server_id}] No scaffoldguided.target_path provided — skipping file write.")

        return cfg_dict  # Return updated config

    except HTTPException as http_exc:
        raise http_exc



def validate_client_inference_input(cfg_dict):
    try:

        if cfg_dict is None:
            logger.error("Received null cfg_dict in request.")
            raise HTTPException(
                status_code=400,
                detail="Missing required input: cfg_dict"
            )

        if "inference" not in cfg_dict or not isinstance(cfg_dict["inference"], dict):
            raise HTTPException(status_code=400, detail="Missing or invalid 'inference' block in cfg_dict")

        if "model_runner" not in cfg_dict["inference"]:
            raise HTTPException(status_code=400, detail="Missing inference.model_runner in cfg_dict")

        # Validate model_runner type
        valid_runners = {"default", "SelfConditioning", "ScaffoldedSampler"}
        if cfg_dict["inference"]["model_runner"] not in valid_runners:
            raise HTTPException(status_code=400, detail=f"Invalid model_runner: {cfg_dict['inference']['model_runner']}")

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

samplers = {
        "default": None,
        "scaffolded": None,
        "self_conditioning": None,
        }

@OpeaComponentRegistry.register("OPEA_OMICS_RFDIFFUSION")
class Opea_RFdiffusion(OpeaComponent):

    def __init__(self, name: str, description: str, config: dict = None):
        super().__init__(name, description, config)
        self.precision = "float32"
        if config.get("bfloat16"):
            self.precision = "bfloat16"
            logger.info(f"Info: Running in bfloat16 mode.")

    async def invoke(self, input: RFDiffusionInput) -> RFDiffusionOutput:

        try:
            total_time = time.time()
            cfg_dict = input.cfg_dict
            validate_client_inference_input(cfg_dict)

            store_input_cfg_dict = store_client_input(cfg_dict,input.port)
            if self.precision != cfg_dict["inference"]["precision"]:
                cfg_dict["inference"]["precision"] = self.precision
                logger.info(f"Warning: Running with precision: {self.precision}")
            try:
                conf = OmegaConf.create(store_input_cfg_dict)
            except Exception as e:
                raise HTTPException(status_code=400, detail=f"Invalid cfg_dict format: {e}")

            load_time = time.time()
            currSampler = sampler_selector(conf, samplers)
            print(f"model_load_time = {time.time() - load_time:.2f} seconds")
            logger.info(f"Model loaded in {time.time() - load_time:.2f} sec")

            infer_time = time.time()
            result = run_inference(currSampler)
            print(f"infer_time = {time.time() - infer_time:.2f} seconds")
            logger.info(f"Inference time: {time.time() - infer_time:.2f} sec")

            if not result:
                raise HTTPException(status_code=500, detail="Inference completed but returned no output")

            logger.info(f"Total request time: {time.time() - total_time:.2f} sec")
            print(f"total run time = {time.time() - total_time:.2f} seconds")

            return RFDiffusionOutput(
                    status="success",
                    out_dict={
                        "write_pdb_string": result["write_pdb_string"],
                        "xt1_traj_pdb_string": result["xt1_traj_pdb_string"],
                        "px0_traj_pdb_string": result["px0_traj_pdb_string"],
                        "metadata": encode_base64_json(result["metadata"])
                    },
                    message="RFdiffusion inference complete."
                )
        except HTTPException as http_exc:
            raise http_exc


    def check_health(self) -> bool:
        return True
