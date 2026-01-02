# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import argparse
import requests
import time
from pathlib import Path
import hydra
from omegaconf import DictConfig, OmegaConf
import sys
import pickle
import base64
import json



# Argument parser setup
parser = argparse.ArgumentParser()
parser.add_argument("--host", default="0.0.0.0", help="Microservice host")
parser.add_argument("--port", default="9000", help="Microservice port")
args, unknown = parser.parse_known_args()
sys.argv = [sys.argv[0]] + unknown  # Clean up for Hydra


def read_pdb_file(pdb_path: Path) -> str:
    """Extract ATOM and HETATM lines from a PDB file."""
    if not pdb_path.exists():
        raise FileNotFoundError(f"File not found: {pdb_path}")
    return "\n".join(
        line for line in pdb_path.read_text().splitlines()
        if line.startswith(("ATOM", "HETATM"))
    )

#Replaces the path to the PDB file stored in cfg[key] with the contents of the file as a string.
def replace_pdb_path_with_string(cfg: DictConfig, key: str) -> None:
    """Replace a PDB path in the config with its reduced string content."""
    try:
        pdb_file_path = OmegaConf.select(cfg, key)
        if pdb_file_path:
            pdb_str = read_pdb_file(Path(pdb_file_path))
            OmegaConf.update(cfg, key, pdb_str)
    except Exception as e:
        print(f"Failed to process PDB at {key}: {e}")

@hydra.main(version_base=None, config_path="./config/", config_name="base")
def main(cfg: DictConfig) -> None:
    # Prepare the request URL
    url = f"http://{args.host}:{args.port}/v1/rfdiffusion"

    # Reduce PDB paths to strings
    replace_pdb_path_with_string(cfg, "inference.input_pdb")
    replace_pdb_path_with_string(cfg, "scaffoldguided.target_path")

    # Convert config to dict for payload
    payload = {"cfg_dict": OmegaConf.to_container(cfg, resolve=True),"port":args.port}

    print("Sending POST to:", url)
    try:
        response = requests.post(url, json=payload)
        print(f"Status: {response.status_code}")

        if response.status_code == 200:
            data = response.json()
            # print("data",data["out_dict"]["write_pdb_string"])
            # Final design PDB

            if data["out_dict"].get("write_pdb_string"):
                Path("output_design.pdb").write_text(data["out_dict"]["write_pdb_string"])
                print("Design PDB written to output_design.pdb")
            else:
                print("No 'write_pdb_string' in response")

            if data["out_dict"].get("xt1_traj_pdb_string"):
                Path("output_Xt-1_traj.pdb").write_text(data["out_dict"]["xt1_traj_pdb_string"])
                print("Xt-1 Trajectory PDB written to output_Xt-1_traj.pdb")
            else:
                print("No 'xt1_traj_pdb_string' in response")

            if data["out_dict"].get("px0_traj_pdb_string"):
                Path("output_pX0_traj.pdb").write_text(data["out_dict"]["px0_traj_pdb_string"])
                print("pX0 Trajectory PDB written to output_pX0_traj.pdb")
            else:
                print("No 'px0_traj_pdb_string' in response")

            if data["out_dict"].get("metadata"):
                decoded_json_str = base64.b64decode(data["out_dict"]["metadata"]).decode("utf-8")
                metadata_obj = json.loads(decoded_json_str)

                with open("output.trb", "wb") as meta_file:
                    pickle.dump(metadata_obj, meta_file)

                print("Metadata written to output.trb")
            else:
                print("No 'metadata' in response")
        else:
            print(f"Server returned error:\n{response.text}")

    except Exception as e:
        print(f"Request failed: {e}")


if __name__ == "__main__":
    main()
