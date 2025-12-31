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
parser.add_argument("--port", default="9005", help="Microservice port")
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
    print(" is enter into fuction replace_pdb_path_with_string")
    """Replace a PDB path in the config with its reduced string content."""
    try:
        # print("cfg")
        # print(key)
        pdb_file_path = OmegaConf.select(cfg, key)
        # print("pdb_file_path",pdb_file_path)
        if pdb_file_path:
            pdb_str = read_pdb_file(Path(pdb_file_path))
            OmegaConf.update(cfg, key, pdb_str)
        
    except Exception as e:
        print(f"Failed to process PDB at {key}: {e}")

@hydra.main(config_path="../../examples/lm-design/conf", config_name="config")
def main(cfg: DictConfig) -> None:
    
    # Prepare the request URL
    url = f"http://{args.host}:{args.port}/v1/esmlmdesign"

    args_no_spaces = [arg.replace(" ", "") for arg in sys.argv[1:]]    
    pdb_fn = cfg.pdb_fn

    if pdb_fn:
        replace_pdb_path_with_string(cfg, "pdb_fn")

    payload = {"cfg_dict": OmegaConf.to_container(cfg, resolve=True), "port": args.port}

    print("Sending POST to:", url)
    try:
        response = requests.post(url, json=payload)
        print(f"Status: {response.status_code}")

        if response.status_code == 200:
            data = response.json()
            decoded_json_str = base64.b64decode(data["results"]).decode("utf-8")
            
            # Save decoded FASTA output to file
            fasta_path = Path("output.fasta")

            with open(fasta_path, "w") as f:
                f.write(decoded_json_str.strip() + "\n")

            print(f"FASTA file saved to: {fasta_path.resolve()}")

        else:
            print(f"Server returned error:\n{response.text}")

    except Exception as e:
        print(f"Request failed: {e}")


if __name__ == "__main__":
    main()