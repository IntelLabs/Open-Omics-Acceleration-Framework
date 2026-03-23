# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import argparse
import requests
import json
import pathlib
from pathlib import Path
import torch
import base64

def main():
    parser = argparse.ArgumentParser(description="Client for ESM Embed microservice")
    parser.add_argument("--host", default="0.0.0.0", help="Microservice host")
    parser.add_argument("--port", default="9002", help="Microservice port")
    parser.add_argument("--fasta_file", required=True, help="Path to input FASTA file")
    parser.add_argument(
        "--num_recycles",
        type=int,
        default=None,
        help="Number of recycles to run. Defaults to number used in training (4).",
    )
    parser.add_argument(
        "--max_tokens_per_batch",
        type=int,
        default=1024,
        help="Maximum number of tokens per gpu forward-pass. This will group shorter sequences together "
        "for batched prediction. Lowering this can help with out of memory issues, if these occur on "
        "short sequences.",
    )
    
    args = parser.parse_args()

    # Construct URL
    url = f"http://{args.host}:{args.port}/v1/esmfold"

    # Build payload with base64-encoded FASTA
    with open(args.fasta_file, "rb") as f:
        fasta_b64 = base64.b64encode(f.read()).decode("utf-8")

    payload = {
        "fasta_base64": fasta_b64,
        "num_recycles": args.num_recycles,
        "max_tokens_per_batch": args.max_tokens_per_batch,
        "port": args.port
    }
    print("Sending POST to:", url)

    try:
        response = requests.post(url, json=payload)
        print(f"Status Code: {response.status_code}")

        if response.status_code == 200:
            data = response.json()
            for label, pdb_string in data["result"].items():
                if isinstance(pdb_string, dict):
                    # Optional: serialize dict as JSON
                    pdb_string = json.dumps(pdb_string, indent=2)
                output_path = Path(f"{label}.pdb")
                output_path.write_text(pdb_string)
                print(f"Saved: {output_path}")

        else:
            print("Server Error:")
            print(response.text)

    except Exception as e:
        print(f"Request failed: {e}")

if __name__ == "__main__":
    main()
