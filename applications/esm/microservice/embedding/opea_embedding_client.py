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
    parser.add_argument("--port", default="9001", help="Microservice port")
    parser.add_argument("--fasta_file", required=True, help="Path to input FASTA file")
    parser.add_argument("--toks_per_batch", type=int, default=4096, help="Tokens per batch")
    parser.add_argument("--repr_layers", nargs="+", type=int, default=[-1], help="Representation layers")
    parser.add_argument("--include", nargs="+", choices=["mean", "per_tok", "bos", "contacts"], default=["mean"], help="Types of output to include")
    args = parser.parse_args()

    # Construct URL
    url = f"http://{args.host}:{args.port}/v1/esmembedding"

    # Build payload with base64-encoded FASTA
    with open(args.fasta_file, "rb") as f:
        fasta_b64 = base64.b64encode(f.read()).decode("utf-8")

    payload = {
        "fasta_base64": fasta_b64,
        "toks_per_batch": args.toks_per_batch,
        "repr_layers": args.repr_layers,
        "include": args.include,
    }
    print("Sending POST to:", url)

    try:
        response = requests.post(url, json=payload)
        print(f"Status Code: {response.status_code}")

        if response.status_code == 200:
            data = response.json()
            for item in data.get("result", []):
                label = item["label"]
                result = {"label": label}
                if "mean_representations" in item:
                    result["mean_representations"] = {
                        int(layer): torch.tensor(tensor_list)
                        for layer, tensor_list in item["mean_representations"].items()
                    }
                if "bos_representations" in item:
                    result["bos_representations"] = {
                        int(layer): torch.tensor(tensor_list)
                        for layer, tensor_list in item["bos_representations"].items()
                    }
                if "representations" in item:
                    result["representations"] = {
                        int(layer): torch.tensor(tensor_list)
                        for layer, tensor_list in item["representations"].items()
                    }
                if "contacts" in item:
                    result["contacts"] = torch.tensor(item["contacts"])

                # Save result per item
                save_path = Path(f"{label}.pt")
                torch.save(result, save_path)
                print(f"Saved: {save_path}")
        else:
            print("Server Error:")
            print(response.text)

    except Exception as e:
        print(f"Request failed: {e}")

if __name__ == "__main__":
    main()
