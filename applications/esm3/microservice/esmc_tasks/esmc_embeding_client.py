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
    parser.add_argument("--port", default="9009", help="Microservice port")
    parser.add_argument("--fasta_file", type=str, help="Path to the FASTA file.")
    parser.add_argument("--sequence", action="store_true", help="Enable sequence logits")
    parser.add_argument("--structure", action="store_true", help="Enable structure logits")
    parser.add_argument("--secondary_structure", action="store_true", help="Enable secondary structure logits")
    parser.add_argument("--sasa", action="store_true", help="Enable SASA logits")
    parser.add_argument("--function", action="store_true", help="Enable function logits")
    parser.add_argument("--residue_annotations", action="store_true", help="Enable residue annotations")
    parser.add_argument("--return_embeddings", action="store_true", help="Return embeddings")
    parser.add_argument("--return_hidden_states", action="store_true", help="Return hidden states")
    parser.add_argument("--ith_hidden_layer", type=int, default=-1, help="Specify which hidden layer to return")
    parser.add_argument("--protein_complex", action="store_true", help="Enable prediction for protein complexes (multi-chain structure) using a multi-chain FASTA file input.")   

    args = parser.parse_args()

    # Construct URL
    url = f"http://{args.host}:{args.port}/v1/esmcembedding"

    # Build payload with base64-encoded FASTA
    with open(args.fasta_file, "rb") as f:
        fasta_b64 = base64.b64encode(f.read()).decode("utf-8")

    payload = {
        "fasta_base64": fasta_b64,
        "sequence": args.sequence,
        "structure": args.structure,
        "secondary_structure": args.secondary_structure,
        "sasa": args.sasa,
        "function": args.function,
        "residue_annotations": args.residue_annotations,
        "return_embeddings": args.return_embeddings,
        "return_hidden_states": args.return_hidden_states,
        "ith_hidden_layer": args.ith_hidden_layer,
        "protein_complex": args.protein_complex,
        "port":args.port

    }
    print("Sending POST to:", url)

    try:
        response = requests.post(url, json=payload)
        print(f"Status Code: {response.status_code}")

        if response.status_code == 200:
            data = response.json()
            encoded_result = data.get("result") # Base64 string
            decoded_json = base64.b64decode(encoded_result).decode("utf-8")
            decoded_data = json.loads(decoded_json)
            
            for item in decoded_data:
                for name, content in item.items():
                    filename = f"{name}.pt"
                    torch.save(content, filename)
                    print(f" Saved {filename}")

        else:
            print("Server Error:")
            print(response.text)

    except Exception as e:
        print(f"Request failed: {e}")

if __name__ == "__main__":
    main()