#!/usr/bin/env python
# encoding: utf-8

# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import base64
import argparse
import requests
import time
import os
import json
from pathlib import Path


def read_pdb_file(pdb_path: Path) -> str:
    if not pdb_path.exists():
        raise FileNotFoundError(f"File not found: {pdb_path}")
    return "\n".join(
        line for line in pdb_path.read_text().splitlines()
        #if line.startswith(("ATOM", "HETATM"))
    )


def encode_file_to_base64(file_path):
    """Encode file content to base64 if file exists"""
    if file_path and os.path.exists(file_path):
        with open(file_path, 'rb') as f:
            return base64.b64encode(f.read()).decode('utf-8')
    return ""


def save_decoded_file(content, save_path):
    """Decode base64 content and save to file."""
    try:
        with open(save_path, "wb") as f:
            f.write(base64.b64decode(content))
        print(f"💾 Saved {save_path}")
    except Exception as e:
        print(f"❌ Failed to save {save_path}: {e}")


def get_base_name(args):
    """Determine base name from input file."""
    path = args.pdb_path or args.jsonl_path or args.path_to_fasta
    return os.path.splitext(os.path.basename(path))[0] if path else "output"


def normalize_to_list(item):
    """Ensure value is always a list."""
    if item is None:
        return []
    return item if isinstance(item, list) else [item]

def get_output_filename(file_type, base_name, index, total):
    """Generate correct filename depending on file_type."""
    print(f"DEBUG: file_type='{file_type}', base_name='{base_name}', index={index}")  # ADD THIS
    # special rule for score_only
    if file_type == "score_only":
        return f"{base_name}_pdb.npz" if index == 0 else f"{base_name}_fasta_{index}.npz"

    # Generate filename based on file_type
    if file_type == "pdb":
        return f"{base_name}_sample_{index+1}.pdb"
    #elif file_type == "fa":
    #    return f"{base_name}_sample_{index+1}.fas"

    elif file_type == "fa" or file_type == "seqs":  # Combined condition
        return f"{base_name}_sample_{index+1}.fa"

    elif file_type == "npz":
        return f"{base_name}_sample_{index+1}.npz"

    elif file_type in ["probs", "scores"]:
        return f"{base_name}_{file_type}_{index+1}.json"

    elif file_type == "conditional_probs_only":
        return f"{base_name}_conditional_probs_only_{index+1}.npz"

    elif file_type == "unconditional_probs_only":
        return f"{base_name}_unconditional_probs_only_{index+1}.npz"
    else:
        print(f"DEBUG: Falling back to .txt for file_type='{file_type}'")
        return f"{base_name}_{file_type}_{index+1}.txt"



def handle_results(result, args):
    """Process all results returned by the API."""
    if "results" not in result:
        print(f"❌ Error: No 'results' key in API response")
        print(f"Full response: {result}")
        return

    base_name = get_base_name(args)
    results = result["results"]

    print(f"Base name: {base_name}")
    print(f"Received file types: {list(results.keys())}")

    # Create output folder if it doesn't exist
    os.makedirs(args.out_folder, exist_ok=True)

    for file_type, encoded_contents in results.items():
        if not encoded_contents:
            print(f"⚠️ Warning: No content for file type '{file_type}'")
            continue

        contents = normalize_to_list(encoded_contents)

        for i, content in enumerate(contents):
            if not content:
                print(f"⚠️ Warning: Empty content at index {i} for '{file_type}'")
                continue

            filename = get_output_filename(file_type, base_name, i, len(contents))
            save_path = os.path.join(args.out_folder, filename)
            save_decoded_file(content, save_path)


def call_api(url, payload, args):
    """Perform API call and process results."""
    try:
        response = requests.post(url, json=payload, timeout=3600)
        print(f"API Response status: {response.status_code}")

        if response.status_code != 200:
            print(f"❌ Inference failed: {response.status_code}")
            print(f"Error: {response.text}")
            return

        result = response.json()
        print("✅ Inference completed successfully")

        handle_results(result, args)

        # Save entire JSON response
        response_path = os.path.join(args.out_folder, "response.json")
        with open(response_path, "w") as f:
            json.dump(result, f, indent=2)
        print(f"💾 Saved full response to {response_path}")

    except requests.exceptions.RequestException as e:
        print(f"❌ Request failed: {e}")
    except Exception as e:
        print(f"❌ Unexpected error: {e}")


def main():
    argparser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    # Required arguments
    argparser.add_argument("--host", type=str, required=True, help="Server host")
    argparser.add_argument("--port", type=int, required=True, help="Server port")
    argparser.add_argument("--out_folder", type=str, required=True, help="Output folder")
    # ProteinMPNN parameters
    argparser.add_argument("--suppress_print", type=int, default=0, help="0 for False, 1 for True")
    argparser.add_argument("--seed", type=int, default=0, help="Random seed")
    argparser.add_argument("--save_score", type=int, default=0, help="Save score to file")
    argparser.add_argument("--save_probs", type=int, default=0, help="Save probabilities to file")
    argparser.add_argument("--score_only", type=int, default=0, help="Only score the structure")
    argparser.add_argument("--path_to_fasta", type=str, default="", help="Path to fasta file")
    argparser.add_argument("--conditional_probs_only", type=int, default=0, help="Calculate conditional probabilities")
    argparser.add_argument("--conditional_probs_only_backbone", type=int, default=0, help="Conditional backbone probabilities")
    argparser.add_argument("--unconditional_probs_only", type=int, default=0, help="Calculate unconditional probabilities")
    argparser.add_argument("--backbone_noise", type=float, default=0.00, help="Noise to add to backbone atoms")
    argparser.add_argument("--num_seq_per_target", type=int, default=1, help="Number of sequences to generate")
    argparser.add_argument("--batch_size", type=int, default=1, help="Batch size")
    argparser.add_argument("--max_length", type=int, default=200000, help="Maximum sequence length")
    argparser.add_argument("--sampling_temp", type=str, default="0.1", help="Sampling temperature")
    argparser.add_argument("--pdb_path", type=str, default='', help="Path to single PDB")
    argparser.add_argument("--pdb_path_chains", type=str, default='', help="Chains to design")
    argparser.add_argument("--jsonl_path", type=str, default='', help="Path to jsonl file")
    argparser.add_argument("--chain_id_jsonl", type=str, default='', help="Path to chain id jsonl")
    argparser.add_argument("--fixed_positions_jsonl", type=str, default='', help="Path to fixed positions jsonl")
    argparser.add_argument("--omit_AAs", type=str, default="X", help="Amino acids to omit")
    argparser.add_argument("--bias_AA_jsonl", type=str, default='', help="Path to bias AA jsonl")
    argparser.add_argument("--bias_by_res_jsonl", type=str, default='', help="Path to bias by residue jsonl")
    argparser.add_argument("--omit_AA_jsonl", type=str, default='', help="Path to omit AA jsonl")
    argparser.add_argument("--pssm_jsonl", type=str, default='', help="Path to PSSM jsonl")
    argparser.add_argument("--pssm_multi", type=float, default=0.0, help="PSMM weight")
    argparser.add_argument("--pssm_threshold", type=float, default=0.0, help="PSMM threshold")
    argparser.add_argument("--pssm_log_odds_flag", type=int, default=0, help="Use log odds flag")
    argparser.add_argument("--pssm_bias_flag", type=int, default=0, help="Use bias flag")
    argparser.add_argument("--tied_positions_jsonl", type=str, default='', help="Path to tied positions jsonl")
    argparser.add_argument("--chains_to_design", type=str, default="A B", help="Chains to design (space-separated)")

    args = argparser.parse_args()

    start_time = time.time()

    # Create output folder
    os.makedirs(args.out_folder, exist_ok=True)

    url = f"http://{args.host}:{args.port}/v1/protein_mpnn"

    # Validate input files
    input_provided = False
    jsonl_base64 = None
    pdb_base64 = None
    fasta_b64 = None

    if args.jsonl_path:
        if not os.path.exists(args.jsonl_path):
            raise FileNotFoundError(f"JSONL file not found: {args.jsonl_path}")
        with open(args.jsonl_path, 'rb') as f:
            content = f.read()
        jsonl_base64 = base64.b64encode(content).decode('utf-8')
        input_provided = True

    if args.pdb_path:
        if not os.path.exists(args.pdb_path):
            raise FileNotFoundError(f"PDB file not found: {args.pdb_path}")
        pdb_str = read_pdb_file(Path(args.pdb_path))
        pdb_base64 = base64.b64encode(pdb_str.encode('utf-8')).decode('utf-8')
        input_provided = True

    if args.path_to_fasta:
        if not os.path.exists(args.path_to_fasta):
            raise FileNotFoundError(f"FASTA file not found: {args.path_to_fasta}")
        with open(args.path_to_fasta, "rb") as f:
            fasta_content = f.read()
        fasta_b64 = base64.b64encode(fasta_content).decode('utf-8')
        input_provided = True

    if not input_provided:
        raise ValueError("No input provided. Use at least one of: --jsonl_path, --pdb_path, or --path_to_fasta")

    # Encode optional files
    chain_id_jsonl_b64 = encode_file_to_base64(args.chain_id_jsonl)
    fixed_positions_jsonl_b64 = encode_file_to_base64(args.fixed_positions_jsonl)
    pssm_jsonl_b64 = encode_file_to_base64(args.pssm_jsonl)
    omit_AA_jsonl_b64 = encode_file_to_base64(args.omit_AA_jsonl)
    bias_AA_jsonl_b64 = encode_file_to_base64(args.bias_AA_jsonl)
    tied_positions_jsonl_b64 = encode_file_to_base64(args.tied_positions_jsonl)
    bias_by_res_jsonl_b64 = encode_file_to_base64(args.bias_by_res_jsonl)

    # Prepare payload
    payload = {
        "suppress_print": args.suppress_print,
        "seed": args.seed,
        "save_score": args.save_score,
        "save_probs": args.save_probs,
        "score_only": args.score_only,
        "path_to_fasta": fasta_b64,
        "conditional_probs_only": args.conditional_probs_only,
        "conditional_probs_only_backbone": args.conditional_probs_only_backbone,
        "unconditional_probs_only": args.unconditional_probs_only,
        "backbone_noise": args.backbone_noise,
        "num_seq_per_target": args.num_seq_per_target,
        "batch_size": args.batch_size,
        "max_length": args.max_length,
        "sampling_temp": args.sampling_temp,
        "pdb_path": pdb_base64,
        "pdb_path_chains": args.pdb_path_chains,
        "jsonl_path": jsonl_base64,
        "chain_id_jsonl": chain_id_jsonl_b64,
        "fixed_positions_jsonl": fixed_positions_jsonl_b64,
        "omit_AAs": args.omit_AAs,
        "bias_AA_jsonl": bias_AA_jsonl_b64,
        "bias_by_res_jsonl": bias_by_res_jsonl_b64,
        "omit_AA_jsonl": omit_AA_jsonl_b64,
        "pssm_jsonl": pssm_jsonl_b64,
        "pssm_multi": args.pssm_multi,
        "pssm_threshold": args.pssm_threshold,
        "pssm_log_odds_flag": args.pssm_log_odds_flag,
        "pssm_bias_flag": args.pssm_bias_flag,
        "tied_positions_jsonl": tied_positions_jsonl_b64,
        "chains_to_design": args.chains_to_design,
        "checkpoint": {}
    }

    # Call API
    call_api(url, payload, args)

    end_time = time.time()
    print(f"⏱️  Total execution time: {end_time - start_time:.2f} seconds")


if __name__ == "__main__":
    main()
