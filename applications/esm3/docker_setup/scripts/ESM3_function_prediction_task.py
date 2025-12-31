# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import os
import torch
import argparse
from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein
from esm.sdk.api import GenerationConfig
import time
import csv
from esm.utils.structure.protein_complex import ProteinComplex

def function_protein(client: ESM3InferenceClient, protein,pdb_file: str,args,output_dir: str = None):
    """Runs protein function prediction and saves the output as a CSV file."""
    print(f"Processing {pdb_file}...")


    with torch.amp.autocast(device_type="cpu", enabled=args.bf16):
        protein.function_annotations = None

        # Inline GenerationConfig inside client.generate()
        protein_with_function = client.generate(
            protein,
            GenerationConfig(
                track="function",
                schedule=args.schedule,
                strategy=args.strategy,
                num_steps=args.num_steps,
                temperature=args.temperature,
                temperature_annealing=args.temperature_annealing,
                top_p=args.top_p,
                condition_on_coordinates_only=args.condition_on_coordinates_only,
            ),
        )

        assert isinstance(protein_with_function, ESMProtein), f"Unexpected output: {protein_with_function}"
        if output_dir:
            output_csv = os.path.join(output_dir, f"{os.path.basename(pdb_file).replace('.pdb', '')}.csv")

            with open(output_csv, "w", newline="") as f:
                writer = csv.writer(f, delimiter="\t")
                writer.writerow(["Label", "Start", "End"])
                if protein_with_function.function_annotations:
                    for annotation in protein_with_function.function_annotations:
                        writer.writerow([annotation.label, annotation.start, annotation.end])

            print(f"Function annotations saved as {output_csv}")
    return protein_with_function.function_annotations

def processing_pdb(client: ESM3InferenceClient, pdb_file: str, args,output_dir: str = None ):
    """Runs protein folding and saves the output as a PDB file in the specified directory."""
    print(f"Processing {pdb_file}...")

    if args.protein_complex:
        protein = ProteinComplex.from_pdb(pdb_file)
        protein = ESMProtein.from_protein_complex(protein)
        result=function_protein(client,protein,pdb_file,args,output_dir)
    else:
        protein = ESMProtein.from_pdb(pdb_file)
        result=function_protein(client,protein,pdb_file,args,output_dir)
    return result

def main(args):
    if not os.path.exists(args.pdb_file) or not args.pdb_file.endswith(".pdb"):
        print(f"Invalid PDB file: {args.pdb_file}")
        return
    os.makedirs(args.output_dir, exist_ok=True)

    client = ESM3InferenceClient() if os.getenv("ESM_API_KEY") else ESM3.from_pretrained("esm3_sm_open_v1", bf16=args.bf16)

    infer_time = time.time()
    processing_pdb(client, args.pdb_file, args,args.output_dir)
    if args.timing:
        print(f"Inference time = {time.time() - infer_time} seconds")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ESM3 protein function annotation on a single PDB file.")
    parser.add_argument("pdb_file", type=str, help="Path to the input PDB file.")
    parser.add_argument("output_dir", type=str, help="Directory to save the output CSV file.")
    parser.add_argument("--bf16", action="store_true", help="Enable bf16 inference.")
    parser.add_argument("--timing", action="store_true", help="Enable timing for inference.")
    parser.add_argument("--schedule", type=str, choices=["cosine", "linear"], default="cosine", help="Schedule type (cosine or linear).")
    parser.add_argument("--strategy", type=str, choices=["random", "entropy"], default="entropy", help="Unmasking strategy (random or entropy).")
    parser.add_argument("--num_steps", type=int, default=1, help="Number of steps for generation.")
    parser.add_argument("--temperature", type=float, default=1.0, help="Sampling temperature.")
    parser.add_argument("--temperature_annealing", action="store_true", help="Enable temperature annealing.")
    parser.add_argument("--top_p", type=float, default=1.0, help="Top-p sampling value.")
    parser.add_argument("--condition_on_coordinates_only", action="store_true", help="Condition only on coordinates.")
    parser.add_argument("--protein_complex", action="store_true", help="Enable prediction for protein complexes (multi-chain structure) using a multi-chain FASTA file input.")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    if args.timing:
        start_time = time.time()
    main(args)
    if args.timing:
        print(f"Complete run time = {time.time() - start_time:.2f} seconds")
