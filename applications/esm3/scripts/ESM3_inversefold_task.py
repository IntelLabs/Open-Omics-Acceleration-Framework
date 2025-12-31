# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import os
import torch
import argparse
import time
from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein, ESMProteinError
from esm.sdk.api import GenerationConfig
from esm.utils.structure.protein_complex import ProteinComplex

def write_fasta(filename, sequences):
    """Writes generated sequences to a FASTA file."""
    with open(filename, "w") as f:
        for i, seq in enumerate(sequences, 1):
            f.write(f">sampled_seq_{i}\n{seq}\n")

def inversefold_protein(client: ESM3InferenceClient, protein ,pdb_file: str, args,output_dir: str = None):
    """Runs inverse folding and saves the output as a FASTA file."""
    print(f"Processing {pdb_file}...")

    with torch.amp.autocast(device_type="cpu", enabled=args.bf16):
        protein.sequence = None
        protein.sasa = None
        protein.function_annotations = None

        if args.batch_run:
            num_sequences = args.num_sequences
            generations = client.batch_generate(
                inputs=[protein] * num_sequences,
                configs=[GenerationConfig(
                    track="sequence",
                    schedule=args.schedule,
                    strategy=args.strategy,
                    num_steps=args.num_steps,
                    temperature=args.temperature,
                    temperature_annealing=args.temperature_annealing,
                    top_p=args.top_p,
                    condition_on_coordinates_only=args.condition_on_coordinates_only
                )] * num_sequences,
            )

            sequences = []
            errors = []
            for i, gen_protein in enumerate(generations):
                if isinstance(gen_protein, ESMProteinError):
                    errors.append((i, gen_protein))
                else:
                    sequences.append(gen_protein.sequence)

            if errors:
                raise RuntimeError(f"Batch generation errors: {errors}")
            if output_dir:
                output_fasta_path = os.path.join(output_dir, f"{os.path.basename(pdb_file).replace('.pdb', '')}_batch.fasta")
                write_fasta(output_fasta_path, sequences)
            return sequences

        else:
            inv_folded_protein = client.generate(
                protein,
                GenerationConfig(
                    track="sequence",
                    schedule=args.schedule,
                    strategy=args.strategy,
                    num_steps=args.num_steps,
                    temperature=args.temperature,
                    temperature_annealing=args.temperature_annealing,
                    top_p=args.top_p,
                    condition_on_coordinates_only=args.condition_on_coordinates_only
                ),
            )

            if isinstance(inv_folded_protein, ESMProteinError):
                raise RuntimeError(f"Error: {str(inv_folded_protein)}")
            if output_dir:
                output_fasta_path = os.path.join(output_dir, os.path.basename(pdb_file).replace(".pdb", ".fasta"))
                write_fasta(output_fasta_path, [inv_folded_protein.sequence])
            return  [inv_folded_protein.sequence]

def processing_pdb(client: ESM3InferenceClient, pdb_file: str, args,output_dir: str = None):
    """Runs protein folding and saves the output as a PDB file in the specified directory."""
    print(f"Processing {pdb_file}...")

    if args.protein_complex:
        protein = ProteinComplex.from_pdb(pdb_file)
        protein = ESMProtein.from_protein_complex(protein)
        result = inversefold_protein(client,protein,pdb_file,args,output_dir)
    else:
        protein = ESMProtein.from_pdb(pdb_file)
        result= inversefold_protein(client,protein,pdb_file,args,output_dir)
    return result

def main(args):
    """Processes a single PDB file."""
    if not os.path.exists(args.pdb_file) or not args.pdb_file.endswith(".pdb"):
        print(f"Invalid PDB file: {args.pdb_file}")
        return

    os.makedirs(args.output_dir, exist_ok=True)

    client = ESM3InferenceClient() if os.getenv("ESM_API_KEY") else ESM3.from_pretrained(
        "esm3_sm_open_v1", bf16=args.bf16
    )

    inter_time = time.time()
    processing_pdb(client, args.pdb_file, args,args.output_dir)
    if args.timing:
        print(f"Inference time = {time.time() - inter_time} seconds")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ESM3 protein inverse folding on a PDB file.")
    parser.add_argument("pdb_file", type=str, help="Path to the input PDB file.")
    parser.add_argument("output_dir", type=str, help="Directory to save the output FASTA file.")
    parser.add_argument("--bf16", action="store_true", help="Enable bf16 inference.")
    parser.add_argument("--timing", action="store_true", help="Enable timing for inference.")
    parser.add_argument("--batch_run", action="store_true", help="Run in batch mode.")
    parser.add_argument("--schedule", type=str, choices=["cosine", "linear"], default="cosine", help="Schedule type (cosine or linear).")
    parser.add_argument("--strategy", type=str, choices=["random", "entropy"], default="entropy", help="Unmasking strategy (random or entropy).")
    parser.add_argument("--num_steps", type=int, default=1, help="Number of steps for generation.")
    parser.add_argument("--temperature", type=float, default=1.0, help="Sampling temperature.")
    parser.add_argument("--temperature_annealing", action="store_true", help="Enable temperature annealing.")
    parser.add_argument("--top_p", type=float, default=1.0, help="Top-p sampling value.")
    parser.add_argument("--condition_on_coordinates_only", action="store_true", help="Condition only on coordinates.")
    parser.add_argument("--num_sequences", type=int, default=8, help="Number of sequences to generate in batch mode.")
    parser.add_argument("--protein_complex", action="store_true", help="Enable prediction for protein complexes (multi-chain structure) using a multi-chain FASTA file input.")

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    if args.timing:
        start_time = time.time()
    main(args)
    if args.timing:
        print(f"Complete run time = {time.time() - start_time} seconds")
