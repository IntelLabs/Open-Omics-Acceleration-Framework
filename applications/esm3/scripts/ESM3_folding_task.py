# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import os
import re
import torch
import argparse
from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein, GenerationConfig
from esm.utils.structure.protein_chain import ProteinChain
from esm.utils.structure.protein_complex import ProteinComplex
from Bio import SeqIO
import numpy as np
import time

def make_protein_chain(record):
    """Create a ProteinChain object from a FASTA record."""
    sequence = str(record.seq)
    L = len(sequence)
    try:
        if "|" in record.description:
            parts = record.description.split('|')
            if len(parts) > 1:
                chain_id = parts[1].strip().split()[-1]
            else:
                chain_id = 'A'
        else:
            chain_id = 'A'
    except Exception as e:
        print(f"Chain ID extraction failed: {e}")
        chain_id = 'A'
    # Creating dummy values for atom37 and other fields
    protein_chain = ProteinChain(
        id=record.id,
        sequence=sequence,
        chain_id=chain_id,
        entity_id=None,  # Can be assigned if needed
        residue_index=np.arange(L),
        insertion_code=np.array([""] * L),  # Empty insertion codes
        atom37_positions=np.zeros((L, 37, 3), dtype=np.float32),  # Dummy atom37 positions
        atom37_mask=np.zeros((L, 37), dtype=np.float32),  # Dummy atom37 mask
        confidence=np.zeros((L,), dtype=np.float32)  # Dummy confidence values
    )
    return protein_chain
# The following function is adapted from:
# https://github.com/facebookresearch/esm (License: MIT)
def read_fasta(path, keep_gaps=True, keep_insertions=True, to_upper=False):
    """Reads a FASTA file and returns a generator of (header, sequence) tuples."""
    with open(path, "r") as f:
        for result in read_alignment_lines(
            f, keep_gaps=keep_gaps, keep_insertions=keep_insertions, to_upper=to_upper
        ):
            yield result

# The following function is adapted from:
# https://github.com/facebookresearch/esm (License: MIT)
def read_alignment_lines(lines, keep_gaps=True, keep_insertions=True, to_upper=False):
    """Parses the FASTA file lines and processes sequence formatting."""
    seq = desc = None

    def parse(s):
        if not keep_gaps:
            s = re.sub("-", "", s)
        if not keep_insertions:
            s = re.sub("[a-z]", "", s)
        return s.upper() if to_upper else s

    for line in lines:
        if len(line) > 0 and line[0] == ">":
            if seq is not None:
                yield desc, parse(seq)
            desc = line.strip().lstrip(">")
            seq = ""
        else:
            assert isinstance(seq, str)
            seq += line.strip()
    assert isinstance(seq, str) and isinstance(desc, str)
    yield desc, parse(seq)

def fold_protein(client: ESM3InferenceClient, protein, args, output_dir, fasta_file):
    """Performs protein folding using the ESM3 model and saves the result as a PDB file."""

    with torch.amp.autocast(device_type="cpu", enabled=args.bf16):
        protein.coordinates = None
        protein.function_annotations = None
        protein.sasa = None
        sequence_length = len(protein.sequence)
        num_steps = args.num_steps if args.num_steps else max(1, sequence_length // 16)

        folded_protein = client.generate(
                protein,
                GenerationConfig(
                    track="structure",
                    schedule=args.schedule,
                    strategy=args.strategy,
                    num_steps=num_steps,
                    temperature=args.temperature,
                    temperature_annealing=args.temperature_annealing,
                    top_p=args.top_p,
                    condition_on_coordinates_only=args.condition_on_coordinates_only,
                ),
            )

        assert isinstance(folded_protein, ESMProtein), f"ESMProtein was expected but got {protein}"
        
        return folded_protein
    
def save_pdb(folding_outuput,fasta_name,output_dir):
    
    output_pdb_path = os.path.join(output_dir, f"{fasta_name}.pdb")
    folding_outuput.to_pdb(output_pdb_path)
    print(f"Saved folded protein to {output_pdb_path}")
    mean_plddt = folding_outuput.plddt.mean().item()
    print(f"Mean pLDDT Score: {mean_plddt:.4f}")

def processing_fasta(client: ESM3InferenceClient, fasta_file: str, output_dir: str, args):
    """Runs protein folding and saves the output as a PDB file in the specified directory."""
    print(f"Processing {fasta_file}...")

    if args.protein_complex:
        protein_chains = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            protein_chain = make_protein_chain(record)
            protein_chains.append(protein_chain)

        protein = ProteinComplex.from_chains(protein_chains)
        protein = ESMProtein.from_protein_complex(protein)
        folding_outuput = fold_protein(client, protein, args, output_dir, fasta_file)
        fasta_name = os.path.splitext(os.path.basename(fasta_file))[0]
        save_pdb(folding_outuput,fasta_name,output_dir)
    else:
        fasta_entries = sorted(read_fasta(fasta_file), key=lambda header_seq: len(header_seq[1]))
        for entry in fasta_entries:
            header, sequence = entry
            protein = ESMProtein(sequence=sequence)
            folding_outuput = fold_protein(client, protein, args, output_dir, fasta_file)
            fasta_name = os.path.splitext(os.path.basename(header))[0]
            save_pdb(folding_outuput,fasta_name,output_dir)


def main(args):
    """Processes a single FASTA file."""
    if not os.path.exists(args.fasta_file) or not args.fasta_file.endswith(".fasta"):
        print(f"Invalid Fasta file: {args.fasta_file}")
        return

    if os.environ.get("ESM_API_KEY", ""):
        print("ESM_API_KEY found. Using model from Forge API...")
        client = ESM3InferenceClient()
    else:
        print("No ESM_API_KEY found. Loading model locally...")
        client = ESM3.from_pretrained("esm3_sm_open_v1", bf16=args.bf16)

    if args.timing:
        inter_time = time.time()
    processing_fasta(client, args.fasta_file, args.output_dir, args)
    if args.timing:
        print(f"inference time = {time.time() - inter_time} seconds")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ESM3 protein folding for a single FASTA file.")
    parser.add_argument("fasta_file", type=str, help="Path to the input FASTA file.")
    parser.add_argument("output_dir", type=str, help="Directory to save the output PDB file.")
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
        print(f"complete run time = {time.time() - start_time} seconds")
