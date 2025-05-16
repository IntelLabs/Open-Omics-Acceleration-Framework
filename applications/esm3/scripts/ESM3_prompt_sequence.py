import os
import re
import torch
import argparse
import time
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

def read_fasta(path, keep_gaps=True, keep_insertions=True, to_upper=False):
    """Reads a FASTA file and returns a generator of (header, sequence) tuples."""
    with open(path, "r") as f:
        for result in read_alignment_lines(
            f, keep_gaps=keep_gaps, keep_insertions=keep_insertions, to_upper=to_upper
        ):
            yield result

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

def write_fasta(filename, sequences):
    """Writes sequences to a FASTA file."""
    print(f"Writing output to: {filename}")
    with open(filename, "w") as f:
        for i, seq in enumerate(sequences, 1):
            f.write(f">generated_seq_{i}\n{seq}\n")


def prompt_protein(client: ESM3InferenceClient, protein, args, output_dir, fasta_file):
    """Performs protein folding using the ESM3 model and saves the result as a PDB file."""

    with torch.amp.autocast(device_type="cpu", enabled=args.bf16):
        prompt_protein = client.generate(protein, GenerationConfig(
            track="sequence",
            schedule=args.schedule,
            num_steps=args.num_steps,
            strategy=args.strategy,
            temperature=args.temperature,
            top_p=args.top_p,
            condition_on_coordinates_only=args.condition_on_coordinates_only))
        assert isinstance(protein, ESMProtein), "Expected ESMProtein output"
    return prompt_protein


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
        prompt_outuput = prompt_protein(client, protein, args, output_dir, fasta_file)
        fasta_name = os.path.splitext(os.path.basename(fasta_file))[0]
        output_fasta_path = os.path.join(output_dir, f"{fasta_name}.fasta")
        write_fasta(output_fasta_path, [prompt_outuput.sequence])
    else:
        fasta_entries = sorted(read_fasta(fasta_file), key=lambda header_seq: len(header_seq[1]))
        
        for entry in fasta_entries:
            header, sequence = entry
            protein = ESMProtein(sequence=sequence)
            prompt_outuput = prompt_protein(client, protein, args, output_dir, fasta_file)
            output_fasta_path = os.path.join(output_dir, f"{header}.fasta")
            write_fasta(output_fasta_path, [prompt_outuput.sequence])

def main(args):
    """Main execution function."""
    if not os.path.exists(args.fasta_file) or not args.fasta_file.endswith(".fasta"):
        print(f"Invalid Fasta file: {args.fasta_file}")
        return
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Load ESM model
    if os.environ.get("ESM_API_KEY"):
        print("Using ESM API Client...")
        client = ESM3InferenceClient()
    else:
        print("Loading ESM3 model locally...")
        client = ESM3.from_pretrained("esm3_sm_open_v1", bf16=args.bf16)

    # Run inference
    if args.timing:
        infer_time = time.time()
    processing_fasta(client, args.fasta_file, args.output_dir, args)
    if args.timing:
        print(f"Inference time: {time.time() - infer_time:.2f} seconds")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ESM3 protein folding on a single FASTA file.")
    parser.add_argument("fasta_file", type=str, help="Path to the input FASTA file.")
    parser.add_argument("output_dir", type=str, help="Directory to save the output FASTA file.")
    parser.add_argument("--bf16", action="store_true", help="Enable bf16 inference.")
    parser.add_argument("--timing", action="store_true", help="Measure inference time.")
    parser.add_argument("--schedule", type=str, choices=["cosine", "linear"], default="cosine", help="Schedule type.")
    parser.add_argument("--strategy", type=str, choices=["random", "entropy"], default="entropy", help="Unmasking strategy.")
    parser.add_argument("--num_steps", type=int, default=1, help="Number of inference steps.")
    parser.add_argument("--temperature", type=float, default=1.0, help="Sampling temperature.")
    parser.add_argument("--top_p", type=float, default=1.0, help="Top-p sampling value.")
    parser.add_argument("--condition_on_coordinates_only", action="store_true", help="Condition only on coordinates.")
    parser.add_argument("--protein_complex", action="store_true", help="Enable prediction for protein complexes (multi-chain structure) using a multi-chain FASTA file input.")   

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    if args.timing:
        start_time = time.time()
    main(args)
    if args.timing:
        print(f"Complete run time = {time.time() - start_time} seconds")