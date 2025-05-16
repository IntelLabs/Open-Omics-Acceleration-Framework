import os
import re
import torch
import argparse
from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein
from esm.utils.structure.protein_chain import ProteinChain
from esm.utils.structure.protein_complex import ProteinComplex
from Bio import SeqIO
import numpy as np
from esm.sdk.api import LogitsConfig, LogitsOutput
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

def logits_embedding(client: ESM3InferenceClient, protein, args):
    """Performs protein folding using the ESM3 model and saves the result as a PDB file."""

    with torch.amp.autocast(device_type="cpu", enabled=args.bf16):
        
        protein_tensor = client.encode(protein)
        logits_output = client.logits(protein_tensor, 
                                        LogitsConfig(
                                                sequence=args.sequence,
                                                structure=args.structure,
                                                secondary_structure=args.secondary_structure,
                                                sasa=args.sasa,
                                                function=args.function,
                                                residue_annotations=args.residue_annotations,
                                                return_embeddings=args.return_embeddings,
                                                return_hidden_states=args.return_hidden_states,
                                                ith_hidden_layer=args.ith_hidden_layer))
        
        assert isinstance(logits_output, LogitsOutput), f"Expected LogitsOutput but got {logits_output}"
        
          
        return logits_output

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
        logits_output = logits_embedding(client, protein, args)
        fasta_name = os.path.splitext(os.path.basename(fasta_file))[0]
        torch.save(logits_output, os.path.join(output_dir, f"{fasta_name}.pt"))
    else:
        fasta_entries = sorted(read_fasta(fasta_file), key=lambda header_seq: len(header_seq[1]))
        for entry in fasta_entries:
            header, sequence = entry
            protein = ESMProtein(sequence=sequence)
            logits_output = logits_embedding(client, protein, args)
            fasta_name = os.path.splitext(os.path.basename(header))[0]
            torch.save(logits_output, os.path.join(output_dir, f"{fasta_name}.pt"))

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
        infer_time = time.time()
    processing_fasta(client, args.fasta_file, args.output_dir, args)
    if args.timing:
        print(f"inference time = {time.time() - infer_time} seconds")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ESM3 protein inference for a single FASTA file.")
    parser.add_argument("fasta_file", type=str, help="Path to the input FASTA file.")
    parser.add_argument("output_dir", type=str, help="Directory to save the output.")
    parser.add_argument("--bf16", action="store_true", help="Enable bf16 inference.")
    parser.add_argument("--timing", action="store_true", help="Enable timing for inference.")
    parser.add_argument("--sequence", action="store_true", help="Enable sequence logits.")
    parser.add_argument("--structure", action="store_true", help="Enable structure logits.")
    parser.add_argument("--secondary_structure", action="store_true", help="Enable secondary structure logits.")
    parser.add_argument("--sasa", action="store_true", help="Enable SASA logits.")
    parser.add_argument("--function", action="store_true", help="Enable function logits.")
    parser.add_argument("--residue_annotations", action="store_true", help="Enable residue annotation logits.")
    parser.add_argument("--return_embeddings", action="store_true", help="Enable embeddings output.")
    parser.add_argument("--return_hidden_states", action="store_true", help="Enable hidden state output.")
    parser.add_argument("--ith_hidden_layer", type=int, default=-1, help="Specify which hidden layer to return.")
    parser.add_argument("--protein_complex", action="store_true", help="Enable prediction for protein complexes (multi-chain structure) using a multi-chain FASTA file input.")   


    args = parser.parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)
    if args.timing:
        start_time = time.time()
    main(args)
    if args.timing:
        print(f"complete run time = {time.time() - start_time} seconds")
