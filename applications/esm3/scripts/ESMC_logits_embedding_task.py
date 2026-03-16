# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import os
import re
import time
import argparse
import torch
from esm.models.esmc import ESMC
from esm.sdk import client
from esm.sdk.api import (
    ESMCInferenceClient,
    ESMProtein,
    ESMProteinTensor,
    LogitsConfig,
    LogitsOutput,
)
from esm.sdk.forge import ESM3ForgeInferenceClient
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

# Function to read a single FASTA file
def read_fasta(path, keep_gaps=True, keep_insertions=True, to_upper=False):
    with open(path, "r") as f:
        for result in read_alignment_lines(
            f, keep_gaps=keep_gaps, keep_insertions=keep_insertions, to_upper=to_upper
        ):
            yield result

def read_alignment_lines(lines, keep_gaps=True, keep_insertions=True, to_upper=False):
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


def main_esmc(client: ESMCInferenceClient, protein , args):


    with torch.amp.autocast(device_type="cpu", enabled=args.bf16):

        # Use logits endpoint. Using bf16 for inference optimization
        protein_tensor = client.encode(protein)
        assert isinstance(protein_tensor, ESMProteinTensor), f"Expected ESMProteinTensor but got error: {protein_tensor}"

        output = client.logits(
            protein_tensor,
            LogitsConfig(sequence=True, 
                            return_embeddings=True, 
                            return_hidden_states=True,
                            structure=args.structure,
                            secondary_structure=args.secondary_structure,
                            sasa=args.sasa,
                            function=args.function,
                            residue_annotations=args.residue_annotations,
                            ith_hidden_layer=args.ith_hidden_layer),
        )
        assert isinstance(output, LogitsOutput), f"LogitsOutput was expected but got error: {output}"
        assert output.logits is not None and output.logits.sequence is not None
        assert output.embeddings is not None
        assert output.hidden_states is not None
        print(f"Client returned logits with shape: {output.logits.sequence.shape}, embeddings with shape: {output.embeddings.shape}, and hidden states with shape {output.hidden_states.shape}")

        

        # Request a specific hidden layer
        output_layer = client.logits(
            protein_tensor,
            LogitsConfig(return_hidden_states=True, 
                            ith_hidden_layer=1,
                            sequence=args.sequence,
                            structure=args.structure,
                            secondary_structure=args.secondary_structure,
                            sasa=args.sasa,
                            function=args.function,
                            residue_annotations=args.residue_annotations,
                            return_embeddings=args.return_embeddings))
        assert isinstance(output_layer, LogitsOutput), f"LogitsOutput was expected but got error: {output_layer}"
        assert output_layer.hidden_states is not None
        print(f"Client returned hidden states with shape {output_layer.hidden_states.shape}")


        return output


def processing_fasta(client: ESMCInferenceClient, fasta_file: str, output_dir: str, args):
    """Runs protein folding and saves the output as a PDB file in the specified directory."""
    print(f"Processing {fasta_file}...")

    if args.protein_complex:
        protein_chains = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            protein_chain = make_protein_chain(record)
            protein_chains.append(protein_chain)

        protein = ProteinComplex.from_chains(protein_chains)
        protein = ESMProtein.from_protein_complex(protein)
        logits_output = main_esmc(client, protein, args)
        fasta_name = os.path.splitext(os.path.basename(fasta_file))[0]
        torch.save(logits_output, os.path.join(output_dir, f"{fasta_name}.pt"))
    else:
        fasta_entries = sorted(read_fasta(fasta_file), key=lambda header_seq: len(header_seq[1]))
        for entry in fasta_entries:
            header, sequence = entry
            protein = ESMProtein(sequence=sequence)
            logits_output = main_esmc(client, protein, args)
            fasta_name = os.path.splitext(os.path.basename(header))[0]
            torch.save(logits_output, os.path.join(output_dir, f"{fasta_name}.pt"))

def main(fasta_file: str, output_dir: str, model_name: str, args):

    if not os.path.exists(args.fasta_file) or not args.fasta_file.endswith(".fasta"):
        print(f"Invalid Fasta file: {args.fasta_file}")
        return
    os.makedirs(args.output_dir, exist_ok=True)
    print("Model Name of ESMC:", model_name)

    # List of valid model names
    valid_model_names = ["esmc_600m", "esmc_300m"]

    if model_name not in valid_model_names:
        print(f"Error: '{model_name}' is not a valid model name. Choose from {valid_model_names}.")
        exit(1)

    model = ESMC.from_pretrained(model_name, bf16=args.bf16)

    if args.timing:
        infer_time = time.time()

    processing_fasta(model, fasta_file, output_dir, args)

    if args.timing:
        print(f"Inference time = {time.time() - infer_time} seconds")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ESM3 protein folding for a single FASTA file.")
    parser.add_argument("fasta_file", type=str, help="Path to the FASTA file.")
    parser.add_argument("output_dir", type=str, help="Directory to save the output PDB files.")
    parser.add_argument("--model_name", type=str, default="esmc_600m", help="Name of the model for ESMC.")
    parser.add_argument("--sequence", action="store_true", help="Enable sequence logits")
    parser.add_argument("--structure", action="store_true", help="Enable structure logits")
    parser.add_argument("--secondary_structure", action="store_true", help="Enable secondary structure logits")
    parser.add_argument("--sasa", action="store_true", help="Enable SASA logits")
    parser.add_argument("--function", action="store_true", help="Enable function logits")
    parser.add_argument("--residue_annotations", action="store_true", help="Enable residue annotations")
    parser.add_argument("--return_embeddings", action="store_true", help="Return embeddings")
    parser.add_argument("--return_hidden_states", action="store_true", help="Return hidden states")
    parser.add_argument("--ith_hidden_layer", type=int, default=-1, help="Specify which hidden layer to return")
    parser.add_argument("--bf16", action="store_true", help="Enable bf16 inference.")
    parser.add_argument("--timing", action="store_true", help="Enable timing for inference.")
    parser.add_argument("--protein_complex", action="store_true", help="Enable prediction for protein complexes (multi-chain structure) using a multi-chain FASTA file input.")   


    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    if args.timing:
        start_time = time.time()
    main(args.fasta_file, args.output_dir, args.model_name, args)
    if args.timing:
        print(f"Complete run time = {time.time() - start_time} seconds")
