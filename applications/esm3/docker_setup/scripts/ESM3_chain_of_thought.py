# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

import os
import torch
import argparse
import time
import pandas as pd
from esm.models.esm3 import ESM3
from esm.sdk.api import ESM3InferenceClient, ESMProtein, ESMProteinTensor
from esm.sdk.api import GenerationConfig
from esm.utils.types import FunctionAnnotation

def get_sample_protein_from_csv(csv_file: str, sequence_length: int = None) -> ESMProtein:
    df = pd.read_csv(csv_file, sep="\t")
    if not {"Label", "Start", "End"}.issubset(df.columns):
        raise ValueError("CSV file must contain 'Label', 'Start', and 'End' columns")

    function_annotations = [
        FunctionAnnotation(label=row["Label"], start=int(row["Start"]), end=int(row["End"]))
        for _, row in df.iterrows()
    ]

    # Determine sequence length
    sequence_length = sequence_length or max(df["End"], default=100)

    protein = ESMProtein(sequence="_" * sequence_length)
    protein.function_annotations = function_annotations
    return protein

def chain_of_thought(client: ESM3InferenceClient, csv_path: str, args,output_dir: str = None):
    cot_protein = get_sample_protein_from_csv(csv_path, args.sequence_length)
    enable_autocast = args.bf16
    device_type = "cpu"

    with torch.amp.autocast(device_type=device_type, enabled=enable_autocast):
        cot_protein.sequence = "_" * len(cot_protein.sequence)
        cot_protein.coordinates = None
        cot_protein.sasa = None
        cot_protein_tensor = client.encode(cot_protein)

        # Generate different properties using command-line args
        for cot_track in ["secondary_structure", "structure", "sequence"]:
            cot_protein_tensor = client.generate(
                cot_protein_tensor,
                GenerationConfig(
                    track=cot_track,
                    schedule=args.schedule,
                    num_steps=args.num_steps,
                    strategy=args.strategy,
                    temperature=args.temperature,
                    top_p=args.top_p,
                    condition_on_coordinates_only=args.condition_on_coordinates_only
                ),
            )

        assert isinstance(cot_protein_tensor, ESMProteinTensor)
        cot_protein = client.decode(cot_protein_tensor)
        assert isinstance(cot_protein, ESMProtein)
        if output_dir:
            csv_name = os.path.splitext(os.path.basename(csv_path))[0]
            output_pdb_path = os.path.join(output_dir, f"{csv_name}.pdb")
            cot_protein.to_pdb(output_pdb_path)
            print(f"Saved output to {output_pdb_path}")
        return cot_protein.to_pdb_string()

def main(args):

    if not os.path.exists(args.csv_file) or not args.csv_file.endswith(".csv"):
        print(f"Invalid CSV file: {args.csv_file}")
        return
    os.makedirs(args.output_dir, exist_ok=True)

    client = ESM3InferenceClient() if os.getenv("ESM_API_KEY") else ESM3.from_pretrained(
        "esm3_sm_open_v1", bf16=args.bf16
    )

    infer_time = time.time()
    chain_of_thought(client, args.csv_file, args,args.output_dir)
    if args.timing:
        print(f"Inference time = {time.time() - infer_time} seconds")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run ESM3 protein inverse folding on a CSV file.")
    parser.add_argument("csv_file", type=str, help="Path to input CSV file.")
    parser.add_argument("output_dir", type=str, help="Directory to save the output PDB file.")
    parser.add_argument("--bf16", action="store_true", help="Enable bf16 inference.")
    parser.add_argument("--timing", action="store_true", help="Enable timing for inference.")
    parser.add_argument("--schedule", type=str, choices=["cosine", "linear"], default="cosine", help="Schedule type.")
    parser.add_argument("--strategy", type=str, choices=["random", "entropy"], default="entropy", help="Unmasking strategy.")
    parser.add_argument("--num_steps", type=int, default=1, help="Number of steps for generation.")
    parser.add_argument("--temperature", type=float, default=1.0, help="Sampling temperature.")
    parser.add_argument("--top_p", type=float, default=1.0, help="Top-p sampling value.")
    parser.add_argument("--condition_on_coordinates_only", action="store_true", help="Condition only on coordinates.")
    parser.add_argument("--sequence_length", type=int,default=None, help="Custom sequence length (optional).")

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    if args.timing:
        start_time = time.time()
    main(args)
    if args.timing:
        print(f"Complete run time = {time.time() - start_time:.2f} seconds")
