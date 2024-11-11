#!/usr/bin/env python
# encoding: utf-8

import os
from subprocess import run
from argparse import ArgumentParser

def main(argv):
    # Argument parsing
    parser = ArgumentParser(description="Run ProteinMPNN pipeline for monomers")
    parser.add_argument('--input', help="Input data directory", default="ProteinMPNN/inputs/PDB_monomers/pdbs/")
    parser.add_argument('--output', help="Output directory", default="/outputs/example_7_outputs")
    parser.add_argument('--num_seq_per_target', type=int, default=1, help="Number of sequences per target")
    parser.add_argument('--sampling_temp', type=float, default=0.1, help="Sampling temperature")
    parser.add_argument('--unconditional_probs_only', type=int, default=1, help="Use unconditional probabilities only")
    parser.add_argument('--seed', type=int, default=37, help="Random seed")
    parser.add_argument('--batch_size', type=int, default=1, help="Batch size")
    parser.add_argument('--precision', choices=['float32', 'bfloat16'], default='float32', help="Precision type for calculations")
    parser.add_argument('--use_ipex', action='store_true', help="Enable IPEX optimizations")  # IPEX flag
    args = parser.parse_args()

    # Check and create output directory if it doesn't exist
    output_dir = args.output
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Paths for output files
    folder_with_pdbs = args.input
    path_for_parsed_chains = os.path.join(output_dir, "parsed_pdbs.jsonl")

    # Run the parsing script
    a = run([
        'python', '../helper_scripts/parse_multiple_chains.py',
        '--input_path', folder_with_pdbs,
        '--output_path', path_for_parsed_chains
    ])
    assert a.returncode == 0, "Error parsing multiple chains"

    # Run the main ProteinMPNN script
    a = run([
        'python', '../protein_mpnn_run.py',
        '--jsonl_path', path_for_parsed_chains,
        '--out_folder', output_dir,
        '--num_seq_per_target', str(args.num_seq_per_target),
        '--sampling_temp', str(args.sampling_temp),
        '--unconditional_probs_only', str(args.unconditional_probs_only),
        '--seed', str(args.seed),
        '--batch_size', str(args.batch_size),
        '--precision', args.precision
    ]+ (['--use_ipex'] if args.use_ipex else []))
    assert a.returncode == 0, "Error running protein folding script"

if __name__ == "__main__":
    import sys
    main(sys.argv)

