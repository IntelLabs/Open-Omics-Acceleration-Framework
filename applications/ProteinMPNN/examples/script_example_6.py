#!/usr/bin/env python
# encoding: utf-8

import os
from subprocess import run
from argparse import ArgumentParser

def main(argv):
    # Argument parsing
    parser = ArgumentParser(description="Run ProteinMPNN pipeline with homooligomers and tied positions")
    parser.add_argument('--input', help="Input data directory", default="/ProteinMPNN/inputs/PDB_homooligomers/pdbs/")
    parser.add_argument('--output', help="Output directory", default="/outputs/example_6_outputs")
    parser.add_argument('--homooligomer', type=int, default=1, help="Homooligomer flag")
    parser.add_argument('--num_seq_per_target', type=int, default=2, help="Number of sequences per target")
    parser.add_argument('--sampling_temp', type=float, default=0.2, help="Sampling temperature")
    parser.add_argument('--seed', type=int, default=37, help="Random seed")
    parser.add_argument('--batch_size', type=int, default=1, help="Batch size")
    parser.add_argument('--precision', choices=['float32', 'bfloat16'], default='float32', help="Precision type for calculations")
    parser.add_argument('--use_ipex', action='store_true', help="Enable IPEX optimizations")
    args = parser.parse_args()

    # Check and create output directory if it doesn't exist
    output_dir = args.output
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Paths for output files
    folder_with_pdbs = args.input
    path_for_parsed_chains = os.path.join(output_dir, "parsed_pdbs.jsonl")
    path_for_tied_positions = os.path.join(output_dir, "tied_pdbs.jsonl")
    path_for_designed_sequences = os.path.join(output_dir, "temp_0.1")

    # Run the parsing script
    a = run([
        'python', 'helper_scripts/parse_multiple_chains.py',
        '--input_path', folder_with_pdbs,
        '--output_path', path_for_parsed_chains
    ])
    assert a.returncode == 0, "Error parsing multiple chains"

    # Run the tied positions script
    a = run([
        'python', 'helper_scripts/make_tied_positions_dict.py',
        '--input_path', path_for_parsed_chains,
        '--output_path', path_for_tied_positions,
        '--homooligomer', str(args.homooligomer)
    ])
    assert a.returncode == 0, "Error making tied positions dict"

    # Run the main ProteinMPNN script
    a = run([
        'python', 'protein_mpnn_run.py',
        '--jsonl_path', path_for_parsed_chains,
        '--tied_positions_jsonl', path_for_tied_positions,
        '--out_folder', output_dir,
        '--num_seq_per_target', str(args.num_seq_per_target),
        '--sampling_temp', str(args.sampling_temp),
        '--seed', str(args.seed),
        '--batch_size', str(args.batch_size),
        '--precision', args.precision
    ]+ (['--use_ipex'] if args.use_ipex else []))
    assert a.returncode == 0, "Error running protein folding script"

if __name__ == "__main__":
    import sys
    main(sys.argv)
