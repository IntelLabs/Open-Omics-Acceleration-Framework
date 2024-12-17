#!/usr/bin/env python
# encoding: utf-8

import os
from subprocess import run
from argparse import ArgumentParser

def main(argv):
    # Argument parsing
    parser = ArgumentParser(description="Run ProteinMPNN pipeline with multiple chains and positions")
    parser.add_argument('--input', help="Input data directory", default="/ProteinMPNN/inputs/PDB_complexes/pdbs/")
    parser.add_argument('--output', help="Output directory", default="/outputs/example_5_outputs")
    parser.add_argument('--chains_to_design', default="A C", help="Chains to design")
    parser.add_argument('--fixed_positions', default="9 10 11 12 13 14 15 16 17 18 19 20 21 22 23, 10 11 18 19 20 22", help="Fixed positions for the chains")
    parser.add_argument('--tied_positions', default="1 2 3 4 5 6 7 8, 1 2 3 4 5 6 7 8", help="Tied positions for the chains")
    parser.add_argument('--num_seq_per_target', type=int, default=2, help="Number of sequences per target")
    parser.add_argument('--sampling_temp', type=float, default=0.1, help="Sampling temperature")
    parser.add_argument('--seed', type=int, default=37, help="Random seed")
    parser.add_argument('--batch_size', type=int, default=1, help="Batch size")
    parser.add_argument('--precision', choices=['float32', 'bfloat16'], default='float32', help="Precision type for calculations")
    args = parser.parse_args()

    # Check and create output directory if it doesn't exist
    output_dir = args.output
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Paths for output files
    folder_with_pdbs = args.input
    path_for_parsed_chains = os.path.join(output_dir, "parsed_pdbs.jsonl")
    path_for_assigned_chains = os.path.join(output_dir, "assigned_pdbs.jsonl")
    path_for_fixed_positions = os.path.join(output_dir, "fixed_pdbs.jsonl")
    path_for_tied_positions = os.path.join(output_dir, "tied_pdbs.jsonl")

    # Run the parsing script
    a = run([
        'python', 'helper_scripts/parse_multiple_chains.py',
        '--input_path', folder_with_pdbs,
        '--output_path', path_for_parsed_chains
    ])
    assert a.returncode == 0, "Error parsing multiple chains"

    # Run the chain assignment script
    a = run([
        'python', 'helper_scripts/assign_fixed_chains.py',
        '--input_path', path_for_parsed_chains,
        '--output_path', path_for_assigned_chains,
        '--chain_list', args.chains_to_design
    ])
    assert a.returncode == 0, "Error assigning fixed chains"

    # Run the fixed positions script
    a = run([
        'python', 'helper_scripts/make_fixed_positions_dict.py',
        '--input_path', path_for_parsed_chains,
        '--output_path', path_for_fixed_positions,
        '--chain_list', args.chains_to_design,
        '--position_list', args.fixed_positions
    ])
    assert a.returncode == 0, "Error making fixed positions dict"

    # Run the tied positions script
    a = run([
        'python', 'helper_scripts/make_tied_positions_dict.py',
        '--input_path', path_for_parsed_chains,
        '--output_path', path_for_tied_positions,
        '--chain_list', args.chains_to_design,
        '--position_list', args.tied_positions
    ])
    assert a.returncode == 0, "Error making tied positions dict"

    # Run the main ProteinMPNN script
    a = run([
        'python', 'protein_mpnn_run.py',
        '--jsonl_path', path_for_parsed_chains,
        '--chain_id_jsonl', path_for_assigned_chains,
        '--fixed_positions_jsonl', path_for_fixed_positions,
        '--tied_positions_jsonl', path_for_tied_positions,
        '--out_folder', output_dir,
        '--num_seq_per_target', str(args.num_seq_per_target),
        '--sampling_temp', str(args.sampling_temp),
        '--seed', str(args.seed),
        '--batch_size', str(args.batch_size),
        '--precision', args.precision
    ])
    assert a.returncode == 0, "Error running protein folding script"

if __name__ == "__main__":
    import sys
    main(sys.argv)
