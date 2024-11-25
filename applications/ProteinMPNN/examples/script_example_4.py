#!/usr/bin/env python
# encoding: utf-8
import os
from subprocess import run
from argparse import ArgumentParser

def main(argv):
    # Argument parsing
    parser = ArgumentParser(description="Run protein folding with fixed positions on specific chains")
    parser.add_argument('--input', help="Input folder with PDBs", default="ProteinMPNN/inputs/PDB_complexes/pdbs/")
    parser.add_argument('--output', help="Output directory", default="/outputs/example_4_outputs")
    parser.add_argument('--chains_to_design', help="Chains to design (space-separated)", default="A C")
    parser.add_argument('--fixed_positions', help="Fixed positions for the chains", default="1 2 3 4 5 6 7 8 23 25, 10 11 12 13 14 15 16 17 18 19 20 40")
    parser.add_argument('--num_seq_per_target', type=int, default=2, help="Number of sequences per target")
    parser.add_argument('--sampling_temp', type=float, default=0.1, help="Sampling temperature")
    parser.add_argument('--seed', type=int, default=37, help="Random seed")
    parser.add_argument('--batch_size', type=int, default=1, help="Batch size")
    parser.add_argument('--precision', choices=['float32', 'bfloat16'], default='float32', help="Precision type for calculations")
    parser.add_argument('--use_ipex', action='store_true', help="Enable IPEX optimizations")
    args = parser.parse_args()

    # Check and create output directory if it doesn't exist
    output_dir = args.output
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    folder_with_pdbs = args.input
    path_for_parsed_chains = os.path.join(output_dir, "parsed_pdbs.jsonl")
    path_for_assigned_chains = os.path.join(output_dir, "assigned_pdbs.jsonl")
    path_for_fixed_positions = os.path.join(output_dir, "fixed_pdbs.jsonl")

    # Run parsing script
    a = run([
        'python', 'helper_scripts/parse_multiple_chains.py',
        '--input_path', folder_with_pdbs,
        '--output_path', path_for_parsed_chains
    ])
    assert a.returncode == 0, "Error parsing multiple chains"

    # Assign chains to design
    a = run([
        'python', 'helper_scripts/assign_fixed_chains.py',
        '--input_path', path_for_parsed_chains,
        '--output_path', path_for_assigned_chains,
        '--chain_list', args.chains_to_design
    ])
    assert a.returncode == 0, "Error assigning fixed chains"

    # Make fixed positions dictionary
    a = run([
        'python', 'helper_scripts/make_fixed_positions_dict.py',
        '--input_path', path_for_parsed_chains,
        '--output_path', path_for_fixed_positions,
        '--chain_list', args.chains_to_design,
        '--position_list', args.fixed_positions
    ])
    assert a.returncode == 0, "Error making fixed positions dictionary"

    # Run protein folding script
    a = run([
        'python', 'protein_mpnn_run.py',
        '--jsonl_path', path_for_parsed_chains,
        '--chain_id_jsonl', path_for_assigned_chains,
        '--fixed_positions_jsonl', path_for_fixed_positions,
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

