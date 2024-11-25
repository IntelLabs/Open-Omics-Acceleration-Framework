#!/usr/bin/env python
# encoding: utf-8

import os
from subprocess import run
from argparse import ArgumentParser

def main(argv):
    # Argument parsing
    parser = ArgumentParser(description="Run ProteinMPNN pipeline with PSSM inputs")
    parser.add_argument('--pssm_input', help="PSSM input directory", default="ProteinMPNN/inputs/PSSM_inputs")
    parser.add_argument('--input', help="Input data directory", default="ProteinMPNN/inputs/PDB_complexes/pdbs/")
    parser.add_argument('--output', help="Output directory", default="/outputs/example_pssm_outputs")
    parser.add_argument('--chains_to_design', default="A B", help="Chains to design")
    parser.add_argument('--num_seq_per_target', type=int, default=2, help="Number of sequences per target")
    parser.add_argument('--sampling_temp', type=float, default=0.1, help="Sampling temperature")
    parser.add_argument('--seed', type=int, default=37, help="Random seed")
    parser.add_argument('--batch_size', type=int, default=1, help="Batch size")
    parser.add_argument('--pssm_multi', type=float, default=0.3, help="PSSM multiplier")
    parser.add_argument('--pssm_bias_flag', type=int, default=1, help="PSSM bias flag")
    parser.add_argument('--precision', choices=['float32', 'bfloat16'], default='float32', help="Precision type for calculations")
    parser.add_argument('--use_ipex', action='store_true', help="Enable IPEX optimizations")
    args = parser.parse_args()

    # Check and create output directory if it doesn't exist
    output_dir = args.output
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Paths for output files
    folder_with_pdbs = args.input
    pssm_input_path = args.pssm_input
    path_for_parsed_chains = os.path.join(output_dir, "parsed_pdbs.jsonl")
    path_for_assigned_chains = os.path.join(output_dir, "assigned_pdbs.jsonl")
    pssm_output_path = os.path.join(output_dir, "pssm.jsonl")

    # Run parsing script
    a = run([
        'python', 'helper_scripts/parse_multiple_chains.py',
        '--input_path', folder_with_pdbs,
        '--output_path', path_for_parsed_chains
    ])
    assert a.returncode == 0, "Error parsing multiple chains"

    # Run assign fixed chains script
    a = run([
        'python', 'helper_scripts/assign_fixed_chains.py',
        '--input_path', path_for_parsed_chains,
        '--output_path', path_for_assigned_chains,
        '--chain_list', args.chains_to_design
    ])
    assert a.returncode == 0, "Error assigning fixed chains"

    # Run make PSSM input dictionary script
    a = run([
        'python', 'helper_scripts/make_pssm_input_dict.py',
        '--jsonl_input_path', path_for_parsed_chains,
        '--PSSM_input_path', pssm_input_path,
        '--output_path', pssm_output_path
    ])
    assert a.returncode == 0, "Error creating PSSM input dictionary"

    # Run the main ProteinMPNN script with PSSM
    a = run([
        'python', 'protein_mpnn_run.py',
        '--jsonl_path', path_for_parsed_chains,
        '--chain_id_jsonl', path_for_assigned_chains,
        '--out_folder', output_dir,
        '--num_seq_per_target', str(args.num_seq_per_target),
        '--sampling_temp', str(args.sampling_temp),
        '--seed', str(args.seed),
        '--batch_size', str(args.batch_size),
        '--pssm_jsonl', pssm_output_path,
        '--pssm_multi', str(args.pssm_multi),
        '--pssm_bias_flag', str(args.pssm_bias_flag),
        '--precision', args.precision
    ]+ (['--use_ipex'] if args.use_ipex else []))
    assert a.returncode == 0, "Error running protein folding script"

if __name__ == "__main__":
    import sys
    main(sys.argv)

