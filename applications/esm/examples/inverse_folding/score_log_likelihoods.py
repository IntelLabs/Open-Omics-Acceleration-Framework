# Copyright (c) Meta Platforms, Inc. and affiliates.
#
# This source code is licensed under the MIT license found in the
# LICENSE file in the root directory of this source tree.
#
# Scores sequences based on a given structure.
#
# usage:
# score_log_likelihoods.py [-h] [--outpath OUTPATH] [--chain CHAIN] pdbfile seqfile

import argparse
from biotite.sequence.io.fasta import FastaFile, get_sequences
import numpy as np
from pathlib import Path
import torch
import torch.nn.functional as F
from tqdm import tqdm
import time
import esm
import esm.inverse_folding


def score_singlechain_backbone(model, alphabet, args):
    if torch.cuda.is_available() and not args.nogpu:
        model = model.cuda()
        print("Transferred model to GPU")
    coords, native_seq = esm.inverse_folding.util.load_coords(args.pdbfile, args.chain)
    print('Native sequence loaded from structure file:')
    print(native_seq)
    print('\n')

    ll, _ = esm.inverse_folding.util.score_sequence(
            model, alphabet, coords, native_seq)
    print('Native sequence')
    print(f'Log likelihood: {ll:.2f}')
    print(f'Perplexity: {np.exp(-ll):.2f}')

    print('\nScoring variant sequences from sequence file..\n')
    infile = FastaFile()
    infile.read(args.seqfile)
    seqs = get_sequences(infile)
    Path(args.outpath).parent.mkdir(parents=True, exist_ok=True)
    with open(args.outpath, 'w') as fout:
        fout.write('seqid,log_likelihood\n')
        for header, seq in tqdm(seqs.items()):
            ll, _ = esm.inverse_folding.util.score_sequence(
                    model, alphabet, coords, str(seq))
            fout.write(header + ',' + str(ll) + '\n')
    print(f'Results saved to {args.outpath}')


def score_multichain_backbone(model, alphabet, args):
    if torch.cuda.is_available() and not args.nogpu:
        model = model.cuda()
        print("Transferred model to GPU")
    structure = esm.inverse_folding.util.load_structure(args.pdbfile)
    coords, native_seqs = esm.inverse_folding.multichain_util.extract_coords_from_complex(structure)
    target_chain_id = args.chain
    native_seq = native_seqs[target_chain_id]
    print('Native sequence loaded from structure file:')
    print(native_seq)
    print('\n')

    ll, _ = esm.inverse_folding.multichain_util.score_sequence_in_complex(
            model, alphabet, coords, target_chain_id, native_seq)
    print('Native sequence')
    print(f'Log likelihood: {ll:.2f}')
    print(f'Perplexity: {np.exp(-ll):.2f}')

    print('\nScoring variant sequences from sequence file..\n')
    infile = FastaFile()
    infile.read(args.seqfile)
    seqs = get_sequences(infile)
    Path(args.outpath).parent.mkdir(parents=True, exist_ok=True)
    with open(args.outpath, 'w') as fout:
        fout.write('seqid,log_likelihood\n')
        for header, seq in tqdm(seqs.items()):
            ll, _ = esm.inverse_folding.multichain_util.score_sequence_in_complex(
                    model, alphabet, coords, target_chain_id, str(seq))
            fout.write(header + ',' + str(ll) + '\n')
    print(f'Results saved to {args.outpath}')


def main():
    parser = argparse.ArgumentParser(
            description='Score sequences based on a given structure.'
    )
    parser.add_argument(
            'pdbfile', type=str,
            help='input filepath, either .pdb or .cif',
    )
    parser.add_argument(
            'seqfile', type=str,
            help='input filepath for variant sequences in a .fasta file',
    )
    parser.add_argument(
            '--outpath', type=str,
            help='output filepath for scores of variant sequences',
            default='output/sequence_scores.csv',
    )
    parser.add_argument(
            '--chain', type=str,
            help='chain id for the chain of interest', default='A',
    )
    parser.set_defaults(multichain_backbone=False)
    parser.add_argument(
            '--multichain-backbone', action='store_true',
            help='use the backbones of all chains in the input for conditioning'
    )
    parser.add_argument(
            '--singlechain-backbone', dest='multichain_backbone',
            action='store_false',
            help='use the backbone of only target chain in the input for conditioning'
    )

    parser.add_argument("--nogpu", action="store_true", help="Do not use GPU even if available")
    parser.add_argument("--noipex", action="store_true", help="Do not use intel_extension_for_pytorch")
    parser.add_argument("--bf16", action="store_true", help="Use bf16 precision")
    args = parser.parse_args()

    model, alphabet = esm.pretrained.esm_if1_gvp4_t16_142M_UR50()
    model = model.eval()

    if not args.noipex:
        print("use ipex .......................")
        dtype = torch.bfloat16 if args.bf16 else torch.float32
        print("dtype.............",dtype)
        import intel_extension_for_pytorch as ipex
        model = ipex.optimize(model, dtype=dtype)
    if args.noipex and args.bf16:
        print("direct code with bf16")
        model=model.bfloat16()

    enable_autocast = args.bf16
    p0=time.time()
    with torch.cpu.amp.autocast(enable_autocast):
        if args.multichain_backbone:
            print("enable_autocast.............",enable_autocast)
            score_multichain_backbone(model, alphabet, args)
        else:
            print("enable_autocast.............",enable_autocast)
            score_singlechain_backbone(model, alphabet, args)
    print("infernce time",time.time()-p0)


if __name__ == '__main__':
    main()
