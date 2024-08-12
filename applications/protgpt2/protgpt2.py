#!/usr/bin/env python
# encoding: utf-8
import os
import argparse
import time
import torch
import random
import numpy as np
import intel_extension_for_pytorch as ipex
from transformers import pipeline

def make_deterministic(seed=42):
    torch.manual_seed(seed)
    np.random.seed(seed)
    random.seed(seed)


def main():
    # Setting up argument parser
    parser = argparse.ArgumentParser(description='Generate protein sequences using ProtGPT2 with IPEX optimization.')
    parser.add_argument('--max_length', type=int, default=100, help='Maximum length of generated sequence')
    parser.add_argument('--do_sample', type=bool, default=True, help='Whether to sample the output or not')
    parser.add_argument('--top_k', type=int, default=950, help='The number of highest probability vocabulary tokens to keep for top-k-filtering')
    parser.add_argument('--repetition_penalty', type=float, default=1.2, help='The parameter for repetition penalty. 1.0 means no penalty')
    parser.add_argument('--num_return_sequences', type=int, default=1, help='The number of sequences to return')
    parser.add_argument('--eos_token_id', type=int, default=0, help='The id of the end of sequence token')
    parser.add_argument('--dtype', type=str, choices=['float32', 'bfloat16'], default='float32', help='Data type for model optimization')
    parser.add_argument('--iterations', type=int, default=5, help='Number of iterations to run')
    parser.add_argument('--model_dir', type=str, required=True, help='Directory to save or load the model')
    args = parser.parse_args()


    make_deterministic()
    # Setting dtype
    dtype = torch.float32 if args.dtype == 'float32' else torch.bfloat16

    model_dir = args.model_dir
    # Generate sequences using ProtGPT2 with IPEX optimization
    protgpt2 = pipeline('text-generation', model=model_dir, torch_dtype=dtype)
    protgpt2.model = ipex.optimize(protgpt2.model, dtype=dtype)

    tic = time.time()

    for i in range(args.iterations):
        print("Iteration:", i)
        sequences = protgpt2(
            "<|endoftext|>",
            max_length=args.max_length,
            do_sample=args.do_sample,
            top_k=args.top_k,
            repetition_penalty=args.repetition_penalty,
            num_return_sequences=args.num_return_sequences,
            eos_token_id=args.eos_token_id
        )

    #output_path = "/app/output/output_file.txt"

    #with open(output_path, 'w') as f:
    #    for seq in sequences:
    #        f.write(seq['generated_text'] + '\n')

    toc = time.time()
    print('Time taken for', args.iterations, 'iterations:', toc - tic, 'seconds')
    print('Average time per iteration:', (toc - tic) / args.iterations, 'seconds')

    # Printing the first two sequences after all iterations
	for seq in sequences[0:2]:
        print(seq)

if __name__ == "__main__":
    main()

