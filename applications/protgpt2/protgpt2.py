#!/usr/bin/env python
# encoding: utf-8

from transformers import pipeline


protgpt2 = pipeline('text-generation', model="nferruz/ProtGPT2")

sequences = protgpt2("<|endoftext|>", max_length=100, do_sample=True, top_k=950, repetition_penalty=1.2, num_return_sequences=10, eos_token_id=0)

for seq in sequences:
	print(seq):
