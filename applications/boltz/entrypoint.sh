#!/bin/bash

INPUT_DIR="/app/boltz/input"
OUTPUT_DIR="/app/boltz/output"
CUDA_VISIBLE_DEVICES=""

FOUND=0

for INPUT_FILE in "$INPUT_DIR"/*; do
    # Only process .yaml or .fasta files
    if [[ "$INPUT_FILE" == *.yaml || "$INPUT_FILE" == *.fasta ]]; then
        echo "📂 Processing: $INPUT_FILE"
        LD_PRELOAD=/opt/conda/lib/libjemalloc.so:$LD_PRELOAD \
        MALLOC_CONF="oversize_threshold:1,background_thread:true,metadata_thp:auto,dirty_decay_ms:-1,muzzy_decay_ms:-1" \
        boltz predict "$INPUT_FILE" --out_dir "$OUTPUT_DIR" --accelerator "cpu"
        FOUND=1
    fi
done

if [[ $FOUND -eq 0 ]]; then
    echo "❌ No .yaml or .fasta files found in $INPUT_DIR"
    exit 1
fi
