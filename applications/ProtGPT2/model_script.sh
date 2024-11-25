#!/bin/bash

# Directory where the files will be downloaded
model_dir="./model_dir"

# Array of URLs to download
URLS=(
    "https://huggingface.co/nferruz/ProtGPT2/resolve/main/.gitattributes"
	"https://huggingface.co/nferruz/ProtGPT2/resolve/main/config.json"
	"https://huggingface.co/nferruz/ProtGPT2/resolve/main/merges.txt"
	"https://huggingface.co/nferruz/ProtGPT2/resolve/main/ppl-plddt.png"
	"https://huggingface.co/nferruz/ProtGPT2/resolve/main/pytorch_model.bin"
	"https://huggingface.co/nferruz/ProtGPT2/resolve/main/special_tokens_map.json"
	"https://huggingface.co/nferruz/ProtGPT2/resolve/main/tokenizer.json"
	"https://huggingface.co/nferruz/ProtGPT2/resolve/main/vocab.json"

)

# Create the directory if it doesn't exist
mkdir -p "$model_dir"

# Change to the download directory
cd "$model_dir" || exit

# Loop through each URL and download the file if it doesn't exist
for url in "${URLS[@]}"; do
    filename=$(basename "$url")
    if [ -f "$filename" ]; then
        echo "$filename already exists. Skipping download."
    else
        echo "Downloading $filename..."
        wget "$url"
    fi
done

echo "Download process completed."

