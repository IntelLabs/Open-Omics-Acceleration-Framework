#!/bin/bash

# Directory where the files will be downloaded
models="./models"

# Array of URLs to download
URLS=(
    "https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t33_650M_UR50D.pt"
	"https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t33_650M_UR50D-contact-regression.pt"
	"https://dl.fbaipublicfiles.com/fair-esm/models/esm_if1_gvp4_t16_142M_UR50.pt"
	"https://dl.fbaipublicfiles.com/fair-esm/models/esmfold_3B_v1.pt"
	"https://dl.fbaipublicfiles.com/fair-esm/models/esm2_t36_3B_UR50D.pt"
	"https://dl.fbaipublicfiles.com/fair-esm/regression/esm2_t36_3B_UR50D-contact-regression.pt"
)

# Create the directory if it doesn't exist
mkdir -p "$models"

# Change to the download directory
cd "$models" || exit

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
