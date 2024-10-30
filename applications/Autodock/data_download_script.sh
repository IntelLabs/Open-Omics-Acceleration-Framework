#!/bin/bash

# Define variables
url="https://zenodo.org/records/4031961/files/data.zip?download=1"
download_dir="./data"
target_folder="4fev"

# Create a temporary data directory
mkdir -p "$download_dir"

# Download the data.zip file
echo "Downloading data.zip..."
wget -O "$download_dir/data.zip" "$url"

# Unzip the downloaded file
echo "Unzipping data.zip..."
unzip "$download_dir/data.zip" -d "$download_dir"

# Check for nested 'data' folder and move target folder to current directory
if [ -d "$download_dir/data" ]; then
  echo "Nested 'data' folder found."
  
  # Move only the target folder to the current directory
  if [ -d "$download_dir/data/$target_folder" ]; then
    mv "$download_dir/data/$target_folder" ./
    echo "$target_folder folder successfully moved to the current directory."
  else
    echo "$target_folder folder not found inside nested 'data' directory."
  fi

  # Remove the intermediate data directory
  rm -rf "$download_dir/data"
fi

# Clean up zip file and data directory
rm -rf "$download_dir"
echo "Download and extraction complete. '$target_folder' folder is now in the current directory."

