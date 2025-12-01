# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License
# Author: Rahamathullah
# Email shaikx.rahamathullah@intel.com

#!/bin/bash

set -e
trap 'echo "Error on line $LINENO"; exit 1;' ERR

SCRIPT_PATH="${BASH_SOURCE:-$0}"
ABS_SCRIPT_PATH="$(realpath "${SCRIPT_PATH}")"
ABS_DIRECTORY="$(dirname "${ABS_SCRIPT_PATH}")"

CONDA_INSTALL_DIR=$(realpath ./miniforge3)

# Parse command line arguments
while (( "$#" )); do
  case "$1" in
    -p)
      CONDA_INSTALL_DIR=$2
      CONDA_INSTALL_DIR=$(realpath "$CONDA_INSTALL_DIR")
      shift 2
      ;;
    -*|--*=)
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    *)
      echo "Error: Unsupported argument $1" >&2
      exit 1
      ;;
  esac
done

# Check if Miniforge3 exists
if [ ! -d "$CONDA_INSTALL_DIR" ]; then
  echo "Miniforge3 is not installed. Installing..."
  command -v wget >/dev/null 2>&1 || { echo "wget is required but not installed. Exiting."; exit 1; }
  wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
  bash Miniforge3-Linux-x86_64.sh -b -p "$CONDA_INSTALL_DIR"
  echo "Miniforge3 installation complete."
else
  echo "Miniforge3 is already installed at: $CONDA_INSTALL_DIR"
fi

export PATH="$CONDA_INSTALL_DIR/bin:$PATH"

# Clone the ProteinMPNN repository
if [ ! -d "ProteinMPNN" ]; then
  git clone https://github.com/dauparas/ProteinMPNN.git
else
  echo "ProteinMPNN repository already exists, skipping git clone."
fi

cd ProteinMPNN
git checkout 8907e6671bfbfc92303b5f79c4b5e6ce47cdef57
PATCH_FILE="$ABS_DIRECTORY/ProteinMPNN_torch.patch"
if [ -f "$PATCH_FILE" ]; then
  if git apply --reverse --check "$PATCH_FILE" > /dev/null 2>&1; then
    echo "Patch has already been applied. Skipping patch step."
  else
    git apply "$PATCH_FILE"
    echo "Patch applied successfully."
  fi
else
  echo "Error: Patch file not found at $PATCH_FILE" >&2
  exit 1
fi

# Create and activate the Conda environment
#source "$CONDA_INSTALL_DIR/bin/activate"
if conda env list | grep -q "^p_mpnn"; then
	echo "Environment exists. Moving ahead without create the env. If the setup crashes, please remove manually."
    else
	echo "Creating conda env p_mpnn.."
	conda create -n p_mpnn -y python=3.11 pip=24.0
fi

source $CONDA_INSTALL_DIR/bin/activate p_mpnn
#conda activate p_mpnn

conda install -n p_mpnn -y pytorch==2.3.1 torchvision==0.18.1 torchaudio==2.3.1 -c pytorch
pip install numpy==1.26.0

echo "setup complete!"
echo "Note:"
echo "Conda (Miniforge3) is installed at $CONDA_INSTALL_DIR"
echo "To manually activate conda env, do: source $CONDA_INSTALL_DIR/bin/activate p_mpnn"
