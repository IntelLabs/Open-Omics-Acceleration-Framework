#!/bin/bash

set -e
SCRIPT_PATH="${BASH_SOURCE:-$0}"
ABS_SCRIPT_PATH="$(realpath "${SCRIPT_PATH}")"
#echo "Value of ABS_SCRIPT_PATH: ${ABS_SCRIPT_PATH}"
ABS_DIRECTORY="$(dirname "${ABS_SCRIPT_PATH}")"
# Default Conda installation directory
CONDA_INSTALL_DIR=$(realpath /home/miniforge3)

# Parse command line arguments
while (( "$#" )); do
  case "$1" in
    -p)
      CONDA_INSTALL_DIR=$2
      CONDA_INSTALL_DIR=$(realpath "$CONDA_INSTALL_DIR")
      shift 2
      ;;
    -*|--*=) # Unsupported flags
      echo "Error: Unsupported flag $1" >&2
      exit 1
      ;;
    *) # Preserve positional arguments
      echo "Error: Unsupported argument $1" >&2
      exit 1
      ;;
  esac
done

# Clone the RFdiffusion repository if it doesn't exist
if [ ! -d "RFdiffusion" ]; then
  git clone https://github.com/RosettaCommons/RFdiffusion.git
else
  echo "RFdiffusion repository already exists, skipping git clone."
fi

echo "$CONDA_INSTALL_DIR"
# Apply patch (assuming patch file is new_changed.patch and it should be applied in RFdiffusion directory)
cd RFdiffusion
#PATCH_FILE="/home/hgx/omics/hgx/rf_diff/setup_rfdiffusion/new_changes.patch"
PATCH_FILE="$ABS_DIRECTORY/new_changes.patch"
echo $PATCH_FILE
if [ -f "$PATCH_FILE" ]; then
  # Check if the patch is already applied
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

# Create and activate the Conda environment using the YAML file, disabling plugins to avoid errors
CONDA_NO_PLUGINS=true conda env create -f env/SE3nv.yml 

#source setup_conda.sh
conda activate SE3nv


# Install SE3Transformer requirements
cd env/SE3Transformer
pip install --no-cache-dir -r requirements.txt
python setup.py install

# Install the rfdiffusion module
cd ../.. # Change into the root directory of the repository
pip install -e .


