#!/bin/bash

set -e
SCRIPT_PATH="${BASH_SOURCE:-$0}"
ABS_SCRIPT_PATH="$(realpath "${SCRIPT_PATH}")"
#echo "Value of ABS_SCRIPT_PATH: ${ABS_SCRIPT_PATH}"
ABS_DIRECTORY="$(dirname "${ABS_SCRIPT_PATH}")"
# Default Conda installation directory
CONDA_INSTALL_DIR=$(realpath ./miniforge3)

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
# Check if Miniforge3 exists and install if not found
if [ ! -d "$CONDA_INSTALL_DIR" ]; then
  echo "Miniforge3 is not installed. Installing..."
  wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
  bash Miniforge3-Linux-x86_64.sh -b -p "$CONDA_INSTALL_DIR"
  echo "Miniforge3 installation complete."
else
  echo "Miniforge3 is already installed at: $CONDA_INSTALL_DIR"
fi
# Export Conda binary path
export PATH="$CONDA_INSTALL_DIR/bin:$PATH"
# Clone the RFdiffusion repository if it doesn't exist
if [ ! -d "RFdiffusion" ]; then
  git clone https://github.com/RosettaCommons/RFdiffusion.git
else
  echo "RFdiffusion repository already exists, skipping git clone."
fi

echo "$CONDA_INSTALL_DIR"
# Apply patch (assuming patch file is RFdiffusion.patch and it should be applied in RFdiffusion directory)
cd RFdiffusion
git checkout 820bfdfaded8c260b962dc40a3171eae316b6ce0
git log -1
PATCH_FILE="$ABS_DIRECTORY/RFdiffusion.patch"
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
mkdir -p models
cd models
wget https://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt
wget https://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt
wget https://files.ipd.uw.edu/pub/RFdiffusion/60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt
wget https://files.ipd.uw.edu/pub/RFdiffusion/74f51cfb8b440f50d70878e05361d8f0/InpaintSeq_ckpt.pt
wget https://files.ipd.uw.edu/pub/RFdiffusion/76d00716416567174cdb7ca96e208296/InpaintSeq_Fold_ckpt.pt
wget https://files.ipd.uw.edu/pub/RFdiffusion/5532d2e1f3a4738decd58b19d633b3c3/ActiveSite_ckpt.pt
wget https://files.ipd.uw.edu/pub/RFdiffusion/12fc204edeae5b57713c5ad7dcb97d39/Base_epoch8_ckpt.pt
# Optional:
wget https://files.ipd.uw.edu/pub/RFdiffusion/f572d396fae9206628714fb2ce00f72e/Complex_beta_ckpt.pt
# original structure prediction weights
wget https://files.ipd.uw.edu/pub/RFdiffusion/1befcb9b28e2f778f53d47f18b7597fa/RF_structure_prediction_weights.pt
cd ../
# Create and activate the Conda environment using the YAML file, disabling plugins to avoid errors
#CONDA_NO_PLUGINS=true 
conda env create -f env/SE3nv.yml 
source $CONDA_INSTALL_DIR/bin/activate SE3nv
#conda init
#conda activate SE3nv

# Install SE3Transformer requirements
cd env/SE3Transformer
pip install --no-cache-dir -r requirements.txt
python setup.py install

# Install the rfdiffusion module
cd ../.. # Change into the root directory of the repository
pip install -e .

echo "Conda (Miniforge3) is installed at $CONDA_INSTALL_DIR"
echo "To manually activate conda env do: source $CONDA_INSTALL_DIR/bin/activate SE3nv"
