#!/bin/bash

# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License
# Author: Rahamathullah
# Email shaikx.rahamathullah@intel.com
set -e  # Stop on first error

# -------------------------------------
# Default values for all possible args
# -------------------------------------
PDB_FILE=""
CONTIG=""
HOTSPOT_RES=""
NUM_DESIGNS=1
PRECISION="bfloat16"
PARTIAL_T=""
NOISE_SCALE_CA=""
NOISE_SCALE_FRAME=""
OUTPUT_DIR="./output_dir"
DATA_DIR=""

# -------------------------
# Parse CLI arguments
# -------------------------
with_pdb_db() {
while [[ "$#" -gt 0 ]]; do
  case "$1" in
    --input_pdb) PDB_FILE="$2"; shift ;;
    --contig) CONTIG="$2"; shift ;;
    --hotspot_res) HOTSPOT_RES="$2"; shift ;;
    --num_designs) NUM_DESIGNS="$2"; shift ;;
    --precision) PRECISION="$2"; shift ;;
    --partial_T) PARTIAL_T="$2"; shift ;;
    --noise_scale_ca) NOISE_SCALE_CA="$2"; shift ;;
    --noise_scale_frame) NOISE_SCALE_FRAME="$2"; shift ;;
    --output_dir) OUTPUT_DIR="$2"; shift ;;
    --db_af2_path) DATA_DIR="$2"; shift ;;
    *) echo "❌ Unknown argument: $1"; exit 1 ;;
  esac
  shift
done


if [[ -z "$CONTIG" ]]; then
  echo "❌ CONTIG is required. Please provide it using --contig"
  exit 1
fi

if [[ "$PRECISION" != "float32" && "$PRECISION" != "bfloat16" ]]; then
  echo "❌ Invalid precision: '$PRECISION'. Only 'float32' or 'bfloat16' are supported."
  exit 1
fi

if [[ -z "$DATA_DIR" || ! -d "$DATA_DIR" ]]; then
  echo "❌ DATA_DIR not provided or doesn't exist: $DATA_DIR"
  exit 1
fi

# ----------------------
# Create folder structure
# ----------------------
RFD_OUTPUT="$OUTPUT_DIR/rfdiffusion_output"
PMPNN_OUTPUT="$OUTPUT_DIR/proteinmpnn_output"
AF_INPUT="$OUTPUT_DIR/alphafold_input"
AF_OUTPUT="$OUTPUT_DIR"
#AF_OUTPUT="$OUTPUT_DIR/alphafold_output"
AF_LOG="$OUTPUT_DIR/alphafold_logs"

mkdir -p "$RFD_OUTPUT" "$PMPNN_OUTPUT" "$AF_INPUT" "$AF_OUTPUT" "$AF_LOG"
chmod a+w "$RFD_OUTPUT" "$PMPNN_OUTPUT" "$AF_INPUT" "$AF_OUTPUT" "$AF_LOG"

# ---------------------
# Base Docker command
# ---------------------
base_cmd="docker run --rm -it \
  -v $RFD_OUTPUT:/output \
  rfdiffusion:latest \
  python run_inference.py \
    inference.output_prefix=/output/design_result"

# Mount pdb only if provided
if [[ -n "${PDB_FILE:-}" ]]; then
  base_cmd="docker run --rm -it \
    -v $PDB_FILE:/app/input.pdb \
    -v $RFD_OUTPUT:/output \
    rfdiffusion:latest \
    python run_inference.py \
      inference.output_prefix=/output/design_result"
fi

# -------------------------------
# Add only non-empty parameters
# -------------------------------
declare -A PARAMS=(
  ["inference.input_pdb"]="/app/input.pdb"
  ["contigmap.contigs"]="$CONTIG"
  ["ppi.hotspot_res"]="$HOTSPOT_RES"
  ["inference.num_designs"]="$NUM_DESIGNS"
  ["inference.precision"]="$PRECISION"
  ["diffuser.partial_T"]="$PARTIAL_T"
  ["denoiser.noise_scale_ca"]="$NOISE_SCALE_CA"
  ["denoiser.noise_scale_frame"]="$NOISE_SCALE_FRAME"
)


for key in "${!PARAMS[@]}"; do
  val="${PARAMS[$key]}"
  if [[ -n "$val" ]]; then
    # Skip input_pdb if PDB not provided
    if [[ "$key" == "inference.input_pdb" && -z "${PDB_FILE:-}" ]]; then
      continue
    fi
    base_cmd="$base_cmd $key=$val"
  fi
done

# -----------------------------------
# Run
# -----------------------------------
echo "✅ Running RFdiffusion with:"
echo "$base_cmd"
eval "$base_cmd"

MULTICHAIN=false
for pdb in "$RFD_OUTPUT"/*.pdb; do
  chains=$(grep "^ATOM" "$pdb" | cut -c22 | sort | uniq | wc -l)
  if [[ "$chains" -gt 1 ]]; then
    MULTICHAIN=true
    break
  fi
done

# Choose the correct script
if [[ "$MULTICHAIN" == true ]]; then
  echo "🔗 Detected multichain PDB(s), using script_example_2.py"
  mpnn_script="examples/script_example_2.py"
  fasta_subdir="example_2_outputs"
else
  echo "🧬 Single-chain PDB(s) detected, using script_example_1.py"
  mpnn_script="examples/script_example_1.py"
  fasta_subdir="example_1_outputs"
fi

# Run ProteinMPNN
echo "✅ Running ProteinMPNN..."
docker run --rm -it \
  -v "$RFD_OUTPUT":/inputs \
  -v "$PMPNN_OUTPUT":/outputs \
  pmpnn:latest \
  python "$mpnn_script" \
    --input /inputs \
    --num_seq_per_target 1 \
    --sampling_temp 0.1 \
    --seed 37

input_fasta="$PMPNN_OUTPUT/$fasta_subdir/seqs/design_result_0.fa"
output_fasta="$AF_INPUT/design.fa"

if [[ ! -f "$input_fasta" ]]; then
  echo "❌ FASTA file not found: $input_fasta"
  exit 1
fi

# Extract second sequence
sequence=$(awk '/^>/{n++} n==2{getline; print}' "$input_fasta")

if [[ "$sequence" == *"/"* ]]; then
  echo "🔎 Detected multimer sequence"

  IFS='/' read -ra CHAINS <<< "$sequence"

  echo ">chain_A" > "$output_fasta"
  echo "${CHAINS[0]}" >> "$output_fasta"

  echo ">chain_B" >> "$output_fasta"
  echo "${CHAINS[1]}" >> "$output_fasta"

  af_mode="multimer"

else
  echo "🔎 Detected monomer sequence"

  echo ">designed_seq" > "$output_fasta"
  echo "$sequence" >> "$output_fasta"

  af_mode="monomer"
fi

echo "✅ FASTA for AlphaFold2 ($af_mode):"
cat "$output_fasta"

# ----------------------
# Run AlphaFold Preprocessing
# ----------------------
echo "✅ Running AlphaFold Preprocessing ($af_mode)..."

if [[ "$af_mode" == "multimer" ]]; then
  docker run -it --rm --cap-add SYS_NICE \
    -v "$DATA_DIR":/data \
    -v "$AF_INPUT":/samples \
    -v "$AF_OUTPUT":/output \
    -v "$AF_LOG":/Open-Omics-Acceleration-Framework/applications/alphafold/logs \
    alphafold:pre multimer
else
  docker run -it --rm --cap-add SYS_NICE \
    -v "$DATA_DIR":/data \
    -v "$AF_INPUT":/samples \
    -v "$AF_OUTPUT":/output \
    -v "$AF_LOG":/Open-Omics-Acceleration-Framework/applications/alphafold/logs \
    alphafold:pre
fi

# ----------------------
# Run AlphaFold Inference
# ----------------------
echo "✅ Running AlphaFold Inference ($af_mode)..."

if [[ "$af_mode" == "multimer" ]]; then
  docker run -it --rm --cap-add SYS_NICE \
    -v "$DATA_DIR":/data \
    -v "$AF_INPUT":/samples \
    -v "$AF_OUTPUT":/output \
    -v "$AF_LOG":/Open-Omics-Acceleration-Framework/applications/alphafold/logs \
    alphafold:inf multimer
else
  docker run -it --rm --cap-add SYS_NICE \
    -v "$DATA_DIR":/data \
    -v "$AF_INPUT":/samples \
    -v "$AF_OUTPUT":/output \
    -v "$AF_LOG":/Open-Omics-Acceleration-Framework/applications/alphafold/logs \
    alphafold:inf
fi

# ----------------------
# Done!
# ----------------------
#echo "✅ Full pipeline completed successfully!"
#echo "👉 Outputs:"
#echo "   RFdiffusion:       $RFD_OUTPUT"
#echo "   ProteinMPNN:       $PMPNN_OUTPUT"
#echo "   AlphaFold Input:   $output_fasta"
#echo "   AlphaFold Output:  $AF_OUTPUT"
#echo "   Logs:              $AF_LOG"

echo "🧹 Cleaning up all intermediate output folders..."

docker run --rm \
  --user root \
  -v "$(dirname "$RFD_OUTPUT")":/mnt \
  -v "$(dirname "$PMPNN_OUTPUT")":/mnt2 \
  -v "$(dirname "$AF_INPUT")":/mnt3 \
  -v "$(dirname "$AF_LOG")":/mnt4 \
  rfdiffusion:latest \
  bash -c "rm -rf /mnt/$(basename "$RFD_OUTPUT") /mnt2/$(basename "$PMPNN_OUTPUT") /mnt3/$(basename "$AF_INPUT") /mnt4/$(basename "$AF_LOG")"

echo "✅ Cleanup done. All specified output folders deleted completely. available output in $OUTPUT_DIR "

}
#echo "🧹 Cleaning up intermediate outputs..."
#sudo rm -rf "$RFD_OUTPUT" "$PMPNN_OUTPUT" "$AF_INPUT" "$AF_LOG"
#echo "✅ Cleanup done. Only AlphaFold output kept at: $AF_OUTPUT"


###########################
## super pipeline
###########################

# ---------------------------
# Default values
# ---------------------------

name=default_name
pdb=""
contigs=""
hotspot=""
iterations=50
num_designs=1
symmetry="none"
order=1
chains=""
add_potential=""


# ---------------------------
# Parse arguments
# ---------------------------

without_pdb_db() {
while [[ "$#" -gt 0 ]]; do
  case $1 in
    --output_dir) OUTPUT="$2"; shift ;;
    --mode_name) name="$2"; shift ;;
    --input_pdb) pdb="$2"; shift ;;
    --contigs) contigs="$2"; shift ;;
    --hotspot) hotspot="$2"; shift ;;
    --iterations) iterations="$2"; shift ;;
    --num_designs) num_designs="$2"; shift ;;
    --symmetry) symmetry="$2"; shift ;;
    --order) order="$2"; shift ;;
    --chains) chains="$2"; shift ;;
    --add_potential) add_potential="--add_potential" ;; # No shift
    *) echo "❌ Unknown parameter passed: $1"; exit 1 ;;
  esac
  shift
done

# ---------------------------
# Ensure output directory
# ---------------------------
mkdir -p "$OUTPUT"
chmod a+w "$OUTPUT"

base_cmd="docker run --rm -it \
  -v $OUTPUT:/app/RFdiffusion/outputs \
  -w /app/RFdiffusion\
  rfdiffusion:latest \
  python protein.py"

# If pdb is provided, mount it and set the --pdb flag value
if [[ -n "${pdb:-}" ]]; then
  base_cmd="docker run --rm -it \
    -v $pdb:/app/input.pdb \
    -v $OUTPUT:/app/RFdiffusion/outputs \
    -w /app/RFdiffusion \
    rfdiffusion:latest \
    python protein.py"
fi

# ---------------------------
# Parameters dictionary
# ---------------------------

declare -A PARAMS=(
  ["--name"]="$name"
  ["--pdb"]="/app/input.pdb"
  ["--contigs"]="$contigs"
  ["--hotspot"]="$hotspot"
  ["--iterations"]="$iterations"
  ["--num_designs"]="$num_designs"
  ["--symmetry"]="$symmetry"
  ["--order"]="$order"
  ["--chains"]="$chains"
)

# ---------------------------
# Add only non-empty parameters
# ---------------------------
for key in "${!PARAMS[@]}"; do
  val="${PARAMS[$key]}"
  if [[ -n "$val" ]]; then
    # Skip --pdb if pdb not provided
    if [[ "$key" == "pdb" && -z "${pdb:-}" ]]; then
      continue
    fi
    base_cmd="$base_cmd $key=$val"
  fi
done

echo "✅ Running RFdiffusion with:"
echo "$base_cmd"
eval "$base_cmd"

######################################
# ---------------------------
# Prepare files for ProteinMPNN
# ---------------------------

# These are the expected outputs:
RF_PDB="$OUTPUT/${name}_0.pdb"
RF_CONTIG="$OUTPUT/${name}/final_contigs.txt"

# Check if files exist
if [[ ! -f "$RF_PDB" ]]; then
  echo "❌ Error: PDB file not found: $RF_PDB"
  exit 1
fi

if [[ ! -f "$RF_CONTIG" ]]; then
  echo "❌ Error: Contig file not found: $RF_CONTIG"
  exit 1
fi

# Read the contig string (important: should be single line)
#CONTIG_STRING=$(cat "$RF_CONTIG")
CONTIG_STRING=$(cat "$RF_CONTIG" | tr -d "'")

# ---------------------------
#  Run ProteinMPNN
# ---------------------------
docker run -it \
  -v $OUTPUT:/output \
  -v $RF_PDB:/input \
  -w /ColabDesign \
  pmpnn:latest \
  python colabdesign/rf/designability_test.py \
    --pdb /input \
    --contigs "$CONTIG_STRING" \
    --loc /output

FASTA_FILE=$(find "$OUTPUT" -maxdepth 1 -type f \( -name "*.fa" -o -name "*.fasta" \) | head -n 1)

if [[ -z "$FASTA_FILE" ]]; then
    echo "❌ ERROR: No .fa or .fasta file found inside $OUTPUT"
    exit 1
fi

echo "✅ Found FASTA file: $FASTA_FILE"


#echo "Running AlphaFold preprocess..."

# --------------------
# Run AlphaFold Preprocess
# --------------------
docker run -it --rm \
  -v $OUTPUT:/app/ColabDesign/af_output \
  -v $OUTPUT:/input \
  -v $OUTPUT:/fasta \
  colab_af:pre \
  python colabdesign/rf/designability_test.py \
  --pdb /input/${name}_0.pdb \
  --contigs "$CONTIG_STRING" \
  --fasta /fasta/$(basename "$FASTA_FILE") \
  --loc ./


mkdir -p "$OUTPUT/input"
cp "$OUTPUT/design.fasta" "$OUTPUT/input/"

if [ -f "$OUTPUT/input/design.fasta" ]; then
    echo "Successfully copied design.fasta to fasta directory"
else
    echo "Failed to copy design.fasta" >&2
    exit 1
fi


# Prepare directories
mkdir -p "$OUTPUT/pre_output"
cp -r "$OUTPUT/design" "$OUTPUT/pre_output/"

mkdir -p "$OUTPUT/log_output"

mkdir -p alphafold_data/params

#=============================
# Download Alphafold parameters
#=============================
# Check if AlphaFold parameters are already downloaded (look for a known key file inside params)

if [ -f alphafold_data/params/params_model_1.npz ]; then
    echo "AlphaFold parameters already present, skipping download."
else
    echo "Downloading AlphaFold parameters..."
    curl -fsSL https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar | tar x -C alphafold_data/params
fi

# Export environment variables for AlphaFold
export DATA_DIR=$(realpath ./alphafold_data)
echo "DATA_DIR set to $DATA_DIR"

export SAMPLES_DIR=$OUTPUT/input
export OUTPUT_DIR=$OUTPUT/pre_output
export LOG_DIR=$OUTPUT/log_output


# Run AlphaFold inference

#--------------------------
# Alphafold Inference
#--------------------------

echo "Running AlphaFold inference..."
docker run -it --cap-add SYS_NICE -v $DATA_DIR:/data \
    -v $SAMPLES_DIR:/samples \
    -v $OUTPUT_DIR:/output \
    -v $LOG_DIR:/Open-Omics-Acceleration-Framework/applications/alphafold/logs \
    alphafold:inf

# ===============================
# Postprocess in RFdiffusion
# ===============================

FINAL_OUTPUT_DIR="$OUTPUT/final_af_output"
mkdir -p "$FINAL_OUTPUT_DIR"

echo "🛠️ Running postprocess..."
docker run -it \
    --user $(id -u):$(id -g) \
    -v "$OUTPUT":/output \
    rfdiffusion:latest \
    bash -c "
    mkdir -p /output/final_af_output &&
    for i in {1..5}; do
        export PDB_PATH=\"/output/pre_output/design/unrelaxed_model_\${i}_pred_0.pdb\" &&
        export LENGTHS_PATH=\"/output/pre_output/design/intermediates/lengths.npz\" &&
        export OUTPUT_PDB=\"/output/final_af_output/unrelaxed_model_\${i}.pdb\" &&
        echo \"Processing model \$i...\" &&
        python3 /app/RFdiffusion/scripts/postprocess.py
    done
    "
#export OUTPUT_PDB=\"/output/final_af_output/processed_model_\${i}.pdb\" &&
###################################################
#### Deleting RFdiffusion and ProteinMPNN outputs #
###################################################
echo "🧹 Cleaning up RFdiffusion and ProteinMPNN outputs..."
docker run --rm \
  --user root \
  -v "$OUTPUT":/output \
  rfdiffusion:latest \
  bash -c 'find /output -mindepth 1 -maxdepth 1 ! -name "final_af_output" -exec rm -rf {} +
  if [ -d "/output/final_af_output" ]; then
    mv /output/final_af_output/* /output/ 2>/dev/null
    rm -rf /output/final_af_output
  fi
 '
echo "✅ Cleanup complete — only AlphaFold output retained, available in $OUTPUT"

}

# ---------------
# Main
# ---------------
usage() {
  echo "Usage: $0 --prep [with_pdb_db|without_pdb_db] [other args...]"
  exit 1
}

if [[ "$1" != "--prep" || -z "$2" ]]; then
  usage
fi

PIPELINE="$2"
shift 2

case "$PIPELINE" in
  without_pdb_db) without_pdb_db "$@" ;;
  with_pdb_db) with_pdb_db "$@" ;;
  *) echo "❌ Invalid pipeline type: $PIPELINE"; usage ;;
esac
