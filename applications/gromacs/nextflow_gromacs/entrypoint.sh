#!/bin/bash
set -eo pipefail

source /opt/conda/bin/activate grms_env
source /grmcs/gmx_mkl_prefix/bin/GMXRC

if [[ "${1:-}" == "/bin/bash" || "${1:-}" == "bash" ]]; then
    exec /bin/bash
fi

if [[ "$#" -eq 0 ]]; then
    echo "Error: No PDB file provided. Please provide a PDB file as an argument."
    echo "  Run full workflow: docker run -v \$(INPUT_GMX_CPU):/input -v \$(OUTPUT_GMX_CPU):/output -it grms_dock <protein.pdb>"
    echo "  Run specific GROMACS command: docker run -v \$(INPUT_GMX_CPU):/input -v \$(OUTPUT_GMX_CPU):/output -it grms_dock gmx <command>"
    exit 1

elif [[ "$1" == "gmx" ]]; then
    shift
    exec gmx "$@"

else
    pdb_file="$1"

    if [[ ! -f "/input/$pdb_file" ]]; then
        echo "Error: File '/input/$pdb_file' not found."
        exit 1
    fi

    echo "Running full workflow with PDB file: /input/$pdb_file"

    if [[ -n "${CUSTOM_SCRIPT:-}" && -f "/input/${CUSTOM_SCRIPT}" ]]; then
        echo "Using custom script: /input/${CUSTOM_SCRIPT}"
        exec "/input/${CUSTOM_SCRIPT}" "/input/$pdb_file"
    else
        echo "Using default script: /input/run_commands.sh"
        exec /input/run_commands.sh "/input/$pdb_file"
    fi
fi

