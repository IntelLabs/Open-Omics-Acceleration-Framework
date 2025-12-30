#!/usr/bin/env bash
set -euo pipefail

complex_dir=""
ligand=""

for arg in "$@"; do
  case "$arg" in
    complex_dir=*) complex_dir="${arg#complex_dir=}" ;;
    input_file=*)  ligand="${arg#input_file=}" ;;
    *) echo "Unknown arg: $arg" >&2; exit 2 ;;
  esac
done

cd "${complex_dir}"

lig_base="$(basename "${ligand}" .pdbqt)"
cfg="vina_${lig_base}_config.txt"

if [[ ! -f "${cfg}" ]]; then
  echo "Config not found: ${cfg}" >&2
  exit 3
fi

# Docking output directory (per complex)
out_dir="results"
mkdir -p "${out_dir}"

out_file="${out_dir}/${lig_base}_out.pdbqt"

threads="${OMP_NUM_THREADS:-1}"
vina --config "${cfg}" --cpu "${threads}" --out "${out_file}"
