#!/usr/bin/env bash
set -euo pipefail

input_file=""
work_root="/tmp/work_gromacs"

for arg in "$@"; do
  case "$arg" in
    input_file=*) input_file="${arg#input_file=}" ;;
    work_root=*)  work_root="${arg#work_root=}" ;;
    *) echo "Unknown arg: $arg" >&2; exit 2 ;;
  esac
done

[[ -n "${input_file}" ]] || { echo "input_file is required" >&2; exit 3; }
[[ -f "${input_file}" ]] || { echo "PDB not found: ${input_file}" >&2; exit 4; }

pdb_base="$(basename "${input_file}" .pdb)"
job_dir="${work_root}/${pdb_base}"
mkdir -p "${job_dir}"

# Copy protein + workflow files into isolated job dir
cp -f "${input_file}" "${job_dir}/input.pdb"
cp -f /input/run_commands.sh "${job_dir}/run_commands.sh"
cp -f /input/mdtut_*.mdp "${job_dir}/"

cd "${job_dir}"
chmod +x ./run_commands.sh

./run_commands.sh "./input.pdb"

