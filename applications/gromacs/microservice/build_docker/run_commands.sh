#!/bin/bash
source /opt/conda/bin/activate grms_env
source /grmcs/gmx_mkl_prefix/bin/GMXRC

file="$1"
num_cpus=$(nproc)
current_ntmpi=$num_cpus
timestamp=$(date +"%Y%m%d_%H%M%S")

exit_on_error() { echo "Error: $1"; exit 1; }

echo "Preprocessing protein structure..."
grep -v HOH "${file}" > prot_clean.pdb || exit_on_error "Failed to clean PDB file"

echo "Running pdb2gmx..."
# *** explicit outputs so we always have these filenames
gmx pdb2gmx -f prot_clean.pdb -o prot_pros.gro -p topol.top -i posre.itp -water spce -ff amber99sb -ignh \
  || exit_on_error "pdb2gmx failed"

echo "Defining box dimensions..."
gmx editconf -f prot_pros.gro -o prot_box.gro -c -d 1.0 -bt cubic || exit_on_error "editconf failed"

echo "Adding solvent..."
gmx solvate -cp prot_box.gro -cs spc216.gro -o prot_solv.gro -p topol.top || exit_on_error "solvate failed"

echo "Preparing for ion addition..."
gmx grompp -f mdtut_ions.mdp -c prot_solv.gro -p topol.top -o ions.tpr -maxwarn 2 \
  || exit_on_error "grompp for ions failed"

echo "Adding ions..."
# *** don't hard-code a numeric index; use the group name for reliability
printf "SOL\n" | gmx genion -s ions.tpr -o prot_solv_ions.gro -p topol.top -pname NA -nname CL -neutral \
  || exit_on_error "genion failed"

echo "Preparing energy minimization..."
gmx grompp -f mdtut_minim.mdp -c prot_solv_ions.gro -p topol.top -o em.tpr -maxwarn 1 \
  || exit_on_error "grompp for minimization failed"

while [ "$current_ntmpi" -ge 1 ]; do
  echo "Running EM with ntmpi=$current_ntmpi..."
  time /usr/bin/time -v gmx mdrun -ntmpi "$current_ntmpi" -v -deffnm em
  mdrun_exit_code=$?
  if [ $mdrun_exit_code -eq 0 ]; then
    echo "EM stage completed successfully with ntmpi=$current_ntmpi"
    break
  fi
  # *** fix precedence & handle missing file quietly
  if [ -f em.log ] && ( grep -q "There is no domain decomposition" em.log || grep -q "^Fatal error:" em.log ); then
    echo "Domain decomposition failed for ntmpi=$current_ntmpi. Reducing MPI ranks..."
    grep -m1 "^Fatal error:" em.log || true
    rm -f em.log
    current_ntmpi=$(( current_ntmpi / 2 ))
  else
    exit_on_error "Unexpected failure in EM stage. Check em.log for details."
  fi
done

echo "Preparing NVT equilibration..."
gmx grompp -f mdtut_nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr \
  || exit_on_error "grompp for NVT failed"

echo "Running NVT mdrun..."
gmx mdrun -ntmpi "$current_ntmpi" -v -deffnm nvt || exit_on_error "NVT mdrun failed"

echo "Preparing NPT equilibration..."
gmx grompp -f mdtut_npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr \
  || exit_on_error "grompp for NPT failed"

echo "Running NPT mdrun..."
gmx mdrun -ntmpi "$current_ntmpi" -v -deffnm npt || exit_on_error "NPT mdrun failed"

echo "Preparing MD simulation..."
gmx grompp -f mdtut_md.mdp -c npt.gro -t npt.cpt -p topol.top -o md01.tpr \
  || exit_on_error "grompp for MD failed"

echo "Running MD mdrun..."
gmx mdrun -ntmpi "$current_ntmpi" -v -deffnm md01 || exit_on_error "MD mdrun failed"

echo "Simulation completed successfully with ntmpi=$current_ntmpi"

