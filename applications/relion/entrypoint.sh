#!/bin/bash
source ~/intel/oneapi/2025.0/oneapi-vars.sh
APPEXE=/opt/relion_5.0_cpu_benchmark_prefix/bin/relion_refine_mpi
DEFAULT_DATA_DIR=/opt/relion_5.0/relion_benchmark
cd "$DEFAULT_DATA_DIR" || exit 1

NUMMPI=$((16+1))
NUMTHR=8
POOLSIZE=4
NUMITER=25

DEFAULT_MODE="${1:-3d}"  
MODE="${DEFAULT_MODE#-}" 
shift                

case "$MODE" in
  2d)
    echo "➡️ Running 2D Classification"
    mkdir -p 2D
    exec mpirun -n $NUMMPI $APPEXE --i Particles/shiny_2sets.star --dont_combine_weights_via_disc --ctf --tau2_fudge 2 --particle_diameter 360 --K 200 --zero_mask --oversampling 1 --psi_step 6 --offset_range 5 --offset_step 2 --norm --scale --random_seed 0 --pad 2 --o 2D/2D --pool $POOLSIZE --j $NUMTHR --iter $NUMITER --cpu
    ;;
  3d)
    echo "➡️ Running 3D Classification with SYCL"
    mkdir -p 3D
    exec mpirun -n $NUMMPI $APPEXE --i Particles/shiny_2sets.star --ref emd_2660.map:mrc --firstiter_cc --ini_high 60 --dont_combine_weights_via_disc --ctf --tau2_fudge 4 --particle_diameter 360 --K 6 --flatten_solvent --zero_mask --oversampling 1 --healpix_order 2 --offset_range 5 --offset_step 2 --sym C1 --norm --scale --pad 2 --random_seed 0 --o 3D/3D --pool $POOLSIZE --j $NUMTHR --iter $NUMITER --cpu
    ;;
  autorefine)
    echo "➡️ Running 3D AutoRefine"
    mkdir -p 3D_AUTO
    exec mpirun -n $NUMMPI $APPEXE --i Particles/shiny_2sets.star --ref emd_2660.map:mrc --firstiter_cc --ini_high 60 --dont_combine_weights_via_disc --ctf --particle_diameter 360 --flatten_solvent --zero_mask --oversampling 1 --healpix_order 2 --offset_range 5 --offset_step 2 --sym C1 --norm --scale --auto_refine --split_random_halves --auto_local_healpix_order 4 --low_resol_join_halves 40 --random_seed 1 --pad 2 --o 3D_AUTO/3D_AUTO --pool $POOLSIZE --j $NUMTHR --cpu
    ;;
  *)
   echo "➡️ Running custom RELION command: $DEFAULT_MODE $@"
    exec $DEFAULT_MODE "$@"
    ;;
esac

