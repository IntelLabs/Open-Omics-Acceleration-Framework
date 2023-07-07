#!/bin/sh
SCRIPT_PATH="${BASH_SOURCE:-$0}"
ABS_SCRIPT_PATH="$(realpath "${SCRIPT_PATH}")"
ABS_DIRECTORY="$(dirname "${ABS_SCRIPTINPUT_DIR_PATH}")"

INDIR=$INPUT_DIR
OUTDIR=$OUTPUT_DIR
#INDIR=/lfs/lfs12/ashish/testdata/
#OUTDIR=${ABS_DIRECTORY}/output/
INDEX=GRCh38_chr1.fna

N=16
CPUS=24
THREADS=24
SHARDS=14
FILEBASE=HG001_
#BINDING=numa
BINDING=socket
PPN=8

#sh run_pipeline.sh 256 24 24 14 8 
[[ $# -gt 0 ]] && N="$1"
[[ $# -gt 1 ]] && CPUS="$2"
[[ $# -gt 2 ]] && THREADS="$3"
[[ $# -gt 3 ]] && SHARDS="$4"
[[ $# -gt 4 ]] && PPN="$5"
[[ $SHARDS -eq 0 ]] && SHARDS=$THREADS


echo $OUTDIR
mkdir -p ${OUTDIR}

echo Starting run with $N ranks, $CPUS threads,$THREADS threads, $SHARDS shards, $PPN ppn.

mpiexec -bootstrap ssh -bind-to $BINDING -map-by $BINDING --hostfile hostfile -n $N -ppn $PPN python -u test_pipeline_final.py --input $INDIR --output  $OUTDIR $TEMPDIR $REFDIR --index $INDEX --read ${FILEBASE}R1_chr1.fastq.gz ${FILEBASE}R2_chr1.fastq.gz --cpus $CPUS --threads $THREADS --shards $SHARDS |& tee ${OUTDIR}log.txt
