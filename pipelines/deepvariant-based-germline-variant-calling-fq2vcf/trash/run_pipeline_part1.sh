#!/bin/bash

set -e
SCRIPT_PATH="${BASH_SOURCE:-$0}"
ABS_SCRIPT_PATH="$(realpath "${SCRIPT_PATH}")"
ABS_DIRECTORY="$(dirname "${ABS_SCRIPTINPUT_DIR_PATH}")"

INDIR=$INPUT_DIR
OUTDIR=$OUTPUT_DIR
REFDIR=$REF_DIR
#* ranks: Number of mpi process that we want the pipeline to run on
#* threads/shards: parameters to different tools in the pipeline, calculated as below
ppn=$2


Sockets=$(cat compute_config | grep -E '^Socket\(s\)' | awk  '{print $2}')   #2
Cores=$(cat compute_config | grep -E '^Core\(s\)' | awk  '{print $4}')  #56
Thread=$(cat compute_config | grep -E '^Thread' | awk  '{print $4}')  #2

a=$(( $(( ${Cores}*${Thread}*${Sockets} / $ppn )) - 2*${Thread} ))   #24 (Four threads are removed for IO)
b=$(( $(( ${Cores}*${Sockets} )) / $ppn ))   #14

if [ $a -lt 1 ]
then
    echo 'Number of cpus are less to run the pipeline.'
    exit 0
fi

N=$1
PPN=$2
CPUS=$a
THREADS=$a
SHARDS=$b
REF=$(basename "$3")  #Change to your reference file
READ1="$4" #Change your read files
READ2="$5"
BINDING=socket
Container=docker

if [ $# -gt 5 ]
then
        Container="$6"
fi

echo "Output directory: $OUTDIR"
mkdir -p ${OUTDIR}
#It is assumed that if reference file is .gz then it is converted using create_reference_index.sh or pcluster_reference_index.sh script.
file_ext=${REF##*.}

if [ "${file_ext}" = "gz" ]
then
        REF=$(basename "$REF" .gz )
        if ! [ -f $REFDIR/${REF} ]; then
                echo "File $REFDIR/${REF} does not exist."
                exit 0
        fi
fi

echo Starting run with $N ranks, $CPUS threads,$THREADS threads, $SHARDS shards, $PPN ppn.
# -in -sindex are required only once for indexing.
# Todo : Make index creation parameterized.
mpiexec -bootstrap ssh -bind-to $BINDING -map-by $BINDING --hostfile hostfile -n $N -ppn $PPN python -u test_pipeline_part1.py --input $INDIR --output  $OUTDIR $TEMPDIR --refdir $REFDIR --index $REF --read $READ1 $READ2 --cpus $CPUS --threads $THREADS --shards $SHARDS --container_tool "$Container"  2>&1 | tee ${OUTDIR}/log_part1.txt

#echo "Pipeline finished. Output vcf can be found at: $OUTPUT_DIR/output.vcf.gz"
