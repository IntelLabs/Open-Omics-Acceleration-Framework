#!/bin/sh

SCRIPT_PATH="${BASH_SOURCE:-$0}"
ABS_SCRIPT_PATH="$(realpath "${SCRIPT_PATH}")"
ABS_DIRECTORY="$(dirname "${ABS_SCRIPTINPUT_DIR_PATH}")"

INDIR=$INPUT_DIR
OUTDIR=$OUTPUT_DIR

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
REF=$(basename "$3")  #Change to your reference file
READ1=$(basename "$4")  #Change your read files
READ2=$(basename "$5")
READ3=$(basename "$6")
BINDING=socket
mode=$7
echo "mode: "$mode

echo "Output directory: $OUTDIR"
mkdir -p ${OUTDIR}

echo Starting run with $N ranks, $CPUS threads,$THREADS threads, $PPN ppn.
# -in -sindex are required only once for indexing.
# Todo : Make index creation parameterized.
#exec=dist_bwa.py
exec=dist_bwav2.py
mpiexec -bootstrap ssh -bind-to $BINDING -map-by $BINDING --hostfile hostfile -n $N -ppn $PPN python -u $exec --input $INDIR --output  $OUTDIR $TEMPDIR $REFDIR --index $REF --read1 $READ1 --read2 $READ2 --read3 $READ3 --cpus $CPUS --threads $THREADS --keep_unmapped --mode $mode   2>&1 | tee ${OUTDIR}log.txt



#if [ "$mode" == "fqprocess" ]
#then
#    mpiexec -bootstrap ssh -bind-to $BINDING -map-by $BINDING --hostfile hostfile -n $N -ppn $PPN python -u $exec --input $INDIR --output  $OUTDIR $TEMPDIR $REFDIR --index $REF --read1 $READ1 --read2 $READ2 --read3 $READ3 --cpus $CPUS --threads $THREADS --keep_unmapped --mode $mode   2>&1 | tee ${OUTDIR}log.txt
#else
#        mpiexec -bootstrap ssh -bind-to $BINDING -map-by $BINDING --hostfile hostfile -n $N -ppn $PPN python -u $exec --input $INDIR --output  $OUTDIR $TEMPDIR $REFDIR --index $REF --preads $READ1 $READ2 --cpus $CPUS --threads $THREADS --keep_unmapped --mode $mode  2>&1 | tee ${OUTDIR}log.txt
#
#fi
