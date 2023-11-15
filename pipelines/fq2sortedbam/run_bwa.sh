source config

mode="fqprocess"
if [ "$#" == "1" ]
    then
    mode=$1
else
    echo "<exec> <pragzip/flatmode/fqprocess>"
    exit
fi
echo "mode: "$mode
source miniconda3/bin/activate distbwa

echo "localhost" > hostfile
num_nodes=1
#num_nodes=`cat hostfile | wc -l`
#first_ip=`head -n 1 hostfile`
#echo $first_ip
#ssh ${first_ip} 'lscpu' > compute_config
lscpu > compute_config


num_cpus_per_node=$(cat compute_config | grep -E '^CPU\(s\)' | awk  '{print $2}')
num_socket=$(cat compute_config | grep -E '^Socket'| awk  '{print $2}')
num_numa=$(cat compute_config | grep '^NUMA node(s)' | awk '{print $3}')
num_cpus_all_node=`expr ${num_cpus_per_node} \* ${num_nodes}`
threads_per_core=$(cat compute_config | grep -E '^Thread' | awk  '{print $4}')
echo "#############################################"
echo "Total number of CPUs: $num_cpus_all_node"
echo "Number of sockets: "$num_socket
echo "Number of NUMA domains: "$num_numa

num_physical_cores_all_nodes=`expr ${num_cpus_all_node} / ${threads_per_core}`
num_physical_cores_per_nodes=`expr ${num_cpus_per_node} / ${threads_per_core}`
num_physical_cores_per_socket=`expr ${num_physical_cores_all_nodes} / ${num_socket}`
num_physical_cores_per_numa=`expr ${num_physical_cores_all_nodes} / ${num_numa}`
echo "Num physical cores per nodes: "$num_physical_cores_per_nodes
echo "Num physical cores per socket: "$num_physical_cores_per_socket
echo "Num physical cores per numa: "$num_physical_cores_per_numa

th=`expr ${num_physical_cores_per_numa} / 2`  #${num_physical_cores_per_numa}  ##20
if [ $th -le 10 ]
then
    th=${num_physical_cores_per_numa}
fi

while [ $num_physical_cores_per_nodes -gt $th ]
do
    num_physical_cores_per_nodes=`expr $num_physical_cores_per_nodes / 2`
done

num_physical_cores_per_rank=$num_physical_cores_per_nodes
total_num_ranks=`expr ${num_physical_cores_all_nodes} / ${num_physical_cores_per_rank}`

ranks_per_node=`expr ${total_num_ranks} / ${num_nodes}`
echo "Cores per node: "$num_physical_cores_per_nodes
echo "Total number of ranks: "${total_num_ranks}
echo "Ranks per node: "${ranks_per_node}
echo "#############################################"
if [ "$R3" == "" ]
then
    R3="NONE"
fi

INDIR=$INPUT_DIR
OUTDIR=$OUTPUT_DIR

#* ranks: Number of mpi process that we want the pipeline to run on
#* threads/shards: parameters to different tools in the pipeline, calculated as below
ppn=${ranks_per_node}


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

#N=$1
#PPN=$2
N=${total_num_ranks}
PPN=${ranks_per_node}
CPUS=$a
THREADS=$a
#REF=$(basename "$3")  #Change to your reference file
#READ1=$(basename "$4")  #Change your read files
#READ2=$(basename "$5")
#READ3=$(basename "$6")
BINDING=socket
#mode=$7
echo "mode: "$mode

## parameters
READ1=${R1[@]}
READ2=${R2[@]}
READ3=${R3[@]}
echo "reads:"
echo "READ1 $READ1"
echo "READ2 $READ2"
echo "READ3 $READ3"

whitelist=""
read_structure=""
barcode_orientation=""
bam_size=""
params=""
outfile=""
istart=""

[[ -n $WHITELLIST ]] && whitelist="--whitelist $WHITELIST" && echo "Whitelist: $whitelist"
[[ -n $READ_STRUCTURE ]] && read_structure="--read_structure $READ_STRUCTURE" && echo "Read Structure: $read_structure"
[[ -n $BARCODE_ORIENTATION ]] && barcode_orientation="--barcode_orientation $BARCODE_ORIENTATION" && echo "Barcode Orientation: $barcode_orientation"
[[ -n $BAM_SIZE ]] && bam_size="--bam_size $BAM_SIZE" && echo "BAM Size: $bam_size"
[[ -n $PARAMS ]] && params="--params $PARAMS" && echo "params: $params"
[[ -n $OUTFILE ]] && outfile="--outfile $OUTFILE" && echo "outfile: $outfile"
if [ "$ISTART" == "True" ]
then
    istart="--istart"
    echo "istart: $istart"
fi

#echo "Whitelist: $whitelist"
#echo "Read Structure: $read_structure"
#echo "Barcode Orientation: $barcode_orientation"
#echo "BAM Size: $bam_size"
#echo "params: $params"
#echo "outfile: $outfile"


echo "Input directory: $INDIR"
echo "Output directory: $OUTDIR"
mkdir -p ${OUTDIR}

echo Starting run with $N ranks, $CPUS threads,$THREADS threads, $PPN ppn.
# -in -sindex are required only once for indexing.
# Todo : Make index creation parameterized.

exec=dist_bwa.py
mpiexec -bootstrap ssh -bind-to $BINDING -map-by $BINDING --hostfile hostfile -n $N -ppn $PPN python -u $exec --input $INDIR --output  $OUTDIR $TEMPDIR $REFDIR --index $REF --read1 $READ1 --read2 $READ2 --read3 $READ3 --cpus $CPUS --threads $THREADS --keep_unmapped ${whitelist} ${read_structure} ${barcode_orientation} ${bam_size} ${params} ${outfile} ${istart}  --mode $mode   2>&1 | tee ${OUTDIR}log.txt



#if [ "$mode" == "fqprocess" ]
#then
#    mpiexec -bootstrap ssh -bind-to $BINDING -map-by $BINDING --hostfile hostfile -n $N -ppn $PPN python -u $exec --input $INDIR --output  $OUTDIR $TEMPDIR $REFDIR --index $REF --read1 $READ1 --read2 $READ2 --read3 $READ3 --cpus $CPUS --threads $THREADS --keep_unmapped --mode $mode   2>&1 | tee ${OUTDIR}log.txt
#else
#        mpiexec -bootstrap ssh -bind-to $BINDING -map-by $BINDING --hostfile hostfile -n $N -ppn $PPN python -u $exec --input $INDIR --output  $OUTDIR $TEMPDIR $REFDIR --index $REF --preads $READ1 $READ2 --cpus $CPUS --threads $THREADS --keep_unmapped --mode $mode  2>&1 | tee ${OUTDIR}log.txt
#
#fi
