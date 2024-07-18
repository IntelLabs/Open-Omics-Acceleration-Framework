source config

mode="multifq"
if [ "$#" == "2" ] || [ "$#" == "3" ]
    then
        mode=$1
        se_mode=$2
        sso=$3
else
    echo "<exec> <pragzip/flatmode/fqprocessonly/multifq> [se/pe] [sso (for single socket only execution)]"
    exit
fi
echo "mode: "$mode
echo "se/pe: $se_mode"
[[ "$3" == "sso" ]] && echo "single socket only"
source miniconda3/bin/activate distbwa

#echo "localhost" > hostfile
hostname > hostfile
num_nodes=1
#num_nodes=`cat hostfile | wc -l`
#first_ip=`head -n 1 hostfile`
#echo $first_ip
#ssh ${first_ip} 'lscpu' > compute_config
lscpu > compute_config

semode=""
read_type="--read_type short"
[[ "$se_mode" == "se" ]] && semode="--se_mode"
[[ "$READ_TYPE" == "long" ]] && read_type="--read_type long"
[[ "$READ_TYPE" == "long" ]] && semode="--se_mode"  ## resetting to se reads for long reads


num_cpus_per_node=$(cat compute_config | grep -E '^CPU\(s\)' | awk  '{print $2}')
num_numa=$(cat compute_config | grep '^NUMA node(s)' | awk '{print $3}')
if [ "$sso" == "sso" ]
then
    tot_socket=$(cat compute_config | grep -E '^Socket'| awk  '{print $2}')
    num_socket=1
    num_cpus_per_node=$(( $num_cpus_per_node / $tot_socket ))
    num_numa=$(( $num_numa / $tot_socket ))
else   
    num_socket=$(cat compute_config | grep -E '^Socket'| awk  '{print $2}')
fi

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

#th=`expr ${num_physical_cores_per_numa} / 2`  #${num_physical_cores_per_numa}  ##20
th=16
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
#if [ "$R3" == "" ]
#then
#    R3="NONE"
#fi

INDIR=$INPUT_DIR
OUTDIR=$OUTPUT_DIR
mkdir -p $OUTDIR/logs
#* ranks: Number of mpi process that we want the pipeline to run on
#* threads/shards: parameters to different tools in the pipeline, calculated as below
ppn=${ranks_per_node}


Sockets=$(cat compute_config | grep -E '^Socket\(s\)' | awk  '{print $2}')   #2
Cores=$(cat compute_config | grep -E '^Core\(s\)' | awk  '{print $4}')  #56
Thread=$(cat compute_config | grep -E '^Thread' | awk  '{print $4}')  #2

if [ "$sso" == "sso" ]
then
    Sockets=1
fi
echo "Socket: $Sockets"
echo "Cores: $Cores"
echo "Thread: $Thread"

a=$(( $(( ${Cores}*${Thread}*${Sockets} / $ppn )) - 2*${Thread} ))   #24 (Four threads are removed for IO)
b=$(( $(( ${Cores}*${Sockets} )) / $ppn ))   #14
c=$(( $b*${Thread} ))
if [ $a -lt 1 ]
then
    echo 'Number of cpus are less to run the pipeline.'
    exit 0
fi
echo "a: $a"
echo "b: $b"
echo "c: $c"
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
READI1=${I1[@]}

echo "reads:"
echo "READ1 $READ1"
echo "READ2 $READ2"
echo "READ3 $READ3"
echo "READI1 $READI1"
echo "PREFIX $PREFIX"
echo "SUFFIX $SUFFIX"

whitelist=""
read_structure=""
barcode_orientation=""
bam_size=""
sample_id=""
output_format=""
params=""
outfile=""
istart=""

[[ -n $WHITELLIST ]] && whitelist="--whitelist $WHITELIST" && echo "Whitelist: $whitelist"
[[ -n $READ_STRUCTURE ]] && read_structure="--read_structure $READ_STRUCTURE" && echo "Read Structure: $read_structure"
[[ -n $BARCODE_ORIENTATION ]] && barcode_orientation="--barcode_orientation $BARCODE_ORIENTATION" && echo "Barcode Orientation: $barcode_orientation"
[[ -n $BAM_SIZE ]] && bam_size="--bam_size $BAM_SIZE" && echo "BAM Size: $bam_size"
[[ -n $SAMPLE_ID ]] && sample_id="--sample_id $SAMPLE_ID" && echo "sample_id: $sample_id"
[[ -n $OUTPUT_FORMAT ]] && output_format="--output_format $OUTPUT_FORMAT" && echo "output_format: $output_format"
#[[ -n $PARAMS ]] && params="--params $PARAMS" && echo "params: $params"
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
echo "Params" : "${PARAMS}"
mkdir -p ${OUTDIR}

echo Starting run with $N ranks, $CPUS threads,$THREADS threads, $PPN ppn.
# -in -sindex are required only once for indexing.
# Todo : Make index creation parameterized.

# Check if number of ranks equals number of splits
[[ "$R1" != "" ]] && echo "Read1: "$R1;
[[ "$R2" != "" ]] && echo "Read2: "$R2;
[[ "$R3" != "" ]] && echo "Read3: "$R3; 
#echo $R2;
#echo $R3; 

if [ "$1" == "multifq" ]
then
    R1_LEN=`echo $R1 | tr ' ' '\n' | wc -l`
    R3_LEN=`echo $R3 | tr ' ' '\n' | wc -l`

    if [ "$N" != "$R1_LEN" ]; then
        echo "Error: Number of ranks ("$N") does not equal number of splits ("$R1_LEN"). Program failed."
        exit 1
    fi

    # Check if number of R1 and R3 fastq files is equal
    if [ "$R1_LEN" != "$R3_LEN" ]; then
        echo "Error: Number of R1 fastq files doesnt equal number of R3 files. Program failed."
        exit 1
    fi
fi

exec=dist_bwa.py
#mpiexec -bootstrap ssh -bind-to $BINDING -map-by $BINDING --hostfile hostfile -n $N -ppn $PPN python -u $exec --input $INDIR --output  $OUTDIR $TEMPDIR $REFDIR --index $REF --preads $READ1 $READ2 --cpus $CPUS --threads $THREADS --keep_unmapped ${whitelist} ${read_structure} ${barcode_orientation} ${bam_size} ${outfile} ${istart} ${sample_id} ${output_format} --prefix $PREFIX --suffix $SUFFIX --params "${PARAMS}" --mode $mode --read_type $READ_TYPE  2>&1 | tee ${OUTDIR}log.txt


envs=""
envs=" -env I_MPI_PIN_DOMAIN=${c}:compact -env I_MPI_PIN_ORDER=range "
#echo "envs: "$envs
if [ "$mode" == "fqprocess" ]
then
    mpiexec -bootstrap ssh -bind-to $BINDING -map-by $BINDING --hostfile hostfile -n $N -ppn $PPN python -u $exec --input $INDIR --output  $OUTDIR $TEMPDIR $REFDIR --index $REF --read1 $READ1 --read2 $READ2 --read3 $READ3 --cpus $CPUS --threads $THREADS --keep_unmapped --params "${PARAMS}" --mode $mode ${semode} ${outfile} ${read_type}  2>&1 | tee ${OUTDIR}/logs/log.txt
    
elif [ "$mode" == "pragzip" ] || [ "$mode" == "flatmode" ]
then
    #/data/swtools/intel/mpi/2021.12/bin/
    mpiexec  --hostfile hostfile -n $N -ppn $PPN $envs python -u $exec --input $INDIR --output  $OUTDIR $TEMPDIR $REFDIR --index $REF --preads $READ1 $READ2 --cpus $CPUS --threads $THREADS --keep_unmapped --params "${PARAMS}" --mode  $mode ${semode} ${outfile} ${read_type} 2>&1 | tee ${OUTDIR}/logs/log.txt
    #mpiexec  -bootstrap ssh -bind-to $BINDING -map-by $BINDING --hostfile hostfile -n $N -ppn $PPNpython -u $exec --input $INDIR --output  $OUTDIR $TEMPDIR $REFDIR --index $REF --preads $READ1 $READ2 --cpus $CPUS --threads $THREADS --keep_unmapped --mode $mode ${semode} ${read_type} ${read_type}  2>&1 | tee ${OUTDIR}log.txt

else
    mpiexec -bootstrap ssh -bind-to $BINDING -map-by $BINDING --hostfile hostfile -n $N -ppn $PPN python -u $exec --input $INDIR --output  $OUTDIR $TEMPDIR $REFDIR --index $REF --preads $READ1 $READ2 --cpus $CPUS --threads $THREADS --keep_unmapped ${whitelist} ${read_structure} ${barcode_orientation} ${bam_size} ${outfile} ${istart} ${sample_id} ${output_format} --prefix $PREFIX --suffix $SUFFIX --params "${PARAMS}" --mode $mode ${semode} ${read_type} 2>&1 | tee ${OUTDIR}/logs/log.txt
    
fi
