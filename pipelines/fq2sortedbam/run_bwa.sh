#!/usr/bin/bash
set -e

trap 'echo "Error occurred"' ERR

source ./miniforge3/bin/activate fq2bam
echo "conda_prefix: "$CONDA_PREFIX

config=""
sso="None"
if [ "$#" == "2" ] || [ "$#" == "3" ]
    then
        mode=$1
        #se_mode=$2
        config=$2
        [[ "$#" == "3" ]] && sso=$3
else
    echo "<exec> <pragzip/flatmode/fqprocessonly/multifq> <config_file> [sso (for single socket only execution)]"
    exit
fi
echo "##### Note: Currently, this code only supports single node. "
echo "##### I've deliberately disabled distributed runs for now. "
echo "##### Contact: <vasimuddin.md@intel.com>"
echo ""

num_nodes=1
echo "mode: "$mode
echo "config: "$config
#echo "se/pe: $se_mode"
[[ "$sso" == "sso" ]] && echo "single socket only"

#echo "localhost" > hostfile
hostname > hostfile
#semode=""
CONFIG=""
#read_type="--read_type short"
#[[ "$se_mode" == "se" ]] && semode="--se_mode"
[[ "$config" != "" ]] && CONFIG="-y $config"
runmode="--mode $mode"

#lscpu > lscpu.txt
chmod +x hwconfig.py
#ls -lh  hwconfig.py
python hwconfig.py $sso $num_nodes > hwconfig
#python hwconfig.py "sso"   > hwconfig

source hwconfig
#rm lscpu.txt
rm hwconfig

BINDING=socket
mkdir -p logs

exec=dist_bwa.py
#echo $I_MPI_PIN_DOMAIN
#-genv I_MPI_PIN_DOMAIN=$I_MPI_PIN_DOMAIN

#echo $N
#echo $PPN
#echo $exec
#echo $CONFIG
mpiexec -bootstrap ssh -n $N -ppn $PPN -bind-to $BINDING -map-by $BINDING  --hostfile hostfile  python -u $exec --cpus $CPUS --threads $THREADS --keep_unmapped ${runmode} ${CONFIG} --keep_unmapped 2>&1 | tee logs/log.txt
echo "The output log file is at logs/log.txt"

##if [ "$mode" == "fqprocess" ]
##then
##    mpiexec -bootstrap ssh -bind-to $BINDING -map-by $BINDING --hostfile hostfile -n $N -ppn $PPN python -u $exec --cpus $CPUS --threads $THREADS --keep_unmapped ${runmode} ${semode}   2>&1 | tee ${OUTDIR}/logs/log.txt
##
##    #--input $INDIR --output  $OUTDIR $TEMPDIR $REFDIR --index $REF --read1 $READ1 --read2 $READ2 --read3 $READ3 --cpus $CPUS --threads $THREADS --keep_unmapped --params "${PARAMS}" --mode $mode ${semode} ${outfile} ${read_type}  2>&1 | tee ${OUTDIR}/logs/log.txt
##    
##elif [ "$mode" == "pragzip" ] || [ "$mode" == "flatmode" ]
##then
##    #/data/swtools/intel/mpi/2021.12/bin/
##    #mpiexec  --hostfile hostfile -n $N -ppn $PPN $envs python -u $exec --input $INDIR --output  $OUTDIR $TEMPDIR $REFDIR --index $REF --preads $READ1 $READ2 --cpus $CPUS --threads $THREADS --keep_unmapped --params "${PARAMS}" --mode  $mode ${semode} ${outfile} ${read_type} 2>&1 | tee ${OUTDIR}/logs/log.txt
##    N=4
##    mpiexec  -bootstrap ssh -bind-to $BINDING -map-by $BINDING --hostfile hostfile -n $N -ppn $PPN python -u $exec --input $INDIR --output  $OUTDIR $TEMPDIR $REFDIR --index $REF --preads $READ1 $READ2 --cpus $CPUS --threads $THREADS --keep_unmapped --mode $mode ${semode} ${read_type}  2>&1 | tee ${OUTDIR}/logs/log.txt
##
##else
##    mpiexec -bootstrap ssh -bind-to $BINDING -map-by $BINDING --hostfile hostfile -n $N -ppn $PPN python -u $exec --input $INDIR --output  $OUTDIR $TEMPDIR $REFDIR --index $REF --preads $READ1 $READ2 --cpus $CPUS --threads $THREADS --keep_unmapped ${whitelist} ${read_structure} ${barcode_orientation} ${bam_size} ${outfile} ${istart} ${sample_id} ${output_format} --prefix $PREFIX --suffix $SUFFIX --params "${PARAMS}" --mode $mode ${semode} ${read_type} 2>&1 | tee ${OUTDIR}/logs/log.txt
##    
##fi
