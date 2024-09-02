#*************************************************************************************
#                           The MIT License
#
#   Intel OpenOmics - fq2sortedbam pipeline
#   Copyright (C) 2023  Intel Corporation.
#
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the
#   "Software"), to deal in the Software without restriction, including
#   without limitation the rights to use, copy, modify, merge, publish,
#   distribute, sublicense, and/or sell copies of the Software, and to
#   permit persons to whom the Software is furnished to do so, subject to
#   the following conditions:
#
#   The above copyright notice and this permission notice shall be
#   included in all copies or substantial portions of the Software.
#
#   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
#   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
#   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
#   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
#   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
#   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
#   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
#   SOFTWARE.
#
#Authors:  Vasimuddin Md <vasimuddin.md@intel.com>; Babu Pillai <padmanabhan.s.pillai@intel.com>;
#*****************************************************************************************/

#!/usr/bin/bash
set -e

trap 'echo "Error occurred"' ERR

source ./miniforge3/bin/activate fq2bam
echo "Your conda env @: "$CONDA_PREFIX

config=""
sso="None"
if [ "$#" == "2" ] || [ "$#" == "3" ]
    then
        mode=$1
        config=$2
        [[ "$#" == "3" ]] && sso=$3
else
    echo "<exec> <sortedbam/flatmode/fqprocessonly/multifq> <config_file> [sso (for single socket only execution)]"
    exit
fi

#Note: "##### Note: Currently, this code only supports single node. "
#Note: "##### I've deliberately disabled distributed runs for now. "
#Note: "##### Contact: <vasimuddin.md@intel.com>"
#echo ""

num_nodes=1
echo "run mode: "$mode
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
echo "[Info] Running $N ranks, each with $THREADS threads ..."
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
mpiexec -bootstrap ssh -n $N -ppn $PPN -bind-to $BINDING -map-by $BINDING  --hostfile hostfile  python -u $exec --cpus $CPUS --threads $THREADS ${runmode} ${CONFIG} --keep_unmapped  2>&1 | tee logs/log.txt
echo "[Info] The output log file is at logs/log.txt"
