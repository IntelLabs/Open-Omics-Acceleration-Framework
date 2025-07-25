#!/bin/bash
set -e
SCRIPT_PATH="${BASH_SOURCE:-$0}"
ABS_SCRIPT_PATH="$(realpath "${SCRIPT_PATH}")"
#echo "Value of ABS_SCRIPT_PATH: ${ABS_SCRIPT_PATH}"
ABS_DIRECTORY="$(dirname "${ABS_SCRIPT_PATH}")"
#echo "Value of ABS_DIRECTORY: ${ABS_DIRECTORY}"

echo $1 $2 $3 $4
#path ranks bins
total=$(( ($2 * $3) ))
for (( j=0 ; j < $total ; j++ ))
do
        printf -v padded_number "%05d" $j
        echo $padded_number
        ls ${1}/${padded_number}/output.vcf.gz -v >> ${1}/a.txt
        
done
vcf_list=`cat ${1}/a.txt`

${ABS_DIRECTORY}/../../applications/bcftools/bcftools concat $vcf_list > ${1}/${4}.vcf.gz

rm ${1}/a.txt
