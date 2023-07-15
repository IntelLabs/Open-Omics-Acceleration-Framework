SCRIPT_PATH="${BASH_SOURCE:-$0}"
ABS_SCRIPT_PATH="$(realpath "${SCRIPT_PATH}")"
#echo "Value of ABS_SCRIPT_PATH: ${ABS_SCRIPT_PATH}"
ABS_DIRECTORY="$(dirname "${ABS_SCRIPT_PATH}")"
#echo "Value of ABS_DIRECTORY: ${ABS_DIRECTORY}"

echo $1 $2 $3
#path ranks bins
for (( j=0 ; j < $2 ; j++ ))
do
        for ((i=0 ; i < $3 ; i++ ))
        do
                printf -v padded_number "%05d" $(( ($i * $2) + $j ))
                echo $padded_number
                ls ${1}/${padded_number}/*_output.vcf.gz -v >> ${1}/a.txt
        done
done
vcf_list=`cat ${1}/a.txt`

${ABS_DIRECTORY}/../../applications/bcftools/bcftools concat $vcf_list > ${1}/output.vcf.gz

rm ${1}/a.txt
