#!/bin/bash 
set -e
SCRIPT_PATH="${BASH_SOURCE:-$0}"
ABS_SCRIPT_PATH="$(realpath "${SCRIPT_PATH}")"
#echo "Value of ABS_SCRIPT_PATH: ${ABS_SCRIPT_PATH}"
ABS_DIRECTORY="$(dirname "${ABS_SCRIPT_PATH}")"

cd $ABS_DIRECTORY
source config
filename=${REF}
file_ext=${filename##*.}
file_name_without_extension=$(basename "$filename" .gz )


if [ ${file_ext} == 'gz' ]
then
	echo "Refecence file is decompressing..."
        gzip -d ${REF_DIR}/${filename}
	REF=${file_name_without_extension}
fi
ref=${REF_DIR}/${REF}
mkdir -p ${OUTPUT_DIR}

echo "Checking the index files for $ref"
ls ${ref}*

# mem2 index
echo "Creating FM-index for the reference sequence ${ref}"
cd ../../../../applications/bwa-mem2
./bwa-mem2 index $ref &> ${OUTPUT_DIR}/bwa_mem2_index_log
cd - &> /dev/null


# samtool idfai index
echo "Creating fai index for the reference sequence ${ref}"
cd ../../../../applications/samtools
./samtools faidx $ref &> ${OUTPUT_DIR}/samtools_fai_log
cd - &> /dev/null


echo "The list of all index files created."
ls ${ref}*
if [ -z $1 ]
then 
		echo "Index files are created."
	else
		echo "Index files are created release instance by typing: 'scancel $1' "
fi
