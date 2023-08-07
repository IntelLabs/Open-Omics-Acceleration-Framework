source config

ref=${INPUT_DIR}/${REF}

echo "Checking the index files for $ref"
ls ${ref}*

# mem2 index
echo "Creating FM-index for the reference sequence ${ref}"
cd ../../../../applications/bwa-mem2
./bwa-mem2 index $ref &> bwa_mem2_index_log
cd - &> /dev/null


# samtool idfai index
echo "Creating fai index for the reference sequence ${ref}"
cd ../../../../applications/samtools
./samtools faidx $ref &> samtools_fai_log
cd - &> /dev/null


echo "The list of all index files created."
ls ${ref}*



