#!/bin/sh
#docker pull google/deepvariant:1.5.0
lscpu > compute
num_cpus_per_node=$(cat compute | grep -E '^CPU\(s\)' | awk  '{print $2}')

INPUT=~/HG001/
OUTPUT=~/HGOO1/OUTPUT/
echo $OUTPUT
mkdir -p $OUTPUT
#python test_pipe_bwa.py --input $INPUT --output $OUTPUT --index GCA_000001405.15_GRCh38_no_alt_analysis_set.fna --read HG001.novaseq.pcr-free.30x.R1.fastq.gz HG001.novaseq.pcr-free.30x.R2.fastq.gz --cpus 108 --threads 108 --shards 112 
python test_pipe_bwa.py --input $INPUT --output $OUTPUT --index Homo_sapiens_assembly38.fasta --read HG001.novaseq.pcr-free.30x.R1.fastq.gz HG001.novaseq.pcr-free.30x.R2.fastq.gz --cpus $num_cpus_per_node --threads $num_cpus_per_node --shards $num_cpus_per_node # 2>&1 | tee ${OUTPUT}log.txti

echo "Ouput files are inside "$OUTPUT" folder"
