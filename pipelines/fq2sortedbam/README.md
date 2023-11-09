## OpenOmics'  Distributed BWA-MEM2 using mpi4py Python pacakge.
Takes input fastq files and produces sorted BAM file in four stages:
1. fastqprocessing
2. bwa-mem2
3. sorting (partially overalapped with bwa-mem2 compute)
4. concatenation of the bam files produces in previous stages


Installation:
Supported OS: ubuntu

```./install.sh```

Installation does the following:
1. Install miniconda, creates an environment in it, and required packages
2. Installs pre-requisites:

3. Installs bwa-mem2, samtools, and warptools (fastqprocess)


### How to run:
setup "config" file in current folder as:
```
export INPUT_DIR=/cold_storage/omics/genomes-indexes-reads/human/reads/broad/
export OUTPUT_DIR=/data/nfs_home/mvasimud/myworkspace/openomics/aws_all/output/
REF=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
R1=r1.fastq.gz
R2=r2.fastq.gz
R3=r3.fastq.gz
```

INPUT_DIR: directory path containing reference genome, bwa-mem2 index, and input read files
OUTPUT_DIR: directory path for the ouptut SAM/BAM files and intermediate files
REF: reference sequence file
R1, R2, R3: input read files

./autorun.sh <mode>
```./autorun.sh fqprocess```

autorun.sh  automatically figures out the optimal conguration (number of ranks) for the distributed bwa-mem2 run based on the number of cores/sockets/NUMA on the compute platform. Each rank executes bwa-mem2 and downstream processing in parallel on each rank.


dist-bwa-mem2 supports 3 different modes:
1. fqprocess: fqprocess does fastqprocessing before invoking bwa-mem2 on the processed fastq files.    - This mode specifically designed for Broad's requirement.
   - It takes fastq files R1, R2, R3 as input
   - fastqprocess using bam-size parameters split the input files into multiple files. The number of files created should be equal to the number of ranks to get better performance. Moreover all split files with equal size is desirable to get minimum load imbalance across rank and hence better performance.
   - By default the output file name is 'final.sorted.bam', it can be provided as command line parameter to the commandline in run_bwa.sh
   - dist_bwav2.py contains the mpi4py code. it supports options for providing arguments required for fqprocess and bwa-mem2. These parameters can be provided to mpiexec command line in run_bwa.sh

2. pragzip: invokes bwa-mem2 on input gzipped fastq file and produces sorted bam files as output
   - Each bwa-mem2 process figures out its chunk from the input gzipped fastq file(s)
   - multiple sam files from each bwa-mem2 process are sorted using SAMTools
   - The sorted BAM files are merged to produce the final bam file

3. flatmode: Same as pragzip mode but does not sort and merge the output SAM/BAM files from each rank
