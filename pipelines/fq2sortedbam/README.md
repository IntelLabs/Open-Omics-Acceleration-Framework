## OpenOmics'  Distributed BWA-MEM2 pipeline for sorted bam using mpi4py Python package.  
The pipeline takes input fastq files and produces sorted BAM file using the following four stages:  
1. fastqprocessing (from Broad's warptools)  
2. bwa-mem2  
3. sorting (partially overlapped (the communication part) with bwa-mem2 compute)  
4. concatenation of the bam files produced in the previous stages  
  
  
Installation:  
Supported OS: ubuntu 22.04  
  
```./install.sh```    ## install.sh calls basic_setup_ubuntu.sh to install pre-requisites for ubuntu  
  
Installation does the following:  
1. Installs miniconda, creates an environment in it, and installs required packages  
2. Installs pre-requisites:  
   - make, gcc/g++, wget, git, numactl, zlib1g-dev,  libbz2-dev, liblzma-dev  
3. Installs bwa-mem2, samtools, and warptools (fastqprocess)  
  
  
### How to run:  
setup "config" file -- supplies paramters to the pipeline -- in the current folder as:  
```  
export INPUT_DIR=~/input/  
export OUTPUT_DIR=~/output/  
REF=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz  
R1=r1.fastq.gz  
R2=r2.fastq.gz  
R3=r3.fastq.gz  
WHITELIST=""  
READ_STRUCTURE=""  
BARCODE_ORIENTATION=""  
BAM_SIZE="9"  
OUTFILE="final"  
PARAMS=""  
ISTART="False"  
```  
  
INPUT_DIR: directory path containing reference genome, bwa-mem2 index, and input read files  
OUTPUT_DIR: directory path for the output SAM/BAM files and intermediate files  
REF: reference sequence file  
R1, R2, R3: input read files as input to fqprocess  
WHITELIST - 10x genomics whitelist file for fqprocess.  
READ_STRUCTURE - read structure for fqprocess.  
BARCODE_ORIENTATION - barcode orientation for fqprocess.  
BAM_SIZE - threshold on the size of the output files from fqprocess.  
OUTFILE - name of the output sorted bam file, by default it is set to "final".  
PARAMS - parameters to bwa-mem2 mapping, excluding "-t <threads>" as value for -t is set automatically.  
ISTART    - flag for bwa-mem2 index creation. You can use this if your bwa-mem2 index is not already created or present in the INPUT_DIR location provided in config file.  
  
  
./run_bwa.sh <mode>  
```./autorun.sh fqprocess```  
  
run_bwa.sh automatically determines the optimal configuration (number of ranks) for the distributed bwa-mem2 run based on the number of cores/sockets/NUMA on the compute platform. Each rank executes bwa-mem2 and downstream processing in parallel. The fqprocess executes only on rank 0 and writes the processed fastq files in the current working directory as fastq_R1_<i> and fastq_R2_<i>, here "i" represents the numbering of the split files e.g. fastq_R1_0, fastq_R1_1, .... Note that if the number of fastq files generated is greater than number of MPI ranks, then the extra fastq files won't be processed by bwa-mem2. Else, if the number of fastq files generated is less then number of MPI rank then the code breaks.  
  
  
Please note:  
1. To tune bam_size parameter, we are providing 'print_condig.sh' script.  
This script prints the number of ranks that will be created by run_bwa.sh during execution.  
This should helps the user in setting bam_size parameter value as input to run_bwa.sh so that the fqprocess can create number of fastq files equal to number of MPI ranks.  
  
2. Since the fastq files produced by fastqprocess tools are stored in current working directory, these files may interfere with subsequent runs of the pipeline. The user can delete the previous fastq files from previous run.  
  
  
dist-bwa-mem2 supports 3 different modes:  
1. fqprocess: fqprocess does fastqprocessing before invoking bwa-mem2 on the processed fastq files.   - This mode is specifically designed for Broad's requirements.  
   - It takes fastq files R1, R2, R3 as input (provided in config file)  
   - fastqprocess uses bam_size parameters to split the input files into multiple files. The number of files created should be equal to the number of MPI ranks to get better performance. Moreover, all split files with equal size are desirable to get proper load balancing across MPI rank and for better performance.  
   - By default the output file name is 'final.sorted.bam', it can be provided as command line parameter in config file with parameter nam 'OUTFILE'  
   - dist_bwa.py contains the mpi4py code. it supports options for providing arguments required for fqprocess and bwa-mem2. These parameters can be provided in config file  
  
2. pragzip: invokes bwa-mem2 on input gzipped fastq file and produces sorted bam files as output  
   - Each bwa-mem2 process figures out its chunk from the input gzipped fastq file(s)  
   - multiple sam files from each bwa-mem2 process are sorted using SAMTools  
   - The sorted BAM files are merged to produce the final bam file  
  
3. flatmode: Same as pragzip mode but does not sort and merge the output SAM/BAM files from each rank  
  
