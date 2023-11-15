## fq2sortedBAM: OpenOmics' Distributed BWA-MEM2 pipeline for getting sorted BAM output using mpi4py (distributed message passing library) Python package
### Overview:
The pipeline takes input fastq files and produces sorted BAM file using the following four stages:
1. fastqprocessing (from Broad's warptools)
2. bwa-mem2
3. sorting (using samtools, partially overlapped (the communication part) with bwa-mem2 compute)
4. concatenation of the bam files produced in the previous stages (using samtools)

Although the code supports distributed memory parallel execution on a set of compute nodes, right now we've setup the scripts to execute on a single compute node (can have multiple sockets or NUMA domains). We plan to enable the code for distributed parallel execution on multiple compute nodes in the future.
**We've tested the code on GCP C3 instance as well as on on-prem cluster node containing Intel(R) Xeon(R) 4th generation scalable processor**.


### Installation (2-step process):
1. ```./basic_setup_ubuntu.sh```
Notes: Step 1 installs the following pre-requisites. The script is only written for ubuntu. In case of other OSes, the users are instructed to manually install these dependencies. The installaion requires sudo access.
       - make
       - gcc/g++
       - autoconf
       - wget
       - git
       - numactl
       - zlib1g-dev
       - libbz2-dev
       - liblzma-dev
       - libncurses5-dev

2. ```./install.sh```
Notes:  Installs miniconda, creates an environment in it, and installs required packages. It also installs bwa-mem2, samtools, and warptools (fastqprocess).

Important Notes:
1. The code is tested on GCP C3 instance (176 and 88 vCPUs) w/ ubuntu 22.04 w/ gcc 11.3.0 using the install script.
2. Thde code is also tested on an onprem cluster node.

### How to run (2-step process):
1. ```./print_config.sh```   ## this reports the number of required splits of the inpute file. Use this to tune bam_size parameter in step 2 below.

2. Setup "config" file in the current folder. It supplies parameters to the pipeline. The "config" file contains the following parameters:
```
export INPUT_DIR=~/input/
export OUTPUT_DIR=~/output/
REF=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
R1=r1.fastq.gz
R2=r2.fastq.gz
R3=r3.fastq.gz
WHITELIST="whitelist.txt"
READ_STRUCTURE=""
BARCODE_ORIENTATION=""
BAM_SIZE="9"
OUTFILE="final"
PARAMS=""
ISTART="False"
```

INPUT_DIR - directory path containing reference genome, bwa-mem2 index, and input read files
OUTPUT_DIR - directory path for the output SAM/BAM files and intermediate files
REF - reference sequence file
R1, R2, R3 - input read files as input to fqprocess
WHITELIST - 10x genomics whitelist file for fqprocess
READ_STRUCTURE - read structure for fqprocess
BARCODE_ORIENTATION - barcode orientation for fqprocess
BAM_SIZE - threshold on the size of the output files from fqprocess
OUTFILE - name of the output sorted bam file, by default it is set to "final"
PARAMS - parameters to bwa-mem2 mapping, excluding "-t <threads>" as value for -t is set automatically
ISTART    - flag for bwa-mem2 index creation. You can use this if your bwa-mem2 index is not already created or present in the INPUT_DIR location provided in config file.

3. Run
```./run_bwa.sh <mode>```
e.g.: ```./run_bwa.sh fqprocess```

Commentary:
1. run_bwa.sh automatically determines the optimal configuration (number of MPI ranks) for the distributed bwa-mem2 run based on the number of cores/sockets/NUMA on the compute node. Each MPI rank performs bwa-mem2 and downstream processing in parallel.
2. The fastqprocess executes only on rank 0 and writes the processed fastq files in the **current working directory** as fastq\_R1\_\<i\> and fastq\_R2\_\<i\>, here "i" represents the numbering of the split files e.g. fastq\_R1\_0, fastq\_R1\_1, ....
3. Each MPI rank executes on one processed fastq file (paired-end).
4. Note that, if the number of fastq files generated is greater than number of MPI ranks, then the extra fastq files won't be processed by bwa-mem2. Else, if the number of fastq files generated is less than number of MPI ranks then the code breaks.
5. **Please go through the notes below on how to tune bam_size parameter so that number of fastq files produced by fastqprocess matches the number of MPI ranks.**

Important Notes:
1. To tune bam_size parameter, we are providing **'print_config.sh'** script: ```./print_config.sh```.
This script prints the number of MPI ranks that will be created by run_bwa.sh during execution.
**This should help the user in setting "bam_size" parameter value as input to run_bwa.sh so that the fastqprocess can create number of fastq files equal to number of MPI ranks.**

2. Since the fastq files produced by fastqprocess are stored in current working directory, these files may interfere with subsequent runs of the pipeline. The user can delete the previous fastq files from previous runs before commencing the next run.

### Modes:
Distributed bwa-mem2 supports 3 different modes:
1. fqprocess: fqprocess mode does fastqprocessing before invoking bwa-mem2 on the processed fastq files.   - **This mode is specifically designed for Broad's requirements.**
   - It takes fastq files R1, R2, R3 as input (provided in config file)
   - It works with gzipped as well as with un-compressed files
   - fastqprocess uses bam_size parameters to split the input files into multiple files. The number of files created should be equal to the number of  MPI ranks to get better performance. Moreover, all split files with equal size are desirable to get proper load balancing across MPI rank and for better performance.
   - By default the output file name is 'final.sorted.bam', it can be provided as command line parameter in config file with parameter name 'OUTFILE'
   - dist_bwa.py contains the mpi4py code. it supports options for providing arguments required for fqprocess and bwa-mem2. These parameters can be provided in the config file

2. pragzip: invokes bwa-mem2 on input gzipped fastq file and produces sorted bam files as output
   - Splitting of the input files into multiple files on disk is not required.
   - To facilitate using single file for multiple MPI ranks, pragzip mode uses pragzip library to first index the input file. Then, this file is equally divided (or chunked) among the MPI ranks. These equal size chunks are used as input to bwa-mem2 processes -- one chunk per bwa-mem2. The pragzip index of the input files helps in locating and reading these chunks from the input gzip file.
   - This mode needs gzip fastq files
   - Each bwa-mem2 process figures out its chunk from the input gzipped fastq files
   - Multiple SAM files from each bwa-mem2 process are sorted using SAMTools
   - The sorted BAM files are concatenated to produce the final bam file

3. flatmode: Same as pragzip mode but does not sort and concat the output SAM/BAM files from each rank
