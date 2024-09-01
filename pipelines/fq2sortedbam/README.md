## fq2SortedBAM: OpenOmics' Secondary Analysis Pipeline
### Overview:
The pipeline takes input fastq files and produces sorted BAM file through the following stages:
1. Sequence Alignment: bwa-mem2 for short reads, mm2-fast (Minimap2) for long reads (PacBio, ONT)
2. SAMSort (Using SAMTools)

### Modes:
fq2SortedBAM supports 4 different modes:  
1. sortedbam: It takes fastq reads files as input and outputs sorted BAM file  
2. flatmode: It takes fastq reads files as input and outputs multiple (equal to the number of ranks) unsorted SAM files  
3. fqprocessonly: Invokes Broad's fastqprocess tool on raw input files  - **This mode is specifically designed for Broad's requirements. More details are in README.md.old**  
4. multifq: Invokes bwa-mem2 on the processed (e.g by fastqprocess) fastq files - **This mode is specifically designed for Broad's requirements. More details are in README.md.old**  


## Use Docker
### Docker build:  
```
wget https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/blob/main/pipelines/fq2sortedbam/Dockerfile
wget https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/blob/main/pipelines/fq2sortedbam/config.yaml  
cp config.yaml <inputdir>
docker build -t fq2bam .
docker save fq2bam:latest > fq2bam.tar  
```
### 
Docker run:
```
docker load -i fq2bam.tar
Setup <input>/config.yaml keys with appropriate values  
docker run -v <inputdir>:/input <outdir>:/out <refdir>:/refdir <tempdir>:/tempdir fq2bam:latest /app/Open-Omics-Acceleration-Framework/pipelines/fq2sortedbam/run_bwa.sh sortedbam /input/config.yaml
```

## Use Source Code  
### Installation:
```
git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
cd Open-Omics-Acceleration-Framework/blob/main/pipelines/fq2sortedbam/
bash install.sh <onprem/cloud>  
```

### Run:
```
Setup ./config.yaml keys with appropriate values  
bash run_bwa.sh sortedbam ./config.yaml
```
### Notes:  
1. Individual pipeline tools are present Open-Omics-Acceleration-Framework/blob/main/applications    
2. To understand various parameters to these tools, you can access their ```man``` page at this location
3. You can setup the parameters to these tools using ```params``` parameters in config.yaml  

### Runtime setup (config.yaml):  
bwa:
  dindex: True/False WIP HERE
  params: +R "@RG\tID:RG1\tSM:RGSN1"
  rindex: 'True'
dataset:
  index: GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
  input: ./data/
  outfile: short.se.sam
  output: ./out/
  read1: HG001.novaseq.pcr-free.30x.R1_10M.fastq.gz
  read2: ''
  read_type:  short
  refdir: ./data/
  tempdir: ''
fqprocess:
  bam_size: '5'
  barcode_orientation: FIRST_BP_RC
  output_format: FASTQ
  prefix: multiome-practice-may15_arcgtf
  read3: ''
  read_structure: 16C
  readi1: ''
  sample_id: ''
  suffix: trimmed_adapters.fastq.gz
  whitelist: whitelist.txt
mm2:
  params: ' -ax map-hifi '



## [OLD README] fq2sortedBAM: OpenOmics' Distributed BWA-MEM2 pipeline for getting sorted BAM output using mpi4py (distributed message passing library) Python package
### Overview:
The pipeline takes input fastq files and produces sorted BAM file using the following four stages:
1. fastqprocessing (from Broad's warptools)
2. bwa-mem2
3. sorting (using samtools, partially overlapped (the communication part) with bwa-mem2 compute)
4. concatenation of the bam files produced in the previous stages (using samtools)

The code supports distributed memory parallel execution on a set of compute nodes using MPI. However, right now we've setup the scripts to execute on a single compute node (can have multiple sockets or NUMA domains) using multiple MPI ranks. We plan to enable the code for distributed parallel execution on multiple compute nodes in the future.  
**We've tested the code on GCP C3 instance as well as on on-prem cluster node containing Intel(R) Xeon(R) 4th generation scalable processor**.  

MPI:  
MPI is a library for parallel execution where each MPI rank executes one process where processes interact with each other through message passing interface (MPI).  

### Modes:
Distributed bwa-mem2 supports 3 different modes:  
1. fqprocessonly: Invokes Broad's fastqprocess tool on raw input files.  - **This mode is specifically designed for Broad's requirements.**  
   - It takes R1, R2, R3, (I1) fastq files as input
   - It produces processed fastq files. The number of files produced are dependent on the bam_size parameter
   - Input paramters can be provided through config file
   - fastqprocess uses bam_size parameters to split the input files into multiple files.

2. multifq: Invokes bwa-mem2 on the processed (e.g by fastqprocess) fastq files - **This mode is specifically designed for Broad's requirements.**
   - It takes multiple R1 and R3 fastq files as input (provided in config file). e.g. fastq_R1_0.fastq.gz, fastq_R1_1.fastq.gz, fastq_R3_0.fastq.gz, fastq_R3_1.fastq.gz, ...
   - The number of files of R1 and R3 should be equal to the number MPI ranks
   - The number of input files created should be equal to the number of  MPI ranks to get better performance. Moreover, all split files with equal size are desirable to get proper load balancing across MPI rank and for better performance.
   - It works with gzipped as well as with un-compressed files
   - By default the output file name is 'final.sorted.bam', it can be provided as command line parameter in config file with parameter name 'OUTFILE'
   - dist_bwa.py contains the mpi4py code. it supports options for providing arguments required for and bwa-mem2. These parameters can be provided in the config file

3. pragzip: invokes bwa-mem2 on input gzipped fastq file and produces sorted bam files as output
   - This mode uses processesed fastq files as input and hence does not execute fastqprocess
   - This mode needs gzipped fastq files
   - It uses pragzip library to split the input gzipped fastq files into chunks and process those chunks in parallel using multiple instances of bwa-mem2
   - It uses pragzip library to first index the input file, and then this file is equally divided (or chunked) among bwa-mem2 instances. These equal size chunks are used as input to bwa-mem2 -- one chunk per bwa-mem2 instance. The pragzip index of the input files helps in locating and reading these chunks from the input gzip file.
   - Multiple SAM files from each bwa-mem2 instance are sorted using SAMTools
   - The sorted BAM files are concatenated to produce the final bam file  

4. flatmode: Same as pragzip mode but does not sort and concat the output SAM/BAM files from each MPI rank  



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

### How to run (3-step process):
1. ```./print_config.sh```   ## this reports the number of required splits of the inpute file. Use this to tune bam_size parameter in step 2 below.  

2. Setup "config" file in the current folder. It supplies parameters to the pipeline. The "config" file contains the following parameters:  
```
export INPUT_DIR=~/input/
export OUTPUT_DIR=~/output/
export LD_PRELOAD=<absolute_path>/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/libmimalloc.so.2.0:$LD_PRELOAD
REF=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
R1=r1.fastq.gz
R2=r2.fastq.gz
R3=r3.fastq.gz
I1=""
R1PREFIX="fastq_R1"
R3PREFIX="fastq_R3"
WHITELIST="whitelist.txt"
READ_STRUCTURE="16C"
BARCODE_ORIENTATION="FIRST_BP_RC"
BAM_SIZE="9"
SAMPLE_ID=""
OUTPUT_FORMAT="FASTQ"
OUTFILE="final"
PARAMS=""
ISTART="False"
```

INPUT_DIR - directory path containing reference genome, bwa-mem2 index, and input read files  
OUTPUT_DIR - directory path for the output SAM/BAM files and intermediate files  
REF - reference sequence file  
OUTFILE - name of the output sorted bam file, by default it is set to "final"  
PARAMS - parameters to bwa-mem2 mapping, excluding "-t <threads>" as value for -t is set automatically  
ISTART    - flag for bwa-mem2 index creation. You can use this if your bwa-mem2 index is not already created or present in the INPUT_DIR location provided in config file.  

%% multifq mode parameters: set these parameters if you are running multifq mode  
R1PREFIX - fastq R1 file name prefix. e.g fastq_R1_0.fastq.gz, fastq_R1_1.fastq.gz, here "_0", "_1" can be used by MPI ranks to identify their file to process. Use this for fastq input when running multifq mode. **These files should be present in INPUT_DIR**.  
R3PREFIX - fastq R3 file name prefix. e.g fastq_R3_0.fastq.gz, fastq_R3_1.fastq.gz, here "_0", "_1" can be used by MPI ranks to identify their file to process. Use this for fastq input when running multifq mode. **These files should be present in INPUT_DIR**.  

%%fqpocess mode paramters: set these parameters if you are running fastqprocess mode  
R1, R2, R3, I1 - input read files as input to fqprocess. These files shouuld be present at INPUT_DIR location.  
WHITELIST - 10x genomics whitelist file for fqprocess  
READ_STRUCTURE - read structure for fqprocess  
BARCODE_ORIENTATION - barcode orientation for fqprocess  
BAM_SIZE - threshold on the size of the output files from fqprocess  
SAMPLE_ID="" - parameter to fqprocess  
OUTPUT_FORMAT="FASTQ" - parameters to fqprocess for ouptut file format  


3. Run
```./run_bwa.sh <mode>```  
e.g.: ```./run_bwa.sh multifq```  

Commentary:
1. **User needs to execute all the above three steps for "fqprocess" mode to get the desired number of processed fastq files. If the input files are available then for multifq mode just execute step 2 & 3 above.**
2. run_bwa.sh automatically determines the optimal configuration (number of MPI ranks) for the distributed bwa-mem2 run based on the number of cores/sockets/NUMA on the compute node. Each MPI rank performs bwa-mem2 and downstream processing in parallel.
3. The fastqprocess executes only on rank 0 and writes the processed fastq files in the **current working directory** as fastq\_R1\_\<i\> and fastq\_R2\_\<i\>, here "i" represents the numbering of the split files e.g. fastq\_R1\_0, fastq\_R1\_1, ....
4. Each MPI rank executes on one processed fastq file (paired-end).
5. Note that, if the number of fastq files generated is greater than number of MPI ranks, then the extra fastq files won't be processed by bwa- mem2. Else, if the number of fastq files generated is less than number of MPI ranks then the code breaks.
6. **Please go through the notes below on how to tune bam_size parameter so that number of fastq files produced by fastqprocess matches the number  of MPI ranks.**

Important Notes:
1. To tune bam_size parameter, we are providing **'print_config.sh'** script: ```./print_config.sh```.
This script prints the number of MPI ranks that will be created by run_bwa.sh during execution.
**This should help the user in setting "bam_size" parameter value as input to run_bwa.sh so that the fastqprocess can create number of fastq files equal to number of MPI ranks.**
