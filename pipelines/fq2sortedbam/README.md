## fq2SortedBAM: OpenOmics' Genomics Secondary Analysis Pipeline
### Overview:
The pipeline takes input fastq files and produces sorted BAM file through the following stages:
1. Sequence Alignment: bwa-mem2 for short reads, mm2-fast (Accelerated Minimap2) for long reads (PacBio, ONT)
2. SAMSort (Using SAMTools)

### Modes:
fq2SortedBAM supports 4 different modes:  
1. ```sortedbam```: It takes fastq reads files and reference genome as input and outputs sorted BAM file  
2. ```flatmode```: It takes fastq reads files and reference genome as input, and outputs multiple (equal to the number of ranks created) unsorted SAM files  
3. ```fqprocessonly```: Custom mode, not for general use
4. ```multifq```: Custom mode, not for general use  


## Use Docker
### Docker build:  
```
wget https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/releases/download/3.0/Source_code_with_submodules.tar.gz  
tar -xzf Source_code_with_submodules.tar.gz  
cp Open-Omics-Acceleration-Framework/pipelines/fq2sortedbam/Dockerfile .
cp Open-Omics-Acceleration-Framework/pipelines/fq2sortedbam/config.yaml <inputdir>
docker build -t fq2bam .
docker save fq2bam:latest > fq2bam.tar     ## this step is optional  
```

### Setup Input Parameters:
Setup \<inputdir\>/config.yaml (described below) with appropriate values

### Docker run:
```
docker load -i fq2bam.tar      ## optional, if the image is build on the same machine or is already loaded  
docker run -v <inputdir>:/input -v <outdir>:/out -v <refdir>:/refdir -v <tempdir>:/tempdir fq2bam:latest bash run_bwa.sh sortedbam /input/config.yaml
```
Note:  
\<inputdir\>: Location of the local directory containing read files read1 & read2  
\<refdir\>: Location of the local directory containing reference sequence file ref  
\<outdir\>: Location of the local directory for output files SAM/BAM  
\<tempdir\>: Location of the local directory for temporary files (defaults to \<outdir\>)  


## Use Source Code  
### Installation:
```
wget https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/releases/download/3.0/Source_code_with_submodules.tar.gz  
tar -xzf Source_code_with_submodules.tar.gz  
cd Open-Omics-Acceleration-Framework/pipelines/fq2sortedbam/
bash install.sh <onprem/cloud>  ## onprem mode: Manually install the depenendies present in basic_setup_ubuntu.sh as it needs sudo access
```

### Setup Input Parameters:
Setup ./config.yaml (described below)  with appropriate values

### Run:
```
bash run_bwa.sh sortedbam ./config.yaml
```

## General Notes:  
1. Individual pipeline tools are present in applications folder    
2. To understand various parameters to these tools, you can access their ```man``` page  
3. You can setup the parameters of these tools using ```params``` variable in ```config.yaml```    

## Setup config.yaml<sup>1<sup>:  
1. bwa: bwa-mem2 related parameters     
   - dindex: dtype=bool, values="True/False", if True it creates .fai index files for input reads files    
   - params: dtype=string, the command line paramteres to bwa-mem2 mapping run e.g. '+R "@RG\tID:RG1\tSM:RGSN1"'  
   - rindex: dtype=bool, 'values=True/False', if "True" it creates bwa-mem2 index for the reference genome    
  
2. dataset:  data details  
   - __index__: dtype:string, Input reference genome file name e.g. "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz"  
   - **input**: dtype=string, folder location of input reads files e.g. "**/input**"  
   - outfile: dtype=string, output SAM/BAM file name(s) e.g. short.se.sam. Default value: "final_fq2bam"    
   - **output**: dtype=string, folder location of output SAM/BAM files e.g. "**/out**"  
   - **read1**: dtype=string, input reads file1 name e.g. "HG001.novaseq.pcr-free.30x.R1.fastq.gz"  
   - read2: dtype=string, input reads file2 name e.g. "HG001.novaseq.pcr-free.30x.R2.fastq.gz"  
   - **read_type**: dtype:string, values=**short/long**, short read mapping using bwa-mem2, long read mapping using mm2-fast  
   - **refdir**: dtype=string, folder location of reference genome and its index files e.g. "**/refdir**"
   - tempdir: dtype=string, folder location for storing intermedaite files e.g. "/out/". In case of none, output folder is used  
 
3. fqprocess: custom mode variables  
    - bam_size: dtype=int,, value=5  
    - barcode_orientation: dtype=string, values=FIRST_BP_RC  
    - output_format: dtype=string, value=FASTQ  
    - prefix: dtype=string,, value=multiome-practice-may15_arcgtf  
    - read3: dtype=string, value=''  
    - read_structure: dtype=string, value=16C  
    - readi1: dtype=string, value=''  
    - sample_id: dtype=int, value=''  
    - suffix: dtype=string, value=trimmed_adapters.fastq.gz  
    - whitelist: dtype=string, value=whitelist.txt 
  
4. mm2: mm2-fast related parameters  
   - params: dtype=string, the command line paramteres to mm2-fast mapping run e.g.' -ax map-hifi '    

<sup>**1**</sup> **Parameters in bold are mandatory ones**
