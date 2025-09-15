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
```
```bash
docker save fq2bam:latest > fq2bam.tar     ## this step is optional  
docker build --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy -t fq2bam .  
```


### Docker run:
```
docker load -i fq2bam.tar      ## optional, if the image is build on the same machine or is already loaded    
docker run  -v <readsdir>:/readsdir -v <outdir>:/outdir -v <refdir>:/refdir fq2bam:latest python run_fq2sortedbam.py --ref /refdir/<reference> --reads  /readsdir/<read1>  /readsdir/<read2>  --output /outdir/<outBAMfile> 
```
Note:  
\<readsdir\>: Location of the local directory containing read files read1 & read2  
\<refdir\>: Location of the local directory containing reference sequence file ref  
\<outdir\>: Location of the local directory for output files SAM/BAM  


## Use Source Code  
### Installation:
```
wget https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/releases/download/3.0/Source_code_with_submodules.tar.gz  
tar -xzf Source_code_with_submodules.tar.gz  
cd Open-Omics-Acceleration-Framework/pipelines/fq2sortedbam/
bash install.sh <onprem/cloud>  ## onprem mode: Manually install the depenendies present in basic_setup_ubuntu.sh as it needs sudo access
```

### Run:
```
python run_fq2sortedbam.py --ref <reference> --reads <read1> <read2> --output <outBAMfile>
```

## General Notes:  
1. Individual pipeline tools are present in applications folder    
2. To understand various parameters to these tools, you can access their ```man``` page  
3. You can setup the parameters of these tools using ```params``` commandline option
4. Understand all the parameters to fq3sortedbam using "-h" option 
5. fq2sortedbam supports the following aligners: 
   a. DNA short read alignment using bwa-mem2  
   b. DNA long read alignment using mm2-fast (minimap2)  
   c. RNA short read alignment using STAR aligner  
   d. bwa-meth based alignment  
