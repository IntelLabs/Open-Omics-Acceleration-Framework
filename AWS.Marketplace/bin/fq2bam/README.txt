# fq2SortedBAM: OpenOmics' Genomics Secondary Analysis Pipeline
## Overview:
    The pipeline takes input fastq files and produces sorted BAM file through the following stages:
    Sequence Alignment: bwa-mem2 for short reads, mm2-fast (Accelerated Minimap2) for long reads (PacBio, ONT)
    SAMSort (Using SAMTools)

## Example run:
    The following script runs fq2sortedbam pipeline on human genome dataset with small number of reads preset in ~/data/fq2bam/ folder
    ```cd ~/bin/fq2bam/
    bash fq2bam.sh````
    Runtime on m7i.48xlarge w/ default AMI hardware: < 600 sec

    You can follow the below steps to run the pipeline on your own input data (reference genome and read(s) files)

## Source code
``` https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/tree/main/pipelines/fq2sortedbam ````


## Run Modes:
    fq2SortedBAM supports 4 different modes:
    1. sortedbam: It takes fastq reads files and reference genome as input and outputs sorted BAM file
    2. flatmode: It takes fastq reads files and reference genome as input, and outputs multiple (equal to the number of ranks created) unsorted SAM files


## Input parameter setup:
### IO folder locations
    <inputdir>: Location of the local directory containing read files read1 & read2
    <refdir>  : Location of the local directory containing reference sequence file ref
    <outdir>  : Location of the local directory for output files SAM/BAM
    <tempdir> : Location of the local directory for temporary files (defaults to <outdir>)

### Copy config.yaml
    ```cp ~/Open-Omics-Acceleration-Framework/pipelines/fq2sortedbam/config.yaml <inputdir>```

### Setup config.yaml at (<inputdir>/config.yaml)
bwa: bwa-mem2 related parameters
    params: dtype=string, the command line paramteres to bwa-mem2 mapping run e.g. '+R "@RG\tID:RG1\tSM:RGSN1"'
    rindex: dtype=bool, 'values=True/False', if "True" it creates bwa-mem2 index for the reference genome. Use this option if you don't have reference already created.

dataset: data details
    index: dtype:string, Input reference genome file name e.g. "GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" -- should be present in <refdir>
    outfile: dtype=string, output SAM/BAM file name(s) e.g. short.se.sam. Default value: "final_fq2bam"  -- will be created in <outputdir>
    read1: dtype=string, input reads file1 name e.g. "HG001.novaseq.pcr-free.30x.R1.fastq.gz"  --  should be present in <inputdir>
    read2: dtype=string, input reads file2 name e.g. "HG001.novaseq.pcr-free.30x.R2.fastq.gz"  --  should be present in <inputdir>
    read_type: dtype:string, values=short/long, short read mapping using bwa-mem2, long read mapping using mm2-fast


## Docker run:
```docker run -v <inputdir>:/input -v <outdir>:/out -v <refdir>:/refdir -v <tempdir>:/tempdir fq2bam:latest bash run_bwa.sh sortedbam /input/config.yaml```


