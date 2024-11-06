# fq2vcf: OpenOmics Deepvariant based Variant Calling Pipeline
## Overview:
    OpenOmics's fq2vcf is a highly optimized, distributed, deep learning-based short-read germline variant calling pipeline for x86 CPUs. The pipeline comprises of:
    1. bwa-mem2 (a highly optimized version of bwa-mem) for sequence mapping
    2. SortSAM using samtools
    3. An optimized version of DeepVariant tool for Variant Calling


## Example run
    The following script runs fq2vcf pipeline on human genome dataset with small number of reads preset in ~/data/fq2vcf/ folder
    ```cd ~/bin/fq2vcf/
    bash fq2vcf sh````
    Runtime on m7i.48xlarge w/ default AMI hardware: < 600 sec
    You can follow the below steps to run the pipeline on your own input data (reference genome and read(s) files)


## Source code
``` https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/tree/main/pipelines/deepvariant-based-germline-variant-calling-fq2vcf ```


## Input parameter setup
    <readdir>: User provided location of the local directory containing read files R1 & R2
    <refdir> : User provided location of the local directory containing reference sequence file REF
    <outdir> : User provided location of the local directory for output files

### Setup config files
    Update the following fields in the config file at:
    1.  ~/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/config
    2.  ~/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/extra_scripts/config
    
    Update the following three fields in the above two config files:
    R1=HG001_1.fastq.gz ## name of input reads file1
    R2=HG001_2.fastq.gz ## name of input reads file2
    REF=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna ## name of reference sequence (should not be in .gz format)


## Index creation
Create bwa-mem2 index and reference sequence index (one-time task)
```
docker run -v ~/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/config:/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/scripts/aws/config  -v <refdir>:/ref  -it deepvariant:part1 bash create_reference_index.sh
```

## Run:
Run the pipeline using the following two steps:

```
docker run -v ~/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/config:/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/scripts/aws/config -v <readdir>:/reads -v <refdir>:/ref -v <outdir>:/output -it deepvariant:part1 bash run_pipeline_ec2_part1.sh
```

```
docker run -v ~/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/extra_scripts/config:/opt/deepvariant/config  -v <refdir>:/ref -v <outdir>:/output -it deepvariant:part2 bash run_pipeline_ec2_part2.sh
```
