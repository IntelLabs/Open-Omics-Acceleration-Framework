# AlphaFold2: OpenOmics' Accelerated AlphaFold2 
## Overview: 
Given one or more protein sequences, this workflow performs preprocessing (database search and multiple sequence alignment using Open Omics HMMER and HH-suite) and structure prediction through AlphaFold2's Evoformer model (Open Omics AlphaFold2) to output the structure(s) of the protein sequences. 


## Example run
```
cd ~/bin/alphafold2/
bash alphafold2.sh
```
Runtime on m7i.48xlarge w/ default AMI hardware: < 2000 sec
Follow the below instructions if you would like to run the pipeline on your own data

## Source code
``` https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/tree/main/pipelines/alphafold2-based-protein-folding ```


## AlphaFold2 Protein Database
AlphaFold database is location at /af2data directory

## Run with your own input data
export DATA_DIR=/af2data/alphafold_database/
export SAMPLES_DIR=<path-to-input-directory>
export OUTPUT_DIR=<path-to-output-directory>
export LOG_DIR=<path-to-log-directory>


### Run pre-processign step
docker run -it --cap-add SYS_NICE -v $DATA_DIR:/data \
    -v $SAMPLES_DIR:/samples \
    -v $OUTPUT_DIR:/output \
    -v $LOG_DIR:/Open-Omics-Acceleration-Framework/applications/alphafold/logs \
    alphafold:pre

### Run inference step
docker run -it --cap-add SYS_NICE -v $DATA_DIR:/data \
    -v $SAMPLES_DIR:/samples \
    -v $OUTPUT_DIR:/output \
    -v $LOG_DIR:/Open-Omics-Acceleration-Framework/applications/alphafold/logs \
    alphafold:inf
