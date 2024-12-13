# Pipeline overview
Given one or more protein sequences, this workflow performs preprocessing (database search and multiple sequence alignment using Open Omics [HMMER](https://github.com/IntelLabs/hmmer) and [HH-suite](https://github.com/IntelLabs/hh-suite)) and structure prediction through AlphaFold2's Evoformer model ([Open Omics AlphaFold2](https://github.com/IntelLabs/open-omics-alphafold)) to output the structure(s) of the protein sequences. The following block diagram illustrates the pipeline.

<p align="center">
<img src="https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/blob/main/images/alphafold2-protein-folding.jpg"/a></br>
</p> 

# Build a docker image

### Current docker image requires a single socket/dual-socket CPU with 1 or 2 NUMA domains, because it runs multiple inference instances in parallel. It can be easily modified to run at other types of machines.

```bash
cd ~/Open-Omics-Acceleration-Framework/pipelines/alphafold2-based-protein-folding
docker build -t alphafold:pre -f Dockerfile_Pre .   # Build a docker image named alphafold:pre for pre-processing step
docker build -t alphafold:inf -f Dockerfile_Inf .   # Build a docker image named alphafold:inf for inference step

```
# Preparation 
1. Follow the instructions from https://github.com/deepmind/alphafold repo and download the database for alphafold2.
2. Create a samples directory that contains fasta files for input proteins. 
3. Create a output directory where model output will be written.
4. Create a log directory where log will be written.
# Run a docker container
```bash
export DATA_DIR=<path-to-database-directory>
export SAMPLES_DIR=<path-to-input-directory>
export OUTPUT_DIR=<path-to-output-directory>
export LOG_DIR=<path-to-log-directory>


# Run pre-processign step for monomer
docker run -it --cap-add SYS_NICE -v $DATA_DIR:/data \
    -v $SAMPLES_DIR:/samples \
    -v $OUTPUT_DIR:/output \
    -v $LOG_DIR:/Open-Omics-Acceleration-Framework/applications/alphafold/logs \
    alphafold:pre

# Run pre-processign step for multimer
docker run -it --cap-add SYS_NICE -v $DATA_DIR:/data \
    -v $SAMPLES_DIR:/samples \
    -v $OUTPUT_DIR:/output \
    -v $LOG_DIR:/Open-Omics-Acceleration-Framework/applications/alphafold/logs \
    alphafold:pre multimer

# Run inference step for monomer with relexation
docker run -it --cap-add SYS_NICE -v $DATA_DIR:/data \
    -v $SAMPLES_DIR:/samples \
    -v $OUTPUT_DIR:/output \
    -v $LOG_DIR:/Open-Omics-Acceleration-Framework/applications/alphafold/logs \
    alphafold:inf monomer relax

# Run inference step for multimer with relexation
docker run -it --cap-add SYS_NICE -v $DATA_DIR:/data \
    -v $SAMPLES_DIR:/samples \
    -v $OUTPUT_DIR:/output \
    -v $LOG_DIR:/Open-Omics-Acceleration-Framework/applications/alphafold/logs \
    alphafold:inf multimer relax
```

# Running baremetal

To run the optimized alphafold2 without docker (baremetal)
1. Clone the open-omics-alphafold submodule present in the applications directory of this repo.
2. Follow the readme instructions of the submodule for creating conda environment and runnning inference.
