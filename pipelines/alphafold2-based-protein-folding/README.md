
# Build a docker image

### Current docker image requires a dual-socket CPU because it runs multiple inference instances in parallel. However, it can be easily modied to run on a single socket.

```bash
cd ~/Open-Omics-Acceleration-Framework/pipelines/alphafold2/
docker build -t alphafold .           # Build a docker image named alphafold
```
# Preparation 
1. Follow the instructions from https://github.com/deepmind/alphafold repo and download dataset for alphafold2.
2. Create a samples directory that contains fasta files for input proteins. 
3. Create a output directory where model output will be written.
4. Create a log directory where log will be written.
# Run a docker container
```bash
export DATA_DIR=<path-to-database-directory>
export SAMPLES_DIR=<path-to-input-directory>
export OUTPUT_DIR=<path-to-output-directory>
export LOG_DIR=<path-to-log-directory>

docker run -it --cap-add SYS_NICE -v $DATA_DIR:/data \
    -v $SAMPLES_DIR:/samples \
    -v $OUTPUT_DIR:/output \
    -v $LOG_DIR:/open-omics-alphafold/logs \
    alphafold:latest
```

# Running baremetal

To run the optimized alphafold2 without docker (baremetal)
1. Clone the open-omics-alphafold submodule present in the applications directory of Open-Omics-Acceleration-Framework.
2. Follow the readme instructions of the submodule for creating conda environment and runnning inference.
