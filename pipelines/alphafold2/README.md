
# Build a docker image
```bash
cd ~/Open-Omics-Acceleration-Framework/pipelines/alphafold2/
docker build -t alphafold .           # Build a docker image named alphafold
```
# Preparation 
1. Follow the instructions from https://github.com/deepmind/alphafold repo and download dataset for alphafold2.
2. Create a samples directory that contains fasta files for input protiens. 
3. Create a output directory where model output will be written.

# Run a docker container
```bash
export DATA_DIR=<path-to-database-directory>
export SAMPLES_DIR=<path-to-input-directory>
export OUTPUT_DIR=<path-to-output-directory>

docker run -it -v $DATA_DIR:/data \
    -v $SAMPLES_DIR:/samples \
    -v $OUTPUT_DIR:/output \
    alphafold:latest
```

# Running baremetal

To run the optimized alphafold2 without docker (baremetal)
1. Clone the open-omics-alphafold submodule present in the applications directory of Open-Omics-Acceleration-Framework.
2. Follow the readme instructions of the submodule for creating conda environment and runnning inference.
