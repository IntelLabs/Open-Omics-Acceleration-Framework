#sudo mkdir /af2data
#sudo mount /dev/nvme1n1 /af2data



export DATA_DIR=/af2data/alphafold_database/
export SAMPLES_DIR=~/data/alphafold2/input/
export OUTPUT_DIR=~/data/alphafold2/output/
export LOG_DIR=~/data/alphafold2/logs/


# Run pre-processign step
docker run -it --cap-add SYS_NICE -v $DATA_DIR:/data \
    -v $SAMPLES_DIR:/samples \
    -v $OUTPUT_DIR:/output \
    -v $LOG_DIR:/Open-Omics-Acceleration-Framework/applications/alphafold/logs \
    alphafold:pre

# Run inference step
docker run -it --cap-add SYS_NICE -v $DATA_DIR:/data \
    -v $SAMPLES_DIR:/samples \
    -v $OUTPUT_DIR:/output \
    -v $LOG_DIR:/Open-Omics-Acceleration-Framework/applications/alphafold/logs \
    alphafold:inf
