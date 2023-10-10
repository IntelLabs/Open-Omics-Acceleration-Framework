#!/bin/bash

SCRIPT_PATH="${BASH_SOURCE:-$0}"
ABS_SCRIPT_PATH="$(realpath "${SCRIPT_PATH}")"
#echo "Value of ABS_SCRIPT_PATH: ${ABS_SCRIPT_PATH}"
ABS_DIRECTORY="$(dirname "${ABS_SCRIPT_PATH}")"
#echo "Value of ABS_DIRECTORY: ${ABS_DIRECTORY}"


echo "Downloading and setting up miniconda..."
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_23.3.1-0-Linux-x86_64.sh
bash ./Miniconda3-py39_23.3.1-0-Linux-x86_64.sh -b -p ./miniconda3
echo "Downloading and setting up miniconda...DONE"

echo "Seeting up conda env named with given argument"
/data/nfs_home/mvasimud/container/miniconda3/bin/conda env create --name bwa -f environment.yml
echo "Seeting up conda env named new_env...DONE"

echo "Activating conda env..."
source /data/nfs_home/mvasimud/container/miniconda3/bin/activate bwa
echo "localhost" > hostfile

