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
miniconda3/bin/conda env create --name $1 -f environment.yml
echo "Seeting up conda env named new_env...DONE"

echo "Activating conda env..."
source miniconda3/bin/activate $1

# Dependancy on machine
############# docker installation ############
sudo apt update
sudo apt-get install ca-certificates curl gnupg
sudo install -m 0755 -d /etc/apt/keyrings
curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /etc/apt/keyrings/docker.gpg
sudo chmod a+r /etc/apt/keyrings/docker.gpg
echo \
  "deb [arch="$(dpkg --print-architecture)" signed-by=/etc/apt/keyrings/docker.gpg] https://download.docker.com/linux/ubuntu \
  "$(. /etc/os-release && echo "$VERSION_CODENAME")" stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
sudo apt-get install docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin
#To verify docker installation use below command
#sudo docker run hello-world 
##############################################


sudo apt update
sudo apt -y upgrade
sudo apt -y install make
sudo apt -y install autoconf
sudo apt -y install build-essential
sudo apt -y install zlib1g-dev
sudo apt -y install libncurses5-dev
sudo apt -y update
sudo apt -y upgrade
sudo apt -y install libbz2-dev
sudo apt -y install liblzma-dev

# This will save deepvariant images
cd ${ABS_DIRECTORY}applications/deepvariant
docker build -t deepvariant .
#save image(~7 GB) to tar file if you are using multiple nodes.
cd ${ABS_DIRECTORY}/pipelines/deepvariant/
dcoker save -o deepvariant.tar deepvariant:latest   #get image id using 'podman images'
