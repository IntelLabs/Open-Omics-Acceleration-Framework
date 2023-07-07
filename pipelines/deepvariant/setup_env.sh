#!/bin/bash
# Tested on ubuntu 18.04
echo "Downloading and setting up miniconda..."
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_23.3.1-0-Linux-x86_64.sh
bash ./Miniconda3-py39_23.3.1-0-Linux-x86_64.sh -b -p ./miniconda3
echo "Downloading and setting up miniconda...DONE"

echo "Seeting up conda env named new_env"
miniconda3/bin/conda env create -f environment.yml
echo "Seeting up conda env named new_env...DONE"

echo "Activating conda env..."
source miniconda3/bin/activate new_env

# Dependancy on machine
sudo apt update
echo "deb https://download.opensuse.org/repositories/devel:/kubic:/libcontainers:/stable/xUbuntu_18.04/ /" \
  | sudo tee /etc/apt/sources.list.d/devel:kubic:libcontainers:stable.list
wget -qO - https://download.opensuse.org/repositories/devel:/kubic:/libcontainers:/stable/xUbuntu_18.04/Release.key \
  | gpg --dearmor \
  | sudo tee /etc/apt/trusted.gpg.d/kubic_libcontainers.gpg > /dev/null
sudo apt update
sudo apt -y upgrade
sudo apt -y install podman
sudo apt -y install make
sudo apt -y install autoconf
sudo apt -y install build-essential
sudo apt -y install zlib1g-dev
sudo apt -y install libncurses5-dev
sudo apt -y update
sudo apt -y upgrade
sudo apt -y install libbz2-dev
sudo apt -y install liblzma-dev
