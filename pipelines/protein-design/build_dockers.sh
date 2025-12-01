# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License
# Author: Rahamathullah
# Email shaikx.rahamathullah@intel.com

#!/bin/bash

set -e  # Stop if any command fails
git clone https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
cd Open-Omics-Acceleration-Framework/
git checkout db29aec9c3e2eb27c36dd91824dfc346e5deae89

cd applications/RFdiffusion

# Build RFdiffusion image
echo "Running RFdiffusion..."
docker build --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy --build-arg no_proxy="127.0.0.1,localhost,apt.repo.inel.com" -t rfdiffusion:latest .

cd ../
cd ProteinMPNN/
# Build ProteinMPNN image
echo "Running proteinmpnn...."
docker build --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy --build-arg no_proxy="127.0.0.1,localhost,apt.repo.inel.com" -t pmpnn:latest .

cd ../../
cd pipelines/protein-design/

# Build Colab Alphafold pre image
echo "Running colab alphafold pre..."
docker build --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy --build-arg no_proxy="127.0.0.1,localhost,apt.repo.inel.com" -t colab_af:pre .

cd ../
cd alphafold2-based-protein-folding
echo "Running alphafold preprocess..."
docker build --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy --build-arg no_proxy="127.0.0.1,localhost,apt.repo.inel.com" -t alphafold:pre -f Dockerfile_Pre .

echo "Running alphafold inference..."
docker build --build-arg http_proxy=$http_proxy --build-arg https_proxy=$https_proxy --build-arg no_proxy="127.0.0.1,localhost,apt.repo.inel.com" -t alphafold:inf -f Dockerfile_Inf .

echo "All Docker images built successfully!"
