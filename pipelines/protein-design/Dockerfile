# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

# Use a base image for Python and Conda
ARG BASE_IMAGE=condaforge/miniforge3:24.3.0-0
FROM ${BASE_IMAGE} as conda_setup
ENV DEBIAN_FRONTEND=noninteractive

# Stage 2: Set up the main build environment
FROM ubuntu:24.04 as builder
ENV DEBIAN_FRONTEND=noninteractive

# Install necessary build tools and clean up
RUN apt-get update && apt-get install -y --no-install-recommends \
    git build-essential wget vim ca-certificates numactl autoconf automake make curl && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get autoremove -y && \
    apt-get clean

# Non-root user setup
ENV SERVICE_NAME="colab-design"
RUN groupadd --gid 1001 $SERVICE_NAME && \
    useradd -m -g $SERVICE_NAME --shell /bin/false --uid 1001 $SERVICE_NAME

# Copy Conda installation from the conda_setup stage
COPY --from=conda_setup /opt/conda /opt/conda
ENV PATH="/opt/conda/bin:$PATH"

# Set working directory
WORKDIR /app

RUN git clone https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
RUN /app/Open-Omics-Acceleration-Framework
RUN git checkout db29aec9c3e2eb27c36dd91824dfc346e5deae89
#RUN git switch protein_design

# Clone the ColabDesign repository
WORKDIR /app
RUN git clone https://github.com/sokrypton/ColabDesign.git

WORKDIR /app/ColabDesign
RUN git checkout 4c0bc6d67f8f967135ecccc135a26b3bfded25e8
RUN cp /app/Open-Omics-Acceleration-Framework/pipelines/protein-design/alphafold_pre.patch .
RUN patch -p1 --force < alphafold_pre.patch || true

RUN rm -rf /app/Open-Omics-Acceleration-Framework
# Create the destination directory
RUN mkdir -p /app/ColabDesign/params

# Step 1: Download the file first
RUN curl -fLo /app/ColabDesign/params/alphafold_params_2022-12-06.tar \
    https://storage.googleapis.com/alphafold/alphafold_params_2022-12-06.tar

# Step 2: Extract the file
RUN tar -xf /app/ColabDesign/params/alphafold_params_2022-12-06.tar \
    -C /app/ColabDesign/params

WORKDIR /app/ColabDesign

# Create a new Conda environment and install dependencies using Python 3.11
RUN conda create -n colab_env python=3.11 -y && \
    echo "source activate colab_env" >> ~/.bashrc

# Ensure Python uses the correct environment
#ENV PATH="/opt/conda/envs/colab_env/bin:$PATH"
ENV PATH="/opt/conda/envs/colab_env/bin:/opt/conda/bin:$PATH"

#RUN source ~/.bashrc

RUN /opt/conda/envs/colab_env/bin/python -m pip install torch==2.6.0+cpu torchvision==0.21.0+cpu torchaudio==2.6.0+cpu --index-url https://download.pytorch.org/whl/cpu

RUN pip install -e .

# Change ownership of the directory
RUN chown -R $SERVICE_NAME:$SERVICE_NAME /app

# Switch to non-root user
USER $SERVICE_NAME

# Healthcheck disabled
HEALTHCHECK NONE

CMD ["python", "designability_test.py"]
