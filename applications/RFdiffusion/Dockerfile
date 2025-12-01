# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

# Use the official Ubuntu image as a base
ARG FROM_IMAGE=ubuntu:24.04

# Stage 1: Set up Conda environment
ARG BASE_IMAGE=condaforge/miniforge3:24.3.0-0
FROM ${BASE_IMAGE} as conda_setup
ENV DEBIAN_FRONTEND=noninteractive

# Stage 2: Set up the main build environment
FROM ${FROM_IMAGE} as builder
ENV DEBIAN_FRONTEND=noninteractive

# Install necessary build tools and clean up
RUN apt-get update && apt-get install -y --no-install-recommends \
    git build-essential wget vim ca-certificates numactl autoconf automake make && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get autoremove -y && \
    apt-get clean

# Non-root user setup
ENV SERVICE_NAME="rfdiffusion-service"
RUN groupadd --gid 1001 $SERVICE_NAME && \
    useradd -m -g $SERVICE_NAME --shell /bin/false --uid 1001 $SERVICE_NAME 

# Copy Conda installation from the conda_setup stage
COPY --from=conda_setup /opt/conda /opt/conda

WORKDIR /app
RUN git clone https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
WORKDIR /app/Open-Omics-Acceleration-Framework
RUN git switch db29aec9c3e2eb27c36dd91824dfc346e5deae89
#RUN git switch protein_design

WORKDIR /app
RUN git clone https://github.com/RosettaCommons/RFdiffusion.git 

WORKDIR /app/RFdiffusion
# adding the git commit id
RUN git checkout 820bfdfaded8c260b962dc40a3171eae316b6ce0

# ----------------------------
# Apply patch to RFdiffusion
# ----------------------------
RUN cp /app/Open-Omics-Acceleration-Framework/applications/RFdiffusion/RFdiffusion.patch .
RUN git apply RFdiffusion.patch

WORKDIR /app/RFdiffusion/models
RUN  wget https://files.ipd.uw.edu/pub/RFdiffusion/6f5902ac237024bdd0c176cb93063dc4/Base_ckpt.pt && \
     wget https://files.ipd.uw.edu/pub/RFdiffusion/e29311f6f1bf1af907f9ef9f44b8328b/Complex_base_ckpt.pt && \
     wget https://files.ipd.uw.edu/pub/RFdiffusion/60f09a193fb5e5ccdc4980417708dbab/Complex_Fold_base_ckpt.pt && \
     wget https://files.ipd.uw.edu/pub/RFdiffusion/74f51cfb8b440f50d70878e05361d8f0/InpaintSeq_ckpt.pt && \
     wget https://files.ipd.uw.edu/pub/RFdiffusion/76d00716416567174cdb7ca96e208296/InpaintSeq_Fold_ckpt.pt && \
     wget https://files.ipd.uw.edu/pub/RFdiffusion/5532d2e1f3a4738decd58b19d633b3c3/ActiveSite_ckpt.pt && \
     wget https://files.ipd.uw.edu/pub/RFdiffusion/12fc204edeae5b57713c5ad7dcb97d39/Base_epoch8_ckpt.pt && \
     wget https://files.ipd.uw.edu/pub/RFdiffusion/f572d396fae9206628714fb2ce00f72e/Complex_beta_ckpt.pt && \
     wget https://files.ipd.uw.edu/pub/RFdiffusion/1befcb9b28e2f778f53d47f18b7597fa/RF_structure_prediction_weights.pt

WORKDIR /app/RFdiffusion

RUN /opt/conda/bin/conda env create -f env/SE3nv.yml

RUN git clone --branch 5.3.0 https://github.com/jemalloc/jemalloc.git
WORKDIR /app/RFdiffusion/jemalloc
RUN bash autogen.sh --prefix=/opt/conda/envs/SE3nv/ && make install
WORKDIR /app/RFdiffusion
RUN rm -rf jemalloc

# Set up environment activation and PATH
ENV PATH="/opt/conda/envs/SE3nv/bin:/opt/conda/bin:$PATH"
	
# Install dependencies
WORKDIR /app/RFdiffusion/env/SE3Transformer
RUN pip install --no-cache-dir -r requirements.txt
RUN python setup.py install

WORKDIR /app/RFdiffusion
RUN pip install -e .

RUN tar -xvf examples/ppi_scaffolds_subset.tar.gz -C examples/

# Install ipywidgets
RUN pip install ipywidgets

WORKDIR /app/RFdiffusion
RUN git clone https://github.com/sokrypton/ColabDesign.git
WORKDIR /app/RFdiffusion/ColabDesign
RUN git checkout 4c0bc6d67f8f967135ecccc135a26b3bfded25e8
RUN cp /app/Open-Omics-Acceleration-Framework/applications/RFdiffusion/pipeline_design_colab.patch . 
RUN git apply pipeline_design_colab.patch

#RUN patch -p1 < pipeline_design_colab.patch || true
RUN pip install -e .

RUN pip install opea-comps==1.3 && \
    pip install numpy==1.26.4 && \
    pip install netifaces==0.11.0 

WORKDIR /app/RFdiffusion/scripts
RUN cp /app/Open-Omics-Acceleration-Framework/applications/RFdiffusion/postprocess.py .
RUN rm -rf /app/Open-Omics-Acceleration-Framework
RUN chown -R $SERVICE_NAME:$SERVICE_NAME /app

# Switch to non-root user
USER $SERVICE_NAME
HEALTHCHECK NONE

ENV LD_PRELOAD "/opt/conda/envs/SE3nv/lib/libjemalloc.so:$LD_PRELOAD"
ENV MALLOC_CONF "oversize_threshold:1,background_thread:true,metadata_thp:auto,dirty_decay_ms:-1,muzzy_decay_ms:-1" 
CMD python 
