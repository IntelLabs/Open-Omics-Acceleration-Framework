# Use the official Ubuntu image as a base
ARG FROM_IMAGE=ubuntu:24.04

# Stage 1: Set up Conda environment
ARG BASE_IMAGE=condaforge/miniforge3:24.3.0-0
FROM ${BASE_IMAGE} as conda_setup

ENV DEBIAN_FRONTEND=noninteractive

# Stage 2: Set up the main build environment
FROM ${FROM_IMAGE} as builder

ENV DEBIAN_FRONTEND=noninteractive

ARG http_proxy
ENV http_proxy=${http_proxy}

ARG https_proxy
ENV https_proxy=${https_proxy}

ARG no_proxy
ENV no_proxy=${no_proxy}

# Install necessary build tools and clean up
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates git build-essential vim numactl autoconf automake make && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /app
# Non-root user setup
ENV SERVICE_NAME="protgpt2-service"

RUN groupadd --gid 1001 $SERVICE_NAME && \
    useradd -m -g $SERVICE_NAME --shell /bin/false --uid 1001 $SERVICE_NAME

# Copy Conda installation from the conda_setup stage
COPY --from=conda_setup /opt/conda /opt/conda

RUN wget https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/releases/download/3.0/Source_code_with_submodules.tar.gz 
RUN tar -xzf Source_code_with_submodules.tar.gz

ENV PATH="/opt/conda/bin:$PATH"

# Copy environment.yml and create Conda environment
RUN cp /app/Open-Omics-Acceleration-Framework/applications/ProtGPT2/env.yml .
RUN cp /app/Open-Omics-Acceleration-Framework/applications/ProtGPT2/protgpt2.py .
RUN cp /app/Open-Omics-Acceleration-Framework/applications/ProtGPT2/model_script.sh .
RUN rm -rf Open-Omics-Acceleration-Framework
RUN rm -rf Source_code_with_submodules.tar.gz
RUN conda env create -f env.yml

RUN git clone --branch 5.3.0 https://github.com/jemalloc/jemalloc.git
WORKDIR /app/jemalloc
RUN bash autogen.sh --prefix=/opt/conda/envs/protgpt2/ && make install
WORKDIR /app
RUN rm -rf jemalloc

# Set up environment activation and PATH
ENV PATH="/opt/conda/envs/protgpt2/bin:$PATH"

# Swith to Non-root user
RUN chown -R $SERVICE_NAME:$SERVICE_NAME /app
USER $SERVICE_NAME

HEALTHCHECK NONE

ENV LD_PRELOAD "/opt/conda/envs/protgpt2/lib/libjemalloc.so:$LD_PRELOAD" 
ENV MALLOC_CONF "oversize_threshold:1,background_thread:true,metadata_thp:auto,dirty_decay_ms:-1,muzzy_decay_ms:-1" 
CMD python
