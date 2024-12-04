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
    git autoconf build-essential wget vim ca-certificates numactl && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get autoremove -y && \
    apt-get clean

# Non-root user setup
ENV SERVICE_NAME="proteinmpnn-service"
RUN groupadd --gid 1001 $SERVICE_NAME && \
    useradd -m -g $SERVICE_NAME --shell /bin/false --uid 1001 $SERVICE_NAME

WORKDIR /
RUN git clone https://github.com/dauparas/ProteinMPNN.git
WORKDIR /ProteinMPNN
RUN git checkout 8907e6671bfbfc92303b5f79c4b5e6ce47cdef57

# Apply the patch file
RUN wget https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/releases/download/3.0/Source_code_with_submodules.tar.gz 
RUN tar -xzf Source_code_with_submodules.tar.gz
RUN cp /ProteinMPNN/Open-Omics-Acceleration-Framework/applications/ProteinMPNN/ProteinMPNN.patch .
RUN git apply ProteinMPNN.patch
RUN rm -rf Open-Omics-Acceleration-Framework
RUN rm -rf Source_code_with_submodules.tar.gz
ENV PATH="/opt/conda/bin:$PATH"

# Copy Conda installation from the conda_setup stage
COPY --from=conda_setup /opt/conda /opt/conda

RUN conda create -n p_mpnn python=3.11 pip=24.0

# Install PyTorch, Torchvision, Torchaudio
RUN conda install -n p_mpnn -y pytorch==2.3.1 torchvision==0.18.1 torchaudio==2.3.1 -c pytorch
# Install Intel PyTorch extension
RUN /opt/conda/envs/p_mpnn/bin/python -m pip install intel-extension-for-pytorch==2.3.100 oneccl-bind-pt==2.3.0 --extra-index-url https://pytorch-extension.intel.com/release-whl/stable/cpu/us/
RUN /opt/conda/envs/p_mpnn/bin/python -m pip install numpy==1.26.0

# Clone jemalloc source for building
WORKDIR /ProteinMPNN
RUN git clone --branch 5.3.0 https://github.com/jemalloc/jemalloc.git
WORKDIR /ProteinMPNN/jemalloc
RUN bash autogen.sh --prefix=/opt/conda/envs/p_mpnn/ && make install
WORKDIR /ProteinMPNN
RUN rm -rf jemalloc

# Set up environment activation and PATH
ENV PATH="/opt/conda/envs/p_mpnn/bin:/opt/conda/bin:$PATH"

# Ensure all scripts inside ProteinMPNN/examples are executable
RUN chmod +x /ProteinMPNN/examples/*.py
RUN mkdir /outputs

# Change ownership of the directory
RUN chown -R $SERVICE_NAME:$SERVICE_NAME /ProteinMPNN /outputs

# Switch to non-root user
USER $SERVICE_NAME

# Healthcheck disabled
HEALTHCHECK NONE

# Set environment variables for jemalloc
ENV LD_PRELOAD "/opt/conda/envs/p_mpnn/lib/libjemalloc.so:$LD_PRELOAD"
ENV MALLOC_CONF "oversize_threshold:1,background_thread:true,metadata_thp:auto,dirty_decay_ms:-1,muzzy_decay_ms:-1"
# Set the default command to run the application
CMD source ~/.bashrc && ["/bin/bash", "python"] 
