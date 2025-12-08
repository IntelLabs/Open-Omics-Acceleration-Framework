# Stage 1: Conda base
ARG BASE_IMAGE=condaforge/miniforge3:24.3.0-0
FROM ${BASE_IMAGE} as conda_setup
ENV DEBIAN_FRONTEND=noninteractive

# Stage 2: Build environment
FROM ubuntu:24.04 as builder
ENV DEBIAN_FRONTEND=noninteractive

# Install base tools
RUN apt-get update && apt-get install -y --no-install-recommends \
    git build-essential wget vim ca-certificates numactl autoconf automake make curl patch && \
    rm -rf /var/lib/apt/lists/*

# Add non-root user
ENV SERVICE_NAME="proteinmpnn-service"
RUN groupadd --gid 1001 $SERVICE_NAME && \
    useradd -m -g $SERVICE_NAME --shell /bin/bash --uid 1001 $SERVICE_NAME

# Copy Conda from setup stage
COPY --from=conda_setup /opt/conda /opt/conda
ENV PATH="/opt/conda/bin:$PATH"

WORKDIR /

RUN git clone https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
WORKDIR /Open-Omics-Acceleration-Framework
RUN git switch db29aec9c3e2eb27c36dd91824dfc346e5deae89
#RUN git switch protein_design
#RUN mv Open-Omics-Acceleration-Framework /OpenOmics


# ---------------------------
# Apply patch to ProteinMPNN
# ---------------------------
WORKDIR /
RUN git clone https://github.com/dauparas/ProteinMPNN.git
WORKDIR /ProteinMPNN
RUN git checkout 8907e6671bfbfc92303b5f79c4b5e6ce47cdef57
RUN cp /Open-Omics-Acceleration-Framework/applications/ProteinMPNN/ProteinMPNN_torch.patch .
RUN git apply ProteinMPNN_torch.patch

# ---------------------------
# Microservice copy
# ---------------------------
#WORKDIR /
#RUN cp -r /OpenOmics/applications/ProteinMPNN/microservice /ProteinMPNN/

# ---------------------------
# Clone and Patch ColabDesign
# ---------------------------
WORKDIR /
RUN git clone https://github.com/sokrypton/ColabDesign.git

WORKDIR /ColabDesign
RUN git checkout 4c0bc6d67f8f967135ecccc135a26b3bfded25e8
# Apply JAX patch from the same Open-Omics folder
RUN cp /Open-Omics-Acceleration-Framework/applications/ProteinMPNN/ProteinMPNN_jax.patch .
RUN patch -p1 --force < ProteinMPNN_jax.patch || true

# ---------------------------
# Cleanup
# ---------------------------
WORKDIR /
RUN rm -rf /Open-Omics-Acceleration-Framework

# Create shared conda env
RUN conda create -n p_mpnn python=3.11 -y && \
    echo "source activate p_mpnn" >> ~/.bashrc

ENV PATH="/opt/conda/envs/p_mpnn/bin:/opt/conda/bin:$PATH"

# Install requirements via pip (from setup.py)
WORKDIR /ColabDesign
RUN pip install .

# Install torch in the same env
RUN conda install -n p_mpnn -y pytorch==2.3.1 torchvision==0.18.1 torchaudio==2.3.1 -c pytorch

WORKDIR /ProteinMPNN
RUN git clone --branch 5.3.0 https://github.com/jemalloc/jemalloc.git
WORKDIR /ProteinMPNN/jemalloc
RUN bash autogen.sh --prefix=/opt/conda/envs/p_mpnn/ && make install
WORKDIR /ProteinMPNN
RUN rm -rf jemalloc
# Make MPNN scripts executable
RUN chmod +x /ProteinMPNN/examples/*.py

RUN pip install opea-comps==1.3 && \
    pip install numpy==1.26.4 && \
    pip install netifaces==0.11.0

# Prepare output dir
RUN mkdir /outputs && chmod 777 /outputs

# Change ownership
RUN chown -R $SERVICE_NAME:$SERVICE_NAME /ProteinMPNN /outputs

# Change ownership of the directory
RUN chown -R $SERVICE_NAME:$SERVICE_NAME /ColabDesign

# Switch to non-root user
USER $SERVICE_NAME

# Healthcheck disabled
HEALTHCHECK NONE

# Set environment variables for jemalloc
ENV LD_PRELOAD "/opt/conda/envs/p_mpnn/lib/libjemalloc.so:$LD_PRELOAD"
ENV MALLOC_CONF "oversize_threshold:1,background_thread:true,metadata_thp:auto,dirty_decay_ms:-1,muzzy_decay_ms:-1"

CMD ["python","designability_test.py","opea_pmpnn_client.py"]
