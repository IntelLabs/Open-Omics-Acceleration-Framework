# Copyright 2025 Intel Corporation
# SPDX-License-Identifier: MIT License

# Base image specification
ARG FROM_IMAGE=ubuntu:24.04
FROM ${FROM_IMAGE} as builder

ENV DEBIAN_FRONTEND=noninteractive

# Proxy configuration
ARG http_proxy
ENV http_proxy=${http_proxy}

ARG https_proxy
ENV https_proxy=${https_proxy}

ARG no_proxy
ENV no_proxy=${no_proxy}

# Install necessary build tools and clean up
RUN apt-get update && apt-get install -y --no-install-recommends \
    git build-essential wget vim ca-certificates autoconf automake make numactl unzip time iproute2 && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get autoremove -y && \
    apt-get clean

# Create a user and group for running the service
ENV SERVICE_NAME="moflow-base-service"
RUN groupadd --gid 1001 $SERVICE_NAME && \
    useradd -m -g $SERVICE_NAME --shell /bin/false --uid 1001 $SERVICE_NAME

WORKDIR /app
RUN chown -R $SERVICE_NAME:$SERVICE_NAME /app
USER $SERVICE_NAME

# Install Miniforge (Conda)
RUN wget --no-check-certificate -O Miniforge3.sh "https://github.com/conda-forge/miniforge/releases/download/24.3.0-0/Miniforge3-$(uname)-$(uname -m).sh"
RUN bash Miniforge3.sh -b -p "${HOME}/conda" && rm Miniforge3.sh
COPY . /app

# Clone MoFlow repository and apply patch
RUN git clone --recursive https://github.com/calvin-zcx/moflow.git moflow
WORKDIR /app/moflow
RUN git checkout 3611c637260272b3d34a298f221623cb59e01091

RUN git clone https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
WORKDIR /app/moflow/Open-Omics-Acceleration-Framework/
RUN git checkout moflow_v4

#RUN wget https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/releases/download/4.0/Source_code_with_submodules.tar.gz
#RUN tar -xzf Source_code_with_submodules.tar.gz
#WORKDIR /app/moflow/Open-Omics-Acceleration-Framework/

USER root
WORKDIR /app/moflow/
RUN mv /app/moflow/Open-Omics-Acceleration-Framework/applications/moflow/* /app/moflow/ && \
    mv /app/moflow/Open-Omics-Acceleration-Framework/applications/common /app/moflow/ && \
    rm -rf /app/moflow/Open-Omics-Acceleration-Framework

RUN git apply /app/moflow/docker_setup/mflow_change_all.patch
USER $SERVICE_NAME

# Set up Conda environment
RUN ${HOME}/conda/bin/mamba env create -f /app/moflow/docker_setup/env.yml

# Install MoFlow package
WORKDIR /app/moflow
RUN bash -c "source ${HOME}/conda/etc/profile.d/conda.sh && \
             source ${HOME}/conda/etc/profile.d/mamba.sh && \
             mamba activate moflow && \
             pip install torch==2.3.0 torchvision==0.18.0 torchaudio==2.3.0 --index-url https://download.pytorch.org/whl/cpu  && \
             pip install intel-extension-for-pytorch==2.3.0 && \
             pip install opea-comps==1.3 && \
             pip install absl-py"

# Initialization script for Conda environment
RUN echo "#!/bin/bash" > /app/init.sh && \
    echo "source ${HOME}/conda/etc/profile.d/conda.sh" >> /app/init.sh && \
    echo "source ${HOME}/conda/etc/profile.d/mamba.sh" >> /app/init.sh && \
    echo "mamba activate moflow" >> /app/init.sh && \
    chmod +x /app/init.sh && \
    echo "source /app/init.sh" >> ~/.bashrc

WORKDIR /app/moflow

RUN echo '#!/bin/bash' > /app/entrypoint.sh && \
    echo 'if [ -z "$1" ]; then' >> /app/entrypoint.sh && \
    echo '    exec /bin/bash' >> /app/entrypoint.sh && \
    echo 'else' >> /app/entrypoint.sh && \
    echo '    source /app/init.sh' >> /app/entrypoint.sh && \
    echo '    exec "$@"' >> /app/entrypoint.sh && \
    echo 'fi' >> /app/entrypoint.sh && \
    chmod +x /app/entrypoint.sh

ENTRYPOINT ["/app/entrypoint.sh"]
CMD []
HEALTHCHECK NONE
