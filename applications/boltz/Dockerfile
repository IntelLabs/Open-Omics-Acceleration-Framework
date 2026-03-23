############################################
# --- Stage 1: Conda Setup ---
############################################
ARG FROM_IMAGE=ubuntu:24.04
ARG BASE_IMAGE=condaforge/miniforge3:24.3.0-0
FROM ${BASE_IMAGE} as conda_setup

############################################
# --- Stage 2: Build environment ---
############################################
FROM ${FROM_IMAGE} as builder
ENV DEBIAN_FRONTEND=noninteractive

# 1. Install system packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    git build-essential wget vim ca-certificates numactl autoconf automake make time \
    cron procps && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get autoremove -y && \
    apt-get clean 

# 2. Host UID/GID setup
ARG USER_ID=2000
ARG GROUP_ID=2000
ENV SERVICE_NAME="boltz-service"

RUN groupadd --gid ${GROUP_ID} $SERVICE_NAME && \
    useradd -m -g $SERVICE_NAME --shell /bin/bash --uid ${USER_ID} $SERVICE_NAME

# 3. Copy Conda from previous layer
COPY --from=conda_setup /opt/conda /opt/conda
ENV PATH="/opt/conda/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH:-}"
RUN echo "source activate" >> ~/.bashrc

# 4. Install Boltz (FROM LOCAL COPY)
WORKDIR /app
RUN chown -R $SERVICE_NAME:$SERVICE_NAME /app

# COPY the local 'boltz' folder (cloned by the script)
RUN git clone --branch v2.2.1 https://github.com/jwohlwend/boltz.git

WORKDIR /app/boltz
RUN pip install -e .

# 5. Install OPEA & Server Dependencies
RUN pip install opea-comps==1.3 fastapi uvicorn anyio pydantic requests psutil absl-py

# 6. Setup Microservice Files
RUN wget https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/releases/download/4.0/Source_code_with_submodules.tar.gz
RUN tar -xzf Source_code_with_submodules.tar.gz
RUN mv /app/boltz/Open-Omics-Acceleration-Framework/applications/boltz/* /app/boltz/
RUN mv /app/boltz/Open-Omics-Acceleration-Framework/applications/common/ /app/boltz/

# 7. Create Directories and Fix Permissions
ENV BOLTZ_OUTPUT_DIR="/app/boltz_all_outputs"
ENV BOLTZ_CACHE="/app/.boltz_cache"

RUN mkdir -p ${BOLTZ_OUTPUT_DIR} && \
    mkdir -p ${BOLTZ_CACHE} && \
    chown -R ${SERVICE_NAME}:${SERVICE_NAME} /app

WORKDIR /app/boltz/

# 8. Finalize
EXPOSE 8000
USER $SERVICE_NAME

HEALTHCHECK NONE

# ENV LD_PRELOAD "/opt/conda/lib/libjemalloc.so:$LD_PRELOAD"
# ENV MALLOC_CONF "oversize_threshold:1,background_thread:true,metadata_thp:auto,dirty_decay_ms:-1,muzzy_decay_ms:-1" 

# Keep container alive forever
CMD ["sleep", "infinity"]