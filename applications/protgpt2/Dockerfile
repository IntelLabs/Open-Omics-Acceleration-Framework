# Use the official Ubuntu image as a base
ARG FROM_IMAGE=ubuntu:24.04

# Stage 1: Set up Conda environment
ARG BASE_IMAGE=condaforge/miniforge3:24.3.0-0
FROM ${BASE_IMAGE} as conda_setup

ENV DEBIAN_FRONTEND=noninteractive

# Update Conda, install Python 3.11, and clean up
RUN conda install -c conda-forge -y python=3.11 && \
    conda clean --all -f -y

# Stage 2: Set up the main build environment
FROM ${FROM_IMAGE} as builder

ENV DEBIAN_FRONTEND=noninteractive

# Install necessary build tools and clean up
RUN apt-get update && apt-get install -y --no-install-recommends \
    git build-essential vim && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get autoremove -y && \
    apt-get clean

# Copy Conda installation from the conda_setup stage
COPY --from=conda_setup /opt/conda /opt/conda

# Copy environment.yml and create Conda environment
COPY env.yml /app/env.yml
WORKDIR /app
RUN /opt/conda/bin/conda env create -f env.yml

# Set up environment activation and PATH
RUN echo "source /opt/conda/bin/activate protgpt2_env" >> ~/.bashrc
ENV PATH="/opt/conda/envs/protgpt2_env/bin:/opt/conda/bin:$PATH"

# Create output directory
RUN mkdir -p /app/output

# Non-root user setup
ENV SERVICE_NAME="my-service"

RUN groupadd --gid 1001 $SERVICE_NAME && \
    useradd -m -g $SERVICE_NAME --shell /bin/false --uid 1001 $SERVICE_NAME && \
    chown -R $SERVICE_NAME:$SERVICE_NAME /app/output /app

# Copy application code
COPY . /app

# Switch to non-root user
USER $SERVICE_NAME


# Set the default command to run the application
#CMD ["python"]
CMD LD_PRELOAD=/opt/conda/lib/libiomp5.so:/opt/conda/lib/libjemalloc.so:$LD_PRELOAD \
    MALLOC_CONF="oversize_threshold:1,background_thread:true,metadata_thp:auto,dirty_decay_ms:-1,muzzy_decay_ms:-1" \
    python 
