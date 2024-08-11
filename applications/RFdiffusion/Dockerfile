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
    git build-essential wget vim ca-certificates && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get autoremove -y && \
    apt-get clean

# Copy Conda installation from the conda_setup stage
COPY --from=conda_setup /opt/conda /opt/conda

COPY . /app

WORKDIR /app
RUN git clone https://github.com/RosettaCommons/RFdiffusion.git
WORKDIR /app/RFdiffusion

RUN mkdir -p /app/RFdiffusion/models 
#    bash /app/download_models.sh /app/RFdiffusion/models

# Apply the patch file
COPY new_changes.patch /app/RFdiffusion/new_changes.patch
RUN git apply new_changes.patch

RUN /opt/conda/bin/conda env create -f env/SE3nv.yml
# Set up environment activation and PATH
RUN echo "source /opt/conda/bin/activate SE3nv" >> ~/.bashrc
ENV PATH="/opt/conda/envs/SE3nv/bin:/opt/conda/bin:$PATH"
	
# Install dependencies
WORKDIR /app/RFdiffusion/env/SE3Transformer
RUN pip install --no-cache-dir -r requirements.txt
RUN python setup.py install

WORKDIR /app/RFdiffusion
RUN pip install -e .

WORKDIR /app/RFdiffusion
RUN tar -xvf examples/ppi_scaffolds_subset.tar.gz -C examples/

# Set the default command to run the application
#CMD ["python"]
CMD LD_PRELOAD=/opt/conda/lib/libiomp5.so:/opt/conda/lib/libjemalloc.so:$LD_PRELOAD \
    MALLOC_CONF="oversize_threshold:1,background_thread:true,metadata_thp:auto,dirty_decay_ms:-1,muzzy_decay_ms:-1" \
    python \ scripts/run_inference.py
