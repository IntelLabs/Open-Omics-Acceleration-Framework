FROM condaforge/miniforge3:4.10.2-0

ENV DEBIAN_FRONTEND=noninteractive

# Install necessary packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    vim \
    git \
    build-essential \
    ocl-icd-opencl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Non-root user setup
ENV SERVICE_NAME="autodock-service"
RUN groupadd --gid 1001 $SERVICE_NAME && \
    useradd -m -g $SERVICE_NAME --shell /bin/false --uid 1001 $SERVICE_NAME

WORKDIR /opt
RUN chown -R $SERVICE_NAME:$SERVICE_NAME /opt

USER $SERVICE_NAME

# Install dependencies via conda
RUN conda install -c conda-forge \
    python=3.11 \
    mkl=2023.1 \
    dpcpp_linux-64=2023.1 \
    dpcpp-cpp-rt=2023.1 \
    mkl-devel=2023.1

# Set the working directory for AutoDock-GPU
WORKDIR /opt

# Clone the AutoDock-GPU repository and build it
RUN git clone https://github.com/emascarenhas/AutoDock-GPU.git && \
    cd AutoDock-GPU && \
    git checkout v1.4 && \
    mv ../AutoDock-GPU ../AutoDock

WORKDIR /opt/AutoDock

# Build AutoDock-GPU
RUN make DEVICE=CPU NUMWI=64

# Update the PATH environment variable to include AutoDock binaries
ENV PATH="/opt/AutoDock/bin:${PATH}"

# Change ownership to the non-root user
#RUN chown -R $SERVICE_NAME:$SERVICE_NAME /opt

# Switch to the non-root user
#USER $SERVICE_NAME

# Disable health checks
HEALTHCHECK NONE

# Set the default command for the container
CMD ["/bin/bash"]


