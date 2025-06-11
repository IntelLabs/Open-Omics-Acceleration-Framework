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

# Build arguments for host UID/GID
ARG USER_ID=2000
ARG GROUP_ID=2000

ENV SERVICE_NAME="boltz-service"

# Create a user and group with same UID and GID as host
RUN groupadd --gid ${GROUP_ID} $SERVICE_NAME && \
    useradd -m -g $SERVICE_NAME --shell /bin/false --uid ${USER_ID} $SERVICE_NAME

# Copy Conda installation from the conda_setup stage
COPY --from=conda_setup /opt/conda /opt/conda
ENV PATH="/opt/conda/bin:$PATH"
ENV LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH:-}"
RUN echo "source activate" >> ~/.bashrc

RUN git clone --branch 5.3.0 https://github.com/jemalloc/jemalloc.git
WORKDIR /jemalloc
RUN bash autogen.sh --prefix=/opt/conda/ && make install
WORKDIR /
RUN rm -rf jemalloc

WORKDIR /app
RUN chown -R $SERVICE_NAME:$SERVICE_NAME /app
RUN git clone --branch v0.4.1 https://github.com/jwohlwend/boltz.git

WORKDIR /app/boltz
RUN pip install -e .

# Switch to non-root user
RUN chown -R $SERVICE_NAME:$SERVICE_NAME /app/boltz
COPY ./entrypoint.sh entrypoint.sh
RUN chmod +x entrypoint.sh

USER $SERVICE_NAME

# bin bash
# Clone Boltz 1 repository (replace with the actual repo URL)
ENTRYPOINT ["/app/boltz/entrypoint.sh"]

# Default command
CMD ["default"]