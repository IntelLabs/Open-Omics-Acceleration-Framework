FROM condaforge/miniforge3:24.3.0-0

# Choose what this image should include:
# microservice | multiprocess | nextflow | all
ARG FLAVOR=microservice

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libboost-all-dev \
    swig \
    iproute2 \
    vim \
    gcc-8 \
    g++-8 \
    numactl \
    time \
    git \
    curl \
    unzip \
    ca-certificates && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
WORKDIR /opt
RUN curl -fsSL https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/releases/download/3.0/Source_code_with_submodules.tar.gz \
    -o Source_code_with_submodules.tar.gz && \
    tar -xzf Source_code_with_submodules.tar.gz && \
    mv Open-Omics-Acceleration-Framework/application/common /opt/common && \
    rm -rf Source_code_with_submodules.tar.gz Open-Omics-Acceleration-Framework

ENV PYTHONPATH="/opt/common:${PYTHONPATH}"

ENV CC=gcc-8
ENV CXX=g++-8
ENV PATH="/opt/conda/bin:${PATH}"

WORKDIR /opt
RUN git clone https://github.com/ccsb-scripps/AutoDock-Vina.git

WORKDIR /opt/AutoDock-Vina
RUN git checkout v1.2.2

WORKDIR /opt/AutoDock-Vina/build/linux/release
RUN make -j"$(nproc)"

ENV PATH="/opt/AutoDock-Vina/build/linux/release:${PATH}"

# ---- FLAVOR based installs (minimal) ----
RUN if [ "${FLAVOR}" = "microservice" ] || [ "${FLAVOR}" = "all" ]; then \
        /opt/conda/bin/python -m pip install --no-cache-dir "opea-comps" "uvicorn[standard]"; \
    fi && \
    if [ "${FLAVOR}" = "multiprocess" ] || [ "${FLAVOR}" = "all" ]; then \
        /opt/conda/bin/python -m pip install --no-cache-dir "psutil" "absl-py"; \
    fi && \
    if [ "${FLAVOR}" = "nextflow" ] || [ "${FLAVOR}" = "all" ]; then \
        curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o /tmp/awscliv2.zip && \
        unzip /tmp/awscliv2.zip -d /tmp && \
        /tmp/aws/install && \
        rm -rf /tmp/aws /tmp/awscliv2.zip; \
    fi

ENV SERVICE_NAME="autodock-vina-service"

RUN groupadd --gid 1001 ${SERVICE_NAME} && \
    useradd -m -g ${SERVICE_NAME} --shell /bin/false --uid 1001 ${SERVICE_NAME}

# Create dirs (covers all 3 use cases)
RUN mkdir -p /microservice /input /output /workspaces /work && \
    chown -R ${SERVICE_NAME}:${SERVICE_NAME} /opt /microservice /input /output /workspaces /work

USER ${SERVICE_NAME}
ENV PYTHONUNBUFFERED=1
HEALTHCHECK NONE
WORKDIR /work
CMD ["vina","--help"]

