FROM condaforge/miniforge3:4.10.2-0
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libboost-all-dev \
    swig \
    vim \
    gcc-8 \
    g++-8 \
    numactl \
    time && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
ENV CC=gcc-8
ENV CXX=g++-8
WORKDIR /opt
RUN git clone https://github.com/ccsb-scripps/AutoDock-Vina.git
WORKDIR /opt/AutoDock-Vina
RUN git checkout v1.2.2
WORKDIR /opt/AutoDock-Vina/build/linux/release
RUN make -j$(nproc)
ENV SERVICE_NAME="autodock-vina-service"
RUN groupadd --gid 1001 $SERVICE_NAME && \
    useradd -m -g $SERVICE_NAME --shell /bin/false --uid 1001 $SERVICE_NAME
RUN chown -R $SERVICE_NAME:$SERVICE_NAME /opt
USER $SERVICE_NAME
ENV PATH="/opt/AutoDock-Vina/build/linux/release:$PATH"
WORKDIR /input
HEALTHCHECK NONE
CMD ["vina","--help"]

