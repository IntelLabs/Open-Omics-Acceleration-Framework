FROM condaforge/miniforge3

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libboost-all-dev \
    swig \
    gcc-8 \
    g++-8 \
    numactl \
    time \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN conda create -n vina_1.2.2_env python=3.11 -y && \
    echo "source activate vina_1.2.2_env" >> ~/.bashrc

RUN conda install -n vina_1.2.2_env -c conda-forge \
    numpy \
    swig \
    boost-cpp \
    sphinx \
    sphinx_rtd_theme \
    -y

WORKDIR /opt
RUN git clone https://github.com/ccsb-scripps/AutoDock-Vina.git && \
    cd AutoDock-Vina && \
    git checkout v1.2.2 && \
    git switch -c v1.2.2

ENV CC=gcc-8
ENV CXX=g++-8

WORKDIR /opt/AutoDock-Vina/build/linux/release
RUN make -j$(nproc)

ENV PATH="/opt/miniforge3/envs/vina_1.2.2_env/bin:/opt/AutoDock-Vina/build/linux/release:$PATH"

ENV SERVICE_NAME="autodock-vina-service"

RUN groupadd --gid 1001 $SERVICE_NAME && \
    useradd -m -g $SERVICE_NAME --shell /bin/false --uid 1001 $SERVICE_NAME

RUN chown -R $SERVICE_NAME:$SERVICE_NAME /opt

USER $SERVICE_NAME

HEALTHCHECK NONE

WORKDIR /opt/AutoDock-Vina/build/linux/release

CMD ["/bin/bash"]

