FROM condaforge/miniforge3
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    libboost-all-dev \
    swig \
    gcc-8 \
    g++-8 \
    numactl \
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
RUN git clone https://github.com/ccsb-scripps/AutoDock-Vina.git && \
    cd AutoDock-Vina && \
    git checkout v1.2.2 && \
    git switch -c v1.2.2

ENV CC=gcc-8
ENV CXX=g++-8
RUN cd AutoDock-Vina/build/linux/release && \
    make -j$(nproc)

ENV PATH="/root/miniforge3/envs/vina_1.2.2_env/bin:/AutoDock-Vina/build/linux/release:$PATH"

WORKDIR /AutoDock-Vina/build/linux/release
CMD ["/bin/bash"]
