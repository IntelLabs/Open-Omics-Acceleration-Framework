FROM condaforge/miniforge3:4.10.2-0
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    vim \
    git \
    build-essential \
    ocl-icd-opencl-dev \
    clinfo && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
RUN conda install -c conda-forge \
    python=3.10 \
    requests=2.28.2 \
    mkl=2023.1 \
    dpcpp_linux-64=2023.1 \
    dpcpp-cpp-rt=2023.1 \
    mkl-devel=2023.1

ENV LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH}"
WORKDIR /opt
RUN git clone https://github.com/emascarenhas/AutoDock-GPU.git && \
    cd AutoDock-GPU && \
    git checkout v1.4 && \
    cd .. && \
    mv AutoDock-GPU AutoDock
WORKDIR /opt/AutoDock
RUN make DEVICE=CPU NUMWI=64
ENV PATH="/opt/AutoDock/bin:${PATH}"
WORKDIR /opt/AutoDock/bin
CMD ["/bin/bash"]

