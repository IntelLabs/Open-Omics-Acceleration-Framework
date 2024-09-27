FROM condaforge/miniforge3:4.10.2-0
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install -y \
    vim \
    git \
    build-essential \
    ocl-icd-opencl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
RUN conda install -c conda-forge \
    python=3.11 \
    mkl=2023.1 \
    dpcpp_linux-64=2023.1 \
    dpcpp-cpp-rt=2023.1 \
    mkl-devel=2023.1
RUN git clone https://github.com/emascarenhas/AutoDock-GPU.git && \
    cd AutoDock-GPU && \
    git checkout v1.4 && \
    make DEVICE=CPU NUMWI=16
CMD ["/bin/bash"]

