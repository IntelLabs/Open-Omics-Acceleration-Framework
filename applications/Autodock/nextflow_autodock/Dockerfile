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
    mkl-devel=2023.1 && \
    conda clean --all -f -y
RUN pip install --no-cache-dir awscli && \
    conda clean --all -f -y

ENV LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH}"
WORKDIR /opt
ENV SERVICE_NAME="autodock-service"
RUN groupadd --gid 1001 $SERVICE_NAME && \
    useradd -m -g $SERVICE_NAME --shell /bin/false --uid 1001 $SERVICE_NAME && \
    mkdir -p /opt/AutoDock && \
    chown -R $SERVICE_NAME:$SERVICE_NAME /opt/AutoDock
USER $SERVICE_NAME
WORKDIR /opt/AutoDock
RUN git clone https://github.com/emascarenhas/AutoDock-GPU.git . && \
    git checkout v1.4
RUN make DEVICE=CPU NUMWI=64 && \
    rm -rf .git build_temp
ENV PATH="/opt/AutoDock/bin:${PATH}"
HEALTHCHECK NONE
WORKDIR /input
CMD ["autodock_cpu_64wi","--help"]

