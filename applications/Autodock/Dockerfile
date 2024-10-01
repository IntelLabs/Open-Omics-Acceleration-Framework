FROM condaforge/miniforge3:4.10.2-0

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    vim \
    git \
    build-essential \
    ocl-icd-opencl-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV SERVICE_NAME="autodock-service"
RUN groupadd --gid 1001 $SERVICE_NAME && \
    useradd -m -g $SERVICE_NAME --shell /bin/false --uid 1001 $SERVICE_NAME

WORKDIR /opt
RUN chown -R $SERVICE_NAME:$SERVICE_NAME /opt

USER $SERVICE_NAME

RUN conda install -c conda-forge \
    python=3.11 \
    mkl=2023.1 \
    dpcpp_linux-64=2023.1 \
    dpcpp-cpp-rt=2023.1 \
    mkl-devel=2023.1

WORKDIR /opt

RUN git clone https://github.com/emascarenhas/AutoDock-GPU.git && \
    cd AutoDock-GPU && \
    git checkout v1.4 && \
    mv ../AutoDock-GPU ../AutoDock

WORKDIR /opt/AutoDock

RUN make DEVICE=CPU NUMWI=64

ENV PATH="/opt/AutoDock/bin:${PATH}"

HEALTHCHECK NONE

CMD ["/bin/bash"]


