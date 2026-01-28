FROM condaforge/miniforge3:24.7.1-0

ENV DEBIAN_FRONTEND=noninteractive
ENV SERVICE_NAME="autodock-service"

RUN apt-get update && apt-get install -y --no-install-recommends \
    iproute2 \
    vim \
    git \
    build-essential \
    ocl-icd-opencl-dev \
    clinfo \
    unzip \
    wget \
    tar \
    ca-certificates \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Keep conda base healthy: do NOT pin requests here
RUN conda update -n base -c conda-forge --all -y && \
    conda install -n base -c conda-forge \
      python=3.10 \
      mkl=2023.1 \
      dpcpp_linux-64=2023.1 \
      dpcpp-cpp-rt=2023.1 \
      mkl-devel=2023.1 \
    && conda clean --all -f -y

WORKDIR /opt
ENV LD_LIBRARY_PATH="/opt/conda/lib:${LD_LIBRARY_PATH}"

# AWS CLI v2 (zip install) - avoids pip awscli v1 issues
RUN wget -q "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -O /tmp/awscliv2.zip && \
    unzip -q /tmp/awscliv2.zip -d /tmp && \
    /tmp/aws/install && \
    rm -rf /tmp/aws /tmp/awscliv2.zip

RUN wget -q https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/releases/download/3.0/Source_code_with_submodules.tar.gz && \
    tar -xzf Source_code_with_submodules.tar.gz && \
    mv Open-Omics-Acceleration-Framework/application/common /opt/common && \
    rm -rf Source_code_with_submodules.tar.gz Open-Omics-Acceleration-Framework

ENV PYTHONPATH="/opt/common:${PYTHONPATH}"


RUN python -m pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir opea-comps uvicorn[standard]

RUN groupadd --gid 1001 ${SERVICE_NAME} && \
    useradd -m -g ${SERVICE_NAME} --shell /bin/false --uid 1001 ${SERVICE_NAME}

WORKDIR /opt/AutoDock-src
RUN git clone https://github.com/emascarenhas/AutoDock-GPU.git . && \
    git checkout v1.4 && \
    make DEVICE=CPU NUMWI=64 && \
    rm -rf .git build_temp

RUN mkdir -p /opt/AutoDock/bin && \
    cp -v /opt/AutoDock-src/bin/autodock_cpu_64wi /opt/AutoDock/bin/ && \
    chmod +x /opt/AutoDock/bin/autodock_cpu_64wi && \
    rm -rf /opt/AutoDock-src

ENV PATH="/opt/AutoDock/bin:${PATH}"

RUN mkdir -p /input /output && chmod 777 /input /output
RUN mkdir -p /opt/AutoDock/microservice
USER ${SERVICE_NAME}
WORKDIR /input
HEALTHCHECK NONE
CMD ["autodock_cpu_64wi", "--help"]

