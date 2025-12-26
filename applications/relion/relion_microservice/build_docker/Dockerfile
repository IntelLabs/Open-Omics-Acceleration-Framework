FROM ubuntu:22.04

LABEL description="RELION 5.0 built with Intel oneAPI and conda"

ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Kolkata
ENV SERVICE_NAME="relion-service"

# Create non-root user
RUN groupadd --gid 1001 $SERVICE_NAME && \
    useradd -m -g $SERVICE_NAME --shell /bin/bash --uid 1001 $SERVICE_NAME

# Install system packages and Python build dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    iproute2 \
    git \
    cmake \
    vim \
    make \
    tar \
    wget \
    unzip \
    gnupg \
    time \
    ca-certificates \
    tzdata \
    libtiff-dev \
    libx11-dev \
    libpng-dev \
    python3-dev \
    libffi-dev \
    libssl-dev \
    pkg-config \
    gfortran \
    libstdc++-11-dev \
    libgl1 \
    libgl1-mesa-glx \
    libxrender1 \
    build-essential && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN wget "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -O "/tmp/awscliv2.zip" && \
    unzip /tmp/awscliv2.zip -d /tmp && \
    /tmp/aws/install && \
    rm -rf /tmp/aws /tmp/awscliv2.zip

WORKDIR /opt
RUN wget -q https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/releases/download/3.0/Source_code_with_submodules.tar.gz && \
    tar -xzf Source_code_with_submodules.tar.gz && \
    mv Open-Omics-Acceleration-Framework/application/common /opt/common && \
    rm -rf Source_code_with_submodules.tar.gz Open-Omics-Acceleration-Framework

ENV PYTHONPATH="/opt/common:${PYTHONPATH}"
RUN mkdir -p /microservice && chown $SERVICE_NAME:$SERVICE_NAME /microservice
RUN chown -R $SERVICE_NAME:$SERVICE_NAME /opt
USER $SERVICE_NAME

# Clone and patch RELION
RUN git clone --branch ver5.0 https://github.com/3dem/relion.git relion_5.0
WORKDIR /opt/relion_5.0
RUN git pull
COPY relion_env_patch.patch /opt/relion_5.0
RUN git apply relion_env_patch.patch

# Install Intel oneAPI HPC Toolkit
WORKDIR /opt
RUN wget https://registrationcenter-download.intel.com/akdlm/IRC_NAS/b7f71cf2-8157-4393-abae-8cea815509f7/intel-oneapi-hpc-toolkit-2025.0.1.47_offline.sh && \
    chmod +x intel-oneapi-hpc-toolkit-2025.0.1.47_offline.sh && \
    ./intel-oneapi-hpc-toolkit-2025.0.1.47_offline.sh -a --silent --cli --eula accept && \
    rm intel-oneapi-hpc-toolkit-2025.0.1.47_offline.sh

# Conda setup
ENV CONDA_DIR=/opt/conda
ENV HOME=/home/relion-service
ENV XDG_CACHE_HOME=$HOME/.cache
ENV TMPDIR=/home/relion-service/tmp
ENV PATH=$CONDA_DIR/bin:$PATH
ENV PIP_EXTRA_INDEX_URL=https://download.pytorch.org/whl/cpu

RUN mkdir -p $XDG_CACHE_HOME $TMPDIR && \
    chmod -R 777 $XDG_CACHE_HOME $TMPDIR 

# Install Miniforge
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && \
    bash /tmp/miniforge.sh -b -u -p $CONDA_DIR && \
    rm /tmp/miniforge.sh && \
    conda init bash

# Create conda env without .[vis] first
RUN conda env create -f /opt/relion_5.0/environment.yml && \
    conda clean -afy

# Install .[vis] manually after env is created
RUN bash -c "source $CONDA_DIR/etc/profile.d/conda.sh && \
    conda activate relion-5.0 && \
    pip install --verbose /opt/relion_5.0[vis]"

RUN bash -c "source $CONDA_DIR/etc/profile.d/conda.sh && \
    conda activate relion-5.0 && \
    pip install --no-cache-dir \
        'opea-comps' \
        'uvicorn[standard]'"

# Build RELION
WORKDIR /opt/relion_5.0
SHELL ["/bin/bash", "-c"]
RUN mkdir -p /opt/relion_5.0/torch_home
RUN mkdir -p /opt/relion_5.0/build_cpu
WORKDIR /opt/relion_5.0/build_cpu

RUN source /home/relion-service/intel/oneapi/2025.0/oneapi-vars.sh --force && \
    cmake -DCMAKE_C_COMPILER=icx \
          -DCMAKE_CXX_COMPILER=icpx \
          -DMPI_C_COMPILER=mpiicx \
          -DMPI_CXX_COMPILER=mpiicpx \
          -DCUDA=OFF \
          -DALTCPU=ON \
          -DMKLFFT=ON \
          -DGUI=OFF \
          -DFETCH_WEIGHTS=OFF \
          -DCMAKE_BUILD_TYPE=Release \
          -DCMAKE_C_FLAGS="-g -O3 -qopenmp-simd -xCORE-AVX512 -qopt-zmm-usage=high" \
          -DCMAKE_CXX_FLAGS="-DTBB_SUPPRESS_DEPRECATED_MESSAGES -g -O3 -qopenmp-simd -xCORE-AVX512 -qopt-zmm-usage=high" \
          -DCMAKE_EXE_LINKER_FLAGS="-static-intel -static-libgcc -qopenmp-link=static -Wno-unused-command-line-argument" \
          -DTORCH_HOME_PATH=/opt/relion_5.0/torch_home \
          -DPYTHON_EXE_PATH=/opt/conda/envs/relion-5.0/bin/python \
          -DCMAKE_INSTALL_PREFIX=/opt/relion_5.0_cpu_benchmark_prefix .. && \
    make -j"$(nproc)" && \
    make install

WORKDIR /opt/relion_5.0
RUN rm -rf /opt/relion_5.0/build_cpu
RUN echo 'source /home/relion-service/intel/oneapi/setvars.sh --force' >> /home/relion-service/.bashrc
ENV PATH="/opt/relion_5.0_cpu_benchmark_prefix/bin:/home/relion-service/intel/oneapi/compiler/2025.0/bin:/home/relion-service/intel/oneapi/mpi/2025.0/bin:/home/relion-service/intel/oneapi/compiler/2025.0/linux/bin/intel64:$PATH"
ENV LD_LIBRARY_PATH="/home/relion-service/intel/oneapi/compiler/2025.0/lib:/home/relion-service/intel/oneapi/mkl/2025.0/lib/intel64:/home/relion-service/intel/oneapi/2025.0/lib:$LD_LIBRARY_PATH"
RUN mkdir -p /opt/relion_5.0/relion_benchmark
COPY entrypoint.sh ./entrypoint_temp.sh
RUN bash -c "cp ./entrypoint_temp.sh ./entrypoint.sh && chmod u+x ./entrypoint.sh && rm ./entrypoint_temp.sh"

ENTRYPOINT ["/opt/relion_5.0/entrypoint.sh"]

HEALTHCHECK NONE

