FROM ubuntu:20.04

SHELL ["/bin/bash", "-c"]
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Asia/Kolkata
ENV SERVICE_NAME="gromacs-service"

RUN set -euo pipefail && \
    apt-get update && apt-get install -y --no-install-recommends \
    wget \
    gnupg \
    software-properties-common \
    time && \
    wget -qO- https://apt.kitware.com/keys/kitware-archive-latest.asc | gpg --dearmor > /usr/share/keyrings/kitware-archive-keyring.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/kitware-archive-keyring.gpg] https://apt.kitware.com/ubuntu/ focal main" > /etc/apt/sources.list.d/kitware.list && \
    apt-get update && apt-get install -y --no-install-recommends \
    cmake \
    git \
    tar \
    vim \
    g++ \
    gcc \
    make \
    curl \
    build-essential \
    tzdata && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O Miniforge3.sh && \
    bash Miniforge3.sh -b -p /opt/conda && \
    rm -f Miniforge3.sh

RUN chown -R 1001:1001 /opt/conda
RUN groupadd --gid 1001 $SERVICE_NAME && \
    useradd -m -g $SERVICE_NAME --shell /bin/bash --uid 1001 $SERVICE_NAME && \
    mkdir -p /grmcs /input /output && \
    chown -R $SERVICE_NAME:$SERVICE_NAME /grmcs /input /output

ENV PATH="/opt/conda/bin:$PATH"

RUN conda init bash && \
    conda create -y -n grms_env python=3.12 \
    sphinx pygments \
    mkl=2023.2 mkl-devel=2023.2 \
    dpcpp_linux-64=2023.2 dpcpp-cpp-rt=2023.2 \
    gxx_linux-64=12 && \
    conda clean -afy && \
    echo "source /opt/conda/bin/activate grms_env" >> /etc/bash.bashrc && \
    echo "conda activate grms_env" >> /etc/bash.bashrc
RUN source /opt/conda/bin/activate grms_env && \
    pip install --no-cache-dir awscli && \
    conda clean -afy
RUN source /opt/conda/bin/activate grms_env && conda activate grms_env

ENV CC=icx \
    CXX=icpx \
    MKLROOT=/opt/conda/envs/grms_env \
    CFLAGS="-fopenmp -I/opt/conda/envs/grms_env/include" \
    CXXFLAGS="-fopenmp -stdlib=libstdc++ -I/opt/conda/envs/grms_env/include" \
    LDFLAGS="-L/opt/conda/envs/grms_env/lib -lmkl_rt -liomp5 -lpthread -lm -ldl" \
    LD_LIBRARY_PATH="/opt/conda/envs/grms_env/lib:$LD_LIBRARY_PATH"

WORKDIR /grmcs

RUN source /opt/conda/bin/activate grms_env && \
    wget -q https://ftp.gromacs.org/gromacs/gromacs-2024.2.tar.gz && \
    tar -xzvf gromacs-2024.2.tar.gz && \
    rm -f gromacs-2024.2.tar.gz

RUN chown -R $SERVICE_NAME:$SERVICE_NAME /grmcs/gromacs-2024.2
WORKDIR /grmcs/gromacs-2024.2
RUN mkdir -p build && chmod -R 775 build

WORKDIR /grmcs/gromacs-2024.2/build

RUN source /opt/conda/bin/activate grms_env && \
    which icx && icx --version && \
    which icpx && icpx --version

RUN source /opt/conda/bin/activate grms_env && \
    cmake .. -DGMX_BUILD_OWN_FFTW=OFF \
             -DREGRESSIONTEST_DOWNLOAD=OFF \
             -DCMAKE_C_COMPILER=$CC \
             -DCMAKE_CXX_COMPILER=$CXX \
             -DGMX_FFT_LIBRARY=mkl \
             -DMKL_INCLUDE_DIR=$MKLROOT/include \
             -DMKL_LIBRARIES=$MKLROOT/lib/libmkl_rt.so \
             -DCMAKE_INSTALL_PREFIX=/grmcs/gmx_mkl_prefix && \
    make -j10 VERBOSE=1 && \
    make install

WORKDIR /grmcs
RUN rm -rf gromacs-2024.2
ENV PATH="/grmcs/gmx_mkl_prefix/bin:$PATH"
RUN echo "source /grmcs/gmx_mkl_prefix/bin/GMXRC" >> /etc/bash.bashrc

WORKDIR /input
COPY run_commands.sh /input/
RUN chmod +x /input/run_commands.sh

COPY entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh

USER $SERVICE_NAME
WORKDIR /input
ENTRYPOINT ["/entrypoint.sh"]
HEALTHCHECK NONE

