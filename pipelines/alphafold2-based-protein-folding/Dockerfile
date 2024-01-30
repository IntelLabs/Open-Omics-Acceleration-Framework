ARG FROM_IMAGE=ubuntu:20.04
# Install Base miniconda image
ARG BASE_IMAGE=continuumio/miniconda3
FROM ${BASE_IMAGE} as conda_setup

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    git build-essential cmake wget tzdata gcc curl gnupg gnupg2 gnupg1 sudo kalign autoconf numactl time vim \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get autoremove -y \
    && apt-get clean

RUN wget -O- https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB \
| gpg --dearmor | tee /usr/share/keyrings/oneapi-archive-keyring.gpg > /dev/null
RUN echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list
RUN apt-get update && \
  DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    curl ca-certificates gpg-agent software-properties-common intel-basekit intel-hpckit && \
  rm -rf /var/lib/apt/lists/*

# Install Anaconda and PIP dependency
RUN conda update -n base conda
RUN conda install python==3.9
RUN conda install -y -c conda-forge openmm=8.0.0=py39h31b9435_1 pdbfixer=1.9=pyh1a96a4e_0
RUN conda install -y -c bioconda hmmer=3.3.2=h87f3376_2 hhsuite=3.3.0=py39pl5321h67e14b5_5 kalign2=2.04=hec16e2b_3
RUN conda install -y jemalloc=5.2.1=h6a678d5_6
RUN conda install -c intel onednn-cpu-iomp=2023.2.0=intel_49516 mkl-devel=2023.2.0=intel_49495
RUN python -m pip install torch==2.0.1 pybind11==2.11.1
RUN python -m pip install absl-py==2.0.0 biopython==1.81 chex==0.1.84 dm-haiku==0.0.10 dm-tree==0.1.8 immutabledict==3.0.0  ml-collections==0.1.1 numpy==1.26.1 scipy==1.11.3 tensorflow==2.14.0 pandas==2.1.1 psutil==5.9.6 tqdm==4.65.0 joblib==1.3.2 pragzip==0.6.0
RUN python -m pip install jax==0.4.21 jaxlib==0.4.21
RUN python -m pip install intel-extension-for-pytorch==2.0.100
#COPY . /Open-Omics-Acceleration-Framework/
RUN git clone https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
RUN cd /Open-Omics-Acceleration-Framework/applications/alphafold
WORKDIR /Open-Omics-Acceleration-Framework/applications/alphafold
RUN git submodule update --recursive --init -- .
# Compile HHsuite from source.
RUN . /opt/intel/oneapi/compiler/latest/env/vars.sh intel64 \
    && . /opt/intel/oneapi/mkl/latest/env/vars.sh intel64 \
    && git clone --recursive https://github.com/IntelLabs/hh-suite.git \
    && cd hh-suite \
    && mkdir build && cd build \
    && cmake -DCMAKE_INSTALL_PREFIX=`pwd`/release -DCMAKE_CXX_COMPILER="icpx" -DCMAKE_CXX_FLAGS_RELEASE="-O3 -march=native" .. \
    && make -j 4 && make install \
    && ./release/bin/hhblits -h \
    && cd .. && cd .. \

# Compile Hmmer from source.
RUN . /opt/intel/oneapi/compiler/latest/env/vars.sh intel64 \
    && . /opt/intel/oneapi/mkl/latest/env/vars.sh intel64 \
    && git clone --recursive https://github.com/IntelLabs/hmmer.git \
    && cd hmmer \
    && cp easel_makefile.in easel/Makefile.in \
    && cd easel && autoconf && ./configure --prefix=`pwd` && cd .. \
    && autoconf && CC=icx CFLAGS="-O3 -march=native -fPIC" ./configure --prefix=`pwd`/release \
    && make -j 4 && make install \
    && ./release/bin/jackhmmer -h \
    && cd ..

RUN cd tpp-pytorch-extension \
    && CC=gcc python setup.py install \
    && python -c "from tpp_pytorch_extension.alphafold.Alpha_Attention import GatingAttentionOpti_forward" \
    && cd ..


FROM ${FROM_IMAGE} as builder
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    git build-essential cmake wget tzdata gcc curl gnupg gnupg2 gnupg1 sudo kalign autoconf numactl time vim \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get autoremove -y \
    && apt-get clean
COPY --from=conda_setup /opt/conda /opt/conda
COPY  --from=conda_setup /Open-Omics-Acceleration-Framework /Open-Omics-Acceleration-Framework
WORKDIR /Open-Omics-Acceleration-Framework/applications/alphafold

RUN echo "source /opt/conda/bin/activate " > ~/.bashrc
CMD source ~/.bashrc
ENV PATH "/opt/conda/bin:$PATH"
ENV LD_LIBRARY_PATH "/opt/conda/lib:$LD_LIBRARY_PATH"
CMD mkdir weights && mkdir weights/extracted && python extract_params.py --input /data/params/params_model_1.npz --output_dir ./weights/extracted/model_1 \
    && python run_multiprocess_pre.py --root_home=/Open-Omics-Acceleration-Framework/applications/alphafold --data_dir=/data --input_dir=/samples --output_dir=/output --model_name=model_1 \
    && LD_PRELOAD=/opt/conda/lib/libiomp5.so:/opt/conda/lib/libjemalloc.so:$LD_PRELOAD \
    MALLOC_CONF="oversize_threshold:1,background_thread:true,metadata_thp:auto,dirty_decay_ms:-1,muzzy_decay_ms:-1" \
    python run_multiprocess_infer.py --root_condaenv=/opt/conda --root_home=/Open-Omics-Acceleration-Framework/applications/alphafold --data_dir=/data --input_dir=/samples --output_dir=/output --model_names=model_1
