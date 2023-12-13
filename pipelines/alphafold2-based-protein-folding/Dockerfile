# Install Base miniconda image
ARG BASE_IMAGE=continuumio/miniconda3
FROM ${BASE_IMAGE}

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
RUN conda install -y -c conda-forge openmm pdbfixer
RUN conda install -y -c bioconda hmmer hhsuite kalign2
RUN conda install -y jemalloc
RUN conda install -c intel onednn-cpu-iomp mkl-devel

RUN python -m pip install torch==2.0.1 pybind11
RUN python -m pip install absl-py biopython chex dm-haiku dm-tree immutabledict jax ml-collections numpy scipy tensorflow pandas psutil tqdm joblib
RUN python -m pip install jax==0.4.21 jaxlib==0.4.21
RUN python -m pip install intel-extension-for-pytorch==2.0.100

RUN git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
RUN cd /Open-Omics-Acceleration-Framework/applications/alphafold
WORKDIR /Open-Omics-Acceleration-Framework/applications/alphafold

# Compile HHsuite from source.
RUN . /opt/intel/oneapi/compiler/latest/env/vars.sh intel64 \
    && . /opt/intel/oneapi/mkl/latest/env/vars.sh intel64 \
    && git clone --recursive https://github.com/IntelLabs/hh-suite.git \
    && cd hh-suite \
    && mkdir build && cd build \
    && cmake -DCMAKE_INSTALL_PREFIX=`pwd`/release -DCMAKE_CXX_COMPILER="icpx" -DCMAKE_CXX_FLAGS_RELEASE="-O3 -march=icelake-server" .. \
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
    && autoconf && CC=icx CFLAGS="-O3 -march=icelake-server -fPIC" ./configure --prefix=`pwd`/release \
    && make -j 4 && make install \
    && ./release/bin/jackhmmer -h \
    && cd ..

RUN cd tpp-pytorch-extension \
    && CC=gcc python setup.py install \
    && python -c "from tpp_pytorch_extension.alphafold.Alpha_Attention import GatingAttentionOpti_forward" \
    && cd ..

CMD . /opt/intel/oneapi/compiler/latest/env/vars.sh intel64 \
    && . /opt/intel/oneapi/mkl/latest/env/vars.sh intel64 \
    && mkdir weights && mkdir weights/extracted && python extract_params.py --input /data/params/params_model_1.npz --output_dir ./weights/extracted/model_1 \
    && python run_multiprocess_pre.py --root_home=/Open-Omics-Acceleration-Framework/applications/alphafold --data_dir=/data --input_dir=/samples --output_dir=/output --model_name=model_1 \
    && LD_PRELOAD=/opt/conda/lib/libiomp5.so:/opt/conda/lib/libjemalloc.so:$LD_PRELOAD \
    MALLOC_CONF="oversize_threshold:1,background_thread:true,metadata_thp:auto,dirty_decay_ms:-1,muzzy_decay_ms:-1" \
    python run_multiprocess_infer.py --root_condaenv=/opt/conda --root_home=/Open-Omics-Acceleration-Framework/applications/alphafold --data_dir=/data --input_dir=/samples --output_dir=/output --model_name=model_1
    # && time bash online_preproc_baremetal.sh /Open-Omics-Acceleration-Framework/applications/alphafold /data /samples /output \
    # && time bash online_inference_baremetal.sh /opt/conda /Open-Omics-Acceleration-Framework/applications/alphafold /data /samples /output model_1
    








