#ARG FROM_IMAGE=ubuntu:20.04
# PYTHON_VERSION is also set in settings.sh.
#ARG PYTHON_VERSION=3.10

#FROM ${FROM_IMAGE} as builder
FROM ubuntu:20.04
FROM python:3.10

#RUN echo "Acquire::http::proxy \"$http_proxy\";\n" \
#         "Acquire::https::proxy \"$https_proxy\";" > "/etc/apt/apt.conf"
RUN apt-get update 
##&& DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
RUN apt-get install -y build-essential
RUN apt-get install -y git
RUN apt-get install -y vim
RUN apt -y install autoconf
RUN apt -y install numactl
RUN apt -y install zlib1g-dev
RUN apt -y install libncurses5-dev
RUN apt -y install libbz2-dev
RUN apt -y install liblzma-dev
RUN apt-get -qq -y install wget gcc-11 g++-11
RUN unlink /usr/bin/gcc && unlink /usr/bin/g++
RUN ln -s /usr/bin/gcc-11 /usr/bin/gcc && ln -s /usr/bin/g++-11 /usr/bin/g++

#RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
#    git  --reinstall ca-certificates
    
WORKDIR /app
RUN git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git -b dev2
#COPY /scratch/users/mvasimud/workspace/Open-Omics-Acceleration-Framework .
WORKDIR Open-Omics-Acceleration-Framework/pipelines/fq2sortedbam

RUN bash install.sh onprem
ENV PATH="/app/Open-Omics-Acceleration-Framework/pipelines/fq2sortedbam/miniforge3/envs/fq2bam/bin/:/app/Open-Omics-Acceleration-Framework/pipelines/fq2sortedbam/miniforge3/bin/:$PATH"
ENV LD_PRELOAD="/app/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/libmimalloc.so.2.0"

RUN mkdir /input
RUN mkdir /out
RUN mkdir /refdir
RUN mkdir /tempdir

CMD ["/bin/bash"]
SHELL ["source", "~/miniforge3/bin/activate", "fq2bam", "/bin/bash", "-c"]
