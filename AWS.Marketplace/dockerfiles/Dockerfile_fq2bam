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
RUN apt -y install autoconf
RUN apt -y install numactl
RUN apt -y install zlib1g-dev
RUN apt -y install libncurses5-dev
RUN apt -y install libbz2-dev
RUN apt -y install liblzma-dev
#RUN apt -y install git
RUN apt-get -qq -y install wget

#RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
#    git  --reinstall ca-certificates
    
WORKDIR /app
RUN git clone --recursive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git
WORKDIR Open-Omics-Acceleration-Framework/pipelines/fq2sortedbam
RUN git checkout e4cae79

RUN bash install.sh onprem

RUN mkdir /input
RUN mkdir /out
RUN mkdir /refdir
RUN mkdir /tempdir

CMD bash

