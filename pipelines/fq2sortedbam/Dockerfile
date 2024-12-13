FROM ubuntu:20.04
FROM python:3.10

RUN apt-get update 
RUN apt-get install -y build-essential
RUN apt-get install -y git
RUN apt -y install autoconf
RUN apt -y install numactl
RUN apt -y install zlib1g-dev
RUN apt -y install libncurses5-dev
RUN apt -y install libbz2-dev
RUN apt -y install liblzma-dev
RUN apt-get -qq -y install wget

WORKDIR /app
RUN wget https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/releases/download/3.0/Source_code_with_submodules.tar.gz 
RUN tar -xzf Source_code_with_submodules.tar.gz
RUN rm -rf Source_code_with_submodules.tar.gz
WORKDIR Open-Omics-Acceleration-Framework/pipelines/fq2sortedbam

RUN bash install.sh onprem

RUN mkdir /input
RUN mkdir /out
RUN mkdir /refdir
RUN mkdir /tempdir

CMD bash

