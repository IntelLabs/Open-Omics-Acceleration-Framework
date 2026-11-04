# fq2sortedbam aligner image. Uses the local checkout (COPY) instead of
# re-cloning from GitHub, and uses the default gcc instead of gcc-11 (the
# `python:3.10` tag now resolves to a Debian release that no longer ships
# gcc-11).
#
# Build (from the repository root):
#   docker build -f pipelines/fq2sortedbam/Dockerfile -t fq2sortedbam:latest .

FROM ubuntu:20.04
FROM python:3.10

RUN apt-get update
RUN apt-get install -y build-essential
RUN apt-get install -y git
RUN apt-get install -y vim
RUN apt -y install autoconf
RUN apt -y install numactl
RUN apt -y install zlib1g-dev
RUN apt -y install libncurses5-dev
RUN apt -y install libbz2-dev
RUN apt -y install liblzma-dev
RUN apt-get -qq -y install wget

WORKDIR /app/Open-Omics-Acceleration-Framework
COPY . .
# Vendored bwa-mem2 submodule (ext/safestringlib) is missing <stdlib.h>/<ctype.h>
# includes needed for abort()/toupper() under modern GCC; patch before building.
RUN sed -i '/#include <stdio.h>/a #include <stdlib.h>\n#include <ctype.h>' \
    applications/bwa-mem2/ext/safestringlib/safeclib/safeclib_private.h
WORKDIR /app/Open-Omics-Acceleration-Framework/pipelines/fq2sortedbam

RUN bash install.sh onprem
ENV PATH="/app/Open-Omics-Acceleration-Framework/applications/samtools:$PATH"
ENV PATH="/app/Open-Omics-Acceleration-Framework/applications/bwa-mem2:$PATH"
ENV PATH="/app/Open-Omics-Acceleration-Framework/applications/bwa-meth:$PATH"
ENV PATH="/app/Open-Omics-Acceleration-Framework/applications/mm2-fast:$PATH"
ENV PATH="/app/Open-Omics-Acceleration-Framework/applications/STAR/source:$PATH"
ENV PATH="/app/Open-Omics-Acceleration-Framework/applications/samtools/samtools:$PATH"

ENV PATH="/app/Open-Omics-Acceleration-Framework/pipelines/fq2sortedbam/miniforge3/envs/fq2bam/bin/:/app/Open-Omics-Acceleration-Framework/pipelines/fq2sortedbam/miniforge3/bin/:$PATH"
ENV LD_PRELOAD="/app/Open-Omics-Acceleration-Framework/pipelines/deepvariant-based-germline-variant-calling-fq2vcf/libmimalloc.so.2.0"

RUN mkdir /input
RUN mkdir /out
RUN mkdir /refdir
RUN mkdir /tempdir

CMD ["/bin/bash"]
