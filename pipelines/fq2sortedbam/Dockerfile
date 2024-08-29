ARG FROM_IMAGE=ubuntu:22.04
# PYTHON_VERSION is also set in settings.sh.
ARG PYTHON_VERSION=3.10


FROM ${FROM_IMAGE} as builder
###COPY --from=conda_setup /opt/conda /opt/conda
RUN echo "Acquire::http::proxy \"$http_proxy\";\n" \
         "Acquire::https::proxy \"$https_proxy\";" > "/etc/apt/apt.conf"

#RUN echo "http_proxy: "$http_proxy
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    git  --reinstall ca-certificates
    
WORKDIR /app
RUN git clone --recrusive https://github.com/IntelLabs/Open-Omics-Acceleration-Framework.git

WORKDIR Open-Omics-Acceleration-Framework/pipelines/fq2sortedbam

RUN bash install.sh cloud

RUN mkdir /input
RUN mkdir /out
RUN mkdir /refdir
RUN mkdir /tempdir

CMD bash

