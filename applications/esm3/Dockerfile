ARG FROM_IMAGE=ubuntu:25.10

FROM ${FROM_IMAGE} AS builder

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    git build-essential wget vim ca-certificates autoconf automake make numactl && \
    rm -rf /var/lib/apt/lists/* && \
    apt-get autoremove -y && \
    apt-get clean
ENV SERVICE_NAME="esm3-base-service"
RUN groupadd --gid 1001 $SERVICE_NAME && \
	useradd -m -g $SERVICE_NAME --shell /bin/false --uid 1001 $SERVICE_NAME

WORKDIR /app

RUN  chown -R $SERVICE_NAME:$SERVICE_NAME /app
USER $SERVICE_NAME

RUN wget --no-check-certificate "https://github.com/conda-forge/miniforge/releases/download/24.3.0-0/Miniforge3-$(uname)-$(uname -m).sh"
RUN bash Miniforge3-$(uname)-$(uname -m).sh -b -p "${HOME}/conda"

WORKDIR /app
COPY . /app
RUN git clone https://github.com/evolutionaryscale/esm.git
WORKDIR /app/esm
RUN git checkout -b esm d40007ea16850da4fbf60244a9d50c2a94cbef3d
RUN cp /app/esm3_changes.patch /app/esm/esm_changes.patch
RUN git apply /app/esm/esm_changes.patch

WORKDIR /app/esm/
RUN ${HOME}/conda/bin/conda update -y --all
RUN ${HOME}/conda/bin/mamba update -y --all
RUN ${HOME}/conda/bin/mamba env create -f /app/env.yml
ENV PIP_EXTRA_INDEX_URL=https://download.pytorch.org/whl/cpu

WORKDIR /app/esm
RUN bash -c "source ${HOME}/conda/etc/profile.d/conda.sh && \
    source ${HOME}/conda/etc/profile.d/mamba.sh && \
    mamba activate esm3 && \
    pip install --upgrade pip setuptools wheel && \
    pip install ."


RUN echo "#!/bin/bash" > /app/init.sh && \
    echo "source ${HOME}/conda/etc/profile.d/conda.sh" >> /app/init.sh && \
    echo "source ${HOME}/conda/etc/profile.d/mamba.sh" >> /app/init.sh && \
    echo "mamba activate esm3" >> /app/init.sh && \
    chmod +x /app/init.sh && \
    echo "source /app/init.sh" >> ~/.bashrc
    
WORKDIR /app/esm

RUN echo '#!/bin/bash' > /app/entrypoint.sh && \
    echo 'if [ -z "$1" ]; then' >> /app/entrypoint.sh && \
    echo '    exec /bin/bash' >> /app/entrypoint.sh && \
    echo 'else' >> /app/entrypoint.sh && \
    echo '    source /app/init.sh' >> /app/entrypoint.sh && \
    echo '    exec "$@"' >> /app/entrypoint.sh && \
    echo 'fi' >> /app/entrypoint.sh && \
    chmod +x /app/entrypoint.sh

ENTRYPOINT ["/app/entrypoint.sh"]
CMD []
HEALTHCHECK NONE

WORKDIR /app
ENV HF_HOME=/home/esm3-base-service/.cache/torch

RUN mkdir -p /home/esm3-base-service/.cache/torch && \
    ln -s /models /home/esm3-base-service/.cache/torch/hub
