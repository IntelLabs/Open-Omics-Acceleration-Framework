# MIT License

# Copyright (c) 2022 Intel Labs

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# Authors: Narendra Chaudhary <narendra.chaudhary@intel.com>; Sanchit Misra <sanchit.misra@intel.com>

# Install Base miniconda image
ARG FROM_IMAGE=ubuntu:22.04

# Install Base miniforge image
ARG BASE_IMAGE=condaforge/miniforge3:23.1.0-3
FROM ${BASE_IMAGE} as conda_setup

RUN conda update -n base conda
COPY ./environment.yml ./

RUN conda env create --name=single_cell -f ./environment.yml
COPY ./_t_sne.py /opt/conda/lib/python3.8/site-packages/daal4py/sklearn/manifold/_t_sne.py


FROM ${FROM_IMAGE} as builder
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
    git build-essential gcc curl gnupg gnupg2 gnupg1 sudo wget tar ca-certificates -y \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get autoremove -y \
    && apt-get clean \
    && apt update


COPY --from=conda_setup /opt/conda /opt/conda
ENV PATH "/opt/conda/envs/single_cell/bin:$PATH"
RUN echo "source /opt/conda/bin/activate single_cell" >> ~/.bashrc
CMD source ~/.bashrc

# Non-root user setup
ENV SERVICE_NAME="scrna"

RUN groupadd --gid 1001 $SERVICE_NAME && \
    useradd -m -g $SERVICE_NAME --shell /bin/false --uid 1001 $SERVICE_NAME


RUN pip uninstall -y umap-learn
WORKDIR /
RUN wget https://github.com/IntelLabs/Open-Omics-Acceleration-Framework/releases/download/3.0/Source_code_with_submodules.tar.gz 
RUN tar -xzf Source_code_with_submodules.tar.gz

#SHELL ["/bin/bash", "-c", "source activate single_cell"]
WORKDIR /Open-Omics-Acceleration-Framework/lib/tal/applications/UMAP_fast/umap_extend/
RUN python setup.py install

WORKDIR ../umap/
RUN python setup.py install

RUN chown -R $SERVICE_NAME:$SERVICE_NAME /Open-Omics-Acceleration-Framework /opt
# Switch to non-root user
USER $SERVICE_NAME


WORKDIR /Open-Omics-Acceleration-Framework/pipelines/single-cell-RNA-seq-analysis/notebooks


CMD jupyter notebook \
                --no-browser \
                --allow-root \
                --port=8888 \
                --ip=0.0.0.0 \
                --notebook-dir=/Open-Omics-Acceleration-Framework/pipelines/single-cell-RNA-seq-analysis/notebooks \
                --NotebookApp.password="" \
                --NotebookApp.token="" \
                --NotebookApp.password_required=False



