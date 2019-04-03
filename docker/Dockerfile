# rCNV Project: Disease association analyses of rare copy number variation
# Dockerfile
# Copyright (c) 2019 Ryan L. Collins <rlcollins@g.harvard.edu> 
# Distributed under terms of the MIT License

# Base image: miniconda 3
FROM continuumio/miniconda3
MAINTAINER "Ryan Collins <rlcollins@g.harvard.edu>"

# Install make and zlib
RUN apt-get -qqy update \
    && apt-get -qqy install --fix-missing build-essential zlib1g-dev

# Create python environment
RUN conda update -n base conda
RUN conda create --name=rcnv python=3.7
RUN conda install -n rcnv -y -c bioconda numpy scipy pandas scikit-learn
RUN conda install -n rcnv -y -c bioconda pysam bzip2 pybedtools
ENV CONDA_DEFAULT_ENV rcnv
RUN echo "source activate rcnv" >> ~/.bashrc
ENV PATH /opt/conda/envs/rcnv/bin:$PATH

# Install google cloud tools
RUN pip install google-cloud-storage

# Lifted from https://github.com/GoogleCloudPlatform/cloud-sdk-docker/blob/master/debian_slim/Dockerfile
ENV CLOUD_SDK_VERSION=240.0.0
ARG INSTALL_COMPONENTS
RUN apt-get update -qqy && apt-get install -qqy \
        curl \
        gcc \
        apt-transport-https \
        lsb-release \
        openssh-client \
        gnupg \
    && pip install -U crcmod && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb https://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" > /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk=${CLOUD_SDK_VERSION}-0 $INSTALL_COMPONENTS && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image && \
    gcloud --version
VOLUME ["/root/.config"]

# Clone rCNV2 repo
RUN git clone -b master https://c9f6da1fb0dd9d5379c313700a1083f9633be849@github.com/talkowski-lab/rCNV2.git /opt/rCNV2 && \
    cd /opt/rCNV2 && \
    git checkout 3e640a97d3d89ce61b70aaa2010ab84b97d36c80 && \
    cd -

# Developer utilities: add bash aliases, etc (for convenience)
RUN echo "source ~/.bash_aliases" >> ~/.bashrc \
    && echo 'alias l="ls -ltrha"' >> ~/.bash_aliases \
    && echo 'alias less="zless"' >> ~/.bash_aliases \
    && echo 'alias median="/opt/rCNV2/utils/median.py"' >> ~/.bash_aliases

# Launch bash
CMD ["/bin/bash"]

