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
RUN conda config --set ssl_verify no
RUN conda install -n rcnv -y r-base
RUN conda install -n rcnv -y -c bioconda numpy scipy pandas scikit-learn
RUN conda install -n rcnv -y -c bioconda samtools=1.9 pysam=0.15.3
RUN conda install -n rcnv -y -c bioconda pybedtools bedtools=2.27.1
RUN conda install -n rcnv -y -c bioconda bzip2
RUN conda install -n rcnv -y -c bioconda click
RUN conda install -n rcnv -y -c bioconda mysqlclient
RUN conda install -n rcnv -y -c bioconda pybigwig
RUN conda install -n rcnv -y -c conda-forge statsmodels lxml matplotlib
ENV CONDA_DEFAULT_ENV rcnv
RUN echo "source activate rcnv" >> ~/.bashrc
ENV PATH /opt/conda/envs/rcnv/bin:$PATH
RUN pip install obonet

# Install locales
RUN apt-get -qqy update \
    && apt-get -qqy install --fix-missing locales

# Install google cloud tools
RUN pip install google-cloud-storage

# Lifted from https://github.com/GoogleCloudPlatform/cloud-sdk-docker/blob/master/debian_slim/Dockerfile
ENV CLOUD_SDK_VERSION=267.0.0
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

# Install R packages (note: intentionally do this in chunks for less painful re-builds)
RUN apt-get -qqy update \
    && apt-get -qqy install libxrender1 libcairo2-dev libxt-dev libx11-dev
RUN apt-get -qqy upgrade
RUN Rscript -e "sapply(c('optparse','metafor','KernSmooth','metap','MASS','flux'),\
    function(pkg){install.packages(pkg,repos='https://cloud.r-project.org')})"
RUN Rscript -e "sapply(c('beeswarm','Hmisc','EQL','vioplot'),\
    function(pkg){install.packages(pkg,repos='https://cloud.r-project.org')})"
RUN Rscript -e "sapply(c('psych','funr','VennDiagram'),\
    function(pkg){install.packages(pkg,repos='https://cloud.r-project.org')})"
RUN Rscript -e \
    "install.packages('https://cran.r-project.org/src/contrib/Archive/plotrix/plotrix_3.6-6.tar.gz')"

# Install Athena toolkit
RUN git clone https://c9f6da1fb0dd9d5379c313700a1083f9633be849@github.com/talkowski-lab/athena /opt/athena \
    && cd /opt/athena \
    && git checkout b66722aa3b0329d21a8f0c2446159d999c8ddd03 \
    && cd -
RUN cd /opt/athena \
    && pip install -e .

# Clone rCNV2 repo
RUN git clone -b master https://c9f6da1fb0dd9d5379c313700a1083f9633be849@github.com/talkowski-lab/rCNV2.git /opt/rCNV2 && \
    cd /opt/rCNV2 && \
    git checkout e89b402b6b2fa442a5a341c736bb3ddb6550dee7 && \
    cd -

# Developer utilities: add bash aliases, etc (for convenience)
RUN echo "source ~/.bash_aliases" >> ~/.bashrc \
    && echo 'alias l="ls -ltrha"' >> ~/.bash_aliases \
    && echo 'alias less="zless"' >> ~/.bash_aliases \
    && echo 'alias median="/opt/rCNV2/utils/median.py"' >> ~/.bash_aliases
RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    locale-gen
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8

# Launch bash
CMD ["/bin/bash"]

