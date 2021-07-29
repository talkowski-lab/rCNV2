# rCNV Project: Disease association analyses of rare copy number variation
# Dockerfile
# Copyright (c) 2019-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# Distributed under terms of the MIT License

# Base image: miniconda 3
FROM continuumio/miniconda3
MAINTAINER "Ryan Collins <rlcollins@g.harvard.edu>"

# Install make, zlib, unzip, and locales
RUN apt-get -qqy update && \
    apt-get -qqy install --fix-missing build-essential zlib1g-dev unzip locales

# Install other system libraries required for R packages
RUN apt-get -qqy update && \
    apt-get -qqy install software-properties-common && \
    add-apt-repository contrib
RUN apt-get -qqy update && \
    apt-get -qqy install ttf-mscorefonts-installer

# Install UCSC tools (liftOver & bigWigToBedGraph)
RUN wget -P /opt/ http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver && \
    chmod a+x /opt/liftOver
RUN wget -P /opt/ http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigWigToBedGraph && \
    chmod a+x /opt/bigWigToBedGraph
ENV PATH /opt:$PATH

# Developer utilities: add bash aliases, etc (for convenience)
RUN echo "source ~/.bash_aliases" >> ~/.bashrc && \
    echo 'alias l="ls -ltrha"' >> ~/.bash_aliases && \
    echo 'alias less="zless"' >> ~/.bash_aliases && \
    echo 'alias median="/opt/rCNV2/utils/median.py"' >> ~/.bash_aliases
RUN sed -i -e 's/# en_US.UTF-8 UTF-8/en_US.UTF-8 UTF-8/' /etc/locale.gen && \
    locale-gen
ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8

# Clone rCNV2 repo
ARG rCNV2_commit_hash="master"
RUN git clone https://github.com/talkowski-lab/rCNV2.git /opt/rCNV2 && \
    cd /opt/rCNV2 && \
    git checkout "$rCNV2_commit_hash" && \
    cd -

# Create python environment from .yml file
# Note: install packages from two separate .yml files to make Docker layers more manageable
RUN conda update -n base conda && \
    conda clean -y --all
RUN conda env update --file /opt/rCNV2/config/conda_env_config.yml && \
    conda clean -y --all
RUN conda env update --file /opt/rCNV2/config/conda_env_config.r_packages.yml && \
    conda clean -y --all
ENV CONDA_DEFAULT_ENV base
RUN echo "source activate base" >> ~/.bashrc
ENV PATH /opt/conda/bin:$PATH

# Install google cloud tools
RUN pip install google-cloud-storage
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get -qqy clean && \
    apt-get -qqy update --fix-missing && \
    apt-get -qqy dist-upgrade && \
    apt-get -qqy install lsb-release && \
    pip install -q -U --no-cache-dir \
                crcmod \
                wheel && \
    export CLOUD_SDK_REPO="cloud-sdk-$( lsb_release -c -s )" && \
    echo "deb https://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" > /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl -s https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get -qqy update && \
    apt-get -qqy install --no-install-recommends \
                 google-cloud-sdk && \
    gcloud config set core/disable_usage_reporting true && \
    gcloud config set component_manager/disable_update_check true && \
    gcloud config set metrics/environment github_docker_image && \
    apt-get -qqy clean && \
    rm -rf /tmp/* \
           /var/tmp/* \
           /var/cache/apt/* \
           /var/lib/apt/lists/* \
           /usr/share/man/?? \
           /usr/share/man/??_* && \
    conda clean -y --all && \
    gcloud --help

# Install R packages not available through conda
RUN Rscript -e "sapply(c('EQL'),\
    function(pkg){install.packages(pkg,repos='https://cloud.r-project.org')})"
ARG rCNV2_Rpkg_source="/opt/rCNV2/source/rCNV2_0.1.0.tar.gz"
RUN Rscript -e "install.packages('$rCNV2_Rpkg_source', source=TRUE, repos=NULL)"

# Launch bash
CMD ["/bin/bash"]

