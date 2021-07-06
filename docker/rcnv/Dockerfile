# rCNV Project: Disease association analyses of rare copy number variation
# Dockerfile
# Copyright (c) 2019-Present Ryan L. Collins <rlcollins@g.harvard.edu> 
# Distributed under terms of the MIT License

# Base image: rcnv-base
FROM gcr.io/gnomad-wgs-v2-sv/rcnv:latest
MAINTAINER "Ryan Collins <rlcollins@g.harvard.edu>"

# Checkout most recent branch of rCNV2 repo and pull updates
ARG rCNV2_commit_hash="master"
RUN cd /opt/rCNV2 && \
	git checkout "$rCNV2_commit_hash" && \
	git pull && \
	cd -

# Reinstall rCNV R library
ARG rCNV2_Rpkg_source="/opt/rCNV2/source/rCNV2_0.1.0.tar.gz"
RUN Rscript -e "install.packages('$rCNV2_Rpkg_source', source=TRUE, repos=NULL)"

# Launch bash
CMD ["/bin/bash"]
