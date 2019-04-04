#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Segment the reference into sequential regular bins and add key annotations


# Launch docker image
docker run --rm -it talkowski/rcnv


# Download necessary references (note: requires permissions)
gcloud auth login
gsutil cp -r gs://rcnv_project/refs ./


# Set desired bin resolution (in kb)
binsize=50


# Create 50kb sequential bins
athena make-bins -z \
	-x refs/GRCh37.Nmask.bed.gz \
	--buffer "$binsize"000 \
	refs/GRCh37.autosomes.genome \
	"$binsize"000 \
	GRCh37."$binsize"kb_bins.cleaned.bed.gz
gsutil cp GRCh37."$binsize"kb_bins.cleaned.bed.gz \
  gs://rcnv_project/cleaned_data/binned_genome/


# Annotate 