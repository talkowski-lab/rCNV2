#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Segment the reference into sliding windows and add key annotations


# Launch docker image
docker run --rm -it us.gcr.io/broad-dsmap/athena-cloud:latest


# Download necessary references (note: requires permissions)
gcloud auth login
mkdir refs/
gsutil cp -r \
  gs://rcnv_project/refs/GRCh37.Nmask.autosomes.bed.gz \
  gs://rcnv_project/refs/GRCh37.somatic_hypermutable_sites.bed.gz \
  gs://rcnv_project/refs/GRCh37.autosomes.genome \
  gs://rcnv_project/analysis/paper/data/large_segments/wgs_common_cnvs.*.1pct.bed.gz \
  refs/
gsutil cp -r gs://rcnv_project/cleaned_data/control_probesets ./


# Prep reference fasta
# wget ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
# gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
# samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa
# gsutil cp Homo_sapiens.GRCh37.dna.primary_assembly.fa gs://rcnv_project/GRCh37_ref_build/GRCh37.primary_assembly.fa
# gsutil cp Homo_sapiens.GRCh37.dna.primary_assembly.fa.fai gs://rcnv_project/GRCh37_ref_build/GRCh37.primary_assembly.fa.fai
# gsutil cp -r gs://rcnv_project/GRCh37_ref_build/* refs/


# Set desired bin size & step size resolution (in kb)
binsize=200
stepsize=10


# Create sliding windows
athena make-bins -z \
	-x refs/GRCh37.Nmask.autosomes.bed.gz \
  -x refs/GRCh37.somatic_hypermutable_sites.bed.gz \
  -x refs/wgs_common_cnvs.DEL.1pct.bed.gz \
  -x refs/wgs_common_cnvs.DUP.1pct.bed.gz \
  -s ${stepsize}000 \
  --exclusion-list-cov 0.3 \
	refs/GRCh37.autosomes.genome \
	${binsize}000 \
	GRCh37.${binsize}kb_bins_${stepsize}kb_steps.raw.bed.gz
tabix -f GRCh37.${binsize}kb_bins_${stepsize}kb_steps.raw.bed.gz
gsutil cp GRCh37.${binsize}kb_bins_${stepsize}kb_steps.raw.bed.gz* \
  gs://rcnv_project/cleaned_data/binned_genome/

