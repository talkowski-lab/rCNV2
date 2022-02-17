#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2022-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Code to prioritize loci to highlight in rCNV2 manuscript


# Launch docker image & authenticate GCP credentials
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d2"


# Download necessary data (note: requires permissions)
gsutil -m cp -r \
  ${rCNV_bucket}/results/gene_association/* \
  ${rCNV_bucket}/results/segment_association/* \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists \
  ${rCNV_bucket}/analysis/gene_scoring/gene_lists \
  ./
mkdir refs
gsutil -m cp -r \
  ${rCNV_bucket}/analysis/analysis_refs/** \
  ${rCNV_bucket}/analysis/paper/data/large_segments/lit_GDs*bed.gz \
  refs/


# Pre-filtering step: exclude all segments & credsets matching known GDs
for CNV in DEL DUP; do
  zcat segment_association/rCNV.final_segments.loci.bed.gz | fgrep -w $CNV \
  | bedtools intersect -v -r -f 0.2 -a - \
      -b <( zcat refs/lit_GDs.*.bed.gz | fgrep -w $CNV ) \
  | bgzip -c > segments.$CNV.prefiltered.bed.gz
  zcat segment_association/rCNV.final_genes.credible_sets.bed.gz | fgrep -w $CNV \
  | bedtools intersect -v -r -f 0.2 -a - \
      -b <( zcat refs/lit_GDs.*.bed.gz | fgrep -w $CNV ) \
  | bgzip -c > segments.$CNV.credsets.bed.gz
done 


# LoF or missense-constrained genes with no known disease associations
for CNV in DEL DUP; do
  zcat segments.$CNV.credsets.bed.gz | fgrep -w $CNV | cut -f17
done


# Known haploinsufficient genes with duplication associations




