#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master code block for secondary analyses of noncoding association for rCNV2 manuscript


# Launch docker image & authenticate GCP credentials
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d1"


# Download necessary data (note: requires permissions)
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/genome_annotations/rCNV.burden_stats.tsv.gz \
  ${rCNV_bucket}/cleaned_data/genome_annotations/rCNV.crbs.bed.gz \
  ./
mkdir refs/
gsutil -m cp -r \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists \
  ${rCNV_bucket}/refs/REP_state_manifest.tsv \
  refs/


# # Plot observed effect size distributions split by gene set membership
# # Used to justify inclusion of some coding effects in noncoding association test
# fgrep -wvf \
#   refs/gene_lists/HP0000118.HPOdb.genes.list \
#   refs/gene_lists/gnomad.v2.1.1.likely_unconstrained.genes.list \
# | sort -V | uniq \
# > loose_noncoding_whitelist.genes.list
# /opt/rCNV2/analysis/paper/plot


# Plot annotation track distributions & stats
/opt/rCNV2/analysis/paper/plot/noncoding_association/plot_track_stats.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.burden_stats.tsv.gz \
  ${prefix}


# Volcano plots of DEL & DUP track-level burden tests
/opt/rCNV2/analysis/paper/plot/noncoding_association/plot_volcanos.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.burden_stats.tsv.gz \
  ${prefix}


# Plot ChromHMM enrichments as positive controls
/opt/rCNV2/analysis/paper/plot/noncoding_association/chromhmm_enrichments.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.burden_stats.tsv.gz \
  refs/REP_state_manifest.tsv \
  ${prefix}


# Plot CRB distributions & stats
/opt/rCNV2/analysis/paper/plot/noncoding_association/plot_crb_stats.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.crbs.bed.gz \
  ${prefix}


# Copy all plots to final gs:// directory
gsutil -m cp -r \
  ${prefix}.track_stats.*.pdf \
  ${prefix}.*volcano*pdf \
  ${prefix}.chromhmm_enrichment.*.pdf \
  ${prefix}.crb_stats.*.pdf \
  ${rCNV_bucket}/analysis/paper/plots/noncoding_association/
