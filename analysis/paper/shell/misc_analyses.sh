#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Miscellaneous helper analyses for rCNV2 manuscript


# Launch docker image & authenticate GCP credentials
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d1"
export finemap_dist=1000000


# Download necessary data (note: requires permissions)
gsutil -m cp \
  ${rCNV_bucket}/results/gene_association/* \
  ${rCNV_bucket}/results/segment_association/* \
  ${rCNV_bucket}/analysis/gene_burden/fine_mapping/rCNV.*.gene_fine_mapping.gene_stats.merged_no_variation_features.all_genes_from_blocks.tsv \
  ./


# Compute total number of significant loci
for CNV in DEL DUP; do
  zcat rCNV.final_segments.loci.bed.gz \
       rCNV.final_genes.credible_sets.bed.gz \
  | fgrep -w $CNV \
  | cut -f1-3 \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  | wc -l
done


# Compute total number of significant associations
for CNV in DEL DUP; do
  zcat rCNV.final_segments.associations.bed.gz \
       rCNV.final_genes.credible_sets.bed.gz \
  | fgrep -w $CNV \
  | awk -v FS="\t" -v OFS="\t" '{ print $6"_"$1, $2, $3 }' \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  | wc -l
done
