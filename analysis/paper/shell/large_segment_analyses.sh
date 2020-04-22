#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master code block for secondary analyses of large rCNV-associated segments


# Launch docker image & authenticate GCP credentials
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export plot_prefix="rCNV2_analysis_d1"


# Localize all analysis refs, sliding window meta-analysis stats, and large segment results
mkdir refs/ 
gsutil -m cp \
   ${rCNV_bucket}/analysis/analysis_refs/** \
   refs/
mkdir meta_stats/
gsutil -m cp \
  ${rCNV_bucket}/analysis/sliding_windows/**.rCNV.**.sliding_window.meta_analysis.stats.bed.gz \
  meta_stats/
gsutil -m cp \
  ${rCNV_bucket}/results/segment_association/* \
  ./


# Tabix all meta-analysis stats
find meta_stats/ -name "*meta_analysis.stats.bed.gz" | xargs -I {} tabix -f {}


# Compute effect size and max P-value per phenotype per final segment
/opt/rCNV2/analysis/paper/scripts/large_segments/calc_all_seg_stats.py \
  -o rCNV.final_segments.loci.all_sumstats.tsv \
  rCNV.final_segments.loci.bed.gz \
  refs/test_phenotypes.list \
  meta_stats
gzip -f rCNV.final_segments.loci.all_sumstats.tsv


# Collapse overlapping DEL/DUP segments for sake of plotting
while read intervals rid; do
  echo "$intervals" | sed -e 's/\;/\n/g' -e 's/\:\|\-/\t/g' \
  | awk -v OFS="\t" -v rid=$rid '{ print $0, rid }'
done < <( zcat rCNV.final_segments.loci.bed.gz | grep -ve '^#' \
          | awk -v FS="\t" -v OFS="\t" '{ print $(NF-3), $4 }' ) \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
| bedtools merge -i - -d 1000000 -c 4 -o distinct \
| cut -f4 \
> locus_clusters.txt


# Plot master grid summarizing segment association across all phenotypes
cut -f2 refs/test_phenotypes.list > phenotypes.list
/opt/rCNV2/analysis/paper/plot/large_segments/plot_association_grid.R \
  --clusters locus_clusters.txt \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.final_segments.loci.bed.gz \
  rCNV.final_segments.loci.all_sumstats.tsv.gz \
  phenotypes.list \
  ${plot_prefix}.large_segments.association_grid.pdf

