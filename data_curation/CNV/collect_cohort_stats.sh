#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Collect and plot summary data for each array CNV cohort before and after filtering


# Launch docker image
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv:latest


# Copy all raw and filtered CNV data from Google Bucket (note: requires permissions)
gcloud auth login
mkdir ./raw_cnv
gsutil -m cp -r gs://rcnv_project/raw_data/cnv/* ./raw_cnv/
mkdir ./cleaned_cnv
gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/* ./cleaned_cnv/
gsutil cp gs://rcnv_project/analysis/analysis_refs/rCNV_metacohort_sample_counts.txt ./


# Add helper alias for formatting long integers
alias addcom="sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta'"


# Collect and plot raw CNV data
awk -v OFS="\t" '{ print $0, "/raw_cnv/"$1".raw.bed.gz" }' \
  /opt/rCNV2/refs/rCNV_sample_counts.txt \
| fgrep -v "#" \
> raw_cnv.input.txt
/opt/rCNV2/data_curation/CNV/build_cnv_stats_table.py \
  --tsv raw_cnv.stats.txt \
  --html raw_cnv.stats.html.txt \
  raw_cnv.input.txt
/opt/rCNV2/data_curation/CNV/plot_cnv_stats_per_cohort.R \
  raw_cnv.stats.txt \
  raw_cnv.stats.jpg


# Iterate over frequency classes
for freq in rCNV vCNV uCNV; do
  # Collect and plot filtered CNV data per cohort
  awk -v OFS="\t" -v freq=${freq} \
    '{ print $0, "/cleaned_cnv/"$1"."freq".bed.gz" }' \
    /opt/rCNV2/refs/rCNV_sample_counts.txt \
  | fgrep -v "#" \
  > ${freq}.input.txt
  /opt/rCNV2/data_curation/CNV/build_cnv_stats_table.py \
    --tsv ${freq}.stats.txt \
    --html ${freq}.stats.html.txt \
    ${freq}.input.txt
  /opt/rCNV2/data_curation/CNV/plot_cnv_stats_per_cohort.R \
    ${freq}.stats.txt \
    ${freq}.stats.jpg

  # Collect and plot filtered CNV data per metacohort
  awk -v OFS="\t" -v freq=${freq} \
    '{ print $0, "/cleaned_cnv/"$1"."freq".bed.gz" }' \
    rCNV_metacohort_sample_counts.txt \
  | fgrep -v "#" \
  > ${freq}.metacohorts.input.txt
  /opt/rCNV2/data_curation/CNV/build_cnv_stats_table.py \
    --tsv ${freq}.metacohort.stats.txt \
    --html ${freq}.metacohort.stats.html.txt \
    ${freq}.metacohorts.input.txt
  /opt/rCNV2/data_curation/CNV/plot_cnv_stats_per_cohort.R \
    ${freq}.metacohort.stats.txt \
    ${freq}.metacohort.stats.jpg
done


# Collect and plot filtered noncoding CNV subsets per metacohort
for subset in strict loose; do
  awk -v OFS="\t" -v subset=$subset \
    '{ print $0, "/cleaned_cnv/noncoding/"$1".rCNV."subset"_noncoding.bed.gz" }' \
    rCNV_metacohort_sample_counts.txt \
  | fgrep -v "#" | fgrep -v mega \
  > ${subset}_noncoding.metacohorts.input.txt
  /opt/rCNV2/data_curation/CNV/build_cnv_stats_table.py \
    --tsv ${subset}_noncoding.metacohort.stats.txt \
    --html ${subset}_noncoding.metacohort.stats.html.txt \
    ${subset}_noncoding.metacohorts.input.txt
  /opt/rCNV2/data_curation/CNV/plot_cnv_stats_per_cohort.R \
    ${subset}_noncoding.metacohort.stats.txt \
    ${subset}_noncoding.metacohort.stats.jpg
done


# Copy all plots to public bucket for viewing on README
gsutil -m cp ./*jpg gs://rcnv_project/public/
gsutil acl ch -u AllUsers:R gs://rcnv_project/public/*.jpg

