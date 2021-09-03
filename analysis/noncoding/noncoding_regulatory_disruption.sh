#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2021 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Regulatory disruption scoring for noncoding CNVs


# Launch docker image
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv:latest
gcloud auth login


# Copy noncoding CNV data and other references from the project Google Bucket (note: requires permissions)
mkdir cleaned_cnv/
gsutil -m cp -r \
  gs://rcnv_project/cleaned_data/cnv/noncoding/* \
  gs://rcnv_project/cleaned_data/cnv/*rCNV.bed.gz* \
  cleaned_cnv/
mkdir refs/
gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/


# Clone regulatory disruption scoring repo
cd /opt \
&& git clone https://github.com/RuderferLab/CNV_FunctionalAnnotation.git \
&& cd -


# Build regulatory disruption scoring references
cd opt/CNV_FunctionalAnnotation \
&& bin/make_annotation_files \
&& cd -


# Apply regulatory disruption model to noncoding CNVs per cohort
if [ -e /reg_scoring ]; then
  rm -rf /reg_scoring
fi
mkdir /reg_scoring
mkdir /reg_scoring/reg_scoring_results
while read meta cohorts; do

  echo -e "\n\n\nStarting ${meta}...\n"

  # Reformat CNV file
  zcat cleaned_cnv/$meta.rCNV.loose_noncoding.bed.gz \
  | fgrep -v "#" \
  | awk -v FS="\t" -v OFS="\t" '{ print $4, $1, $2, $3, $5 }' \
  | cat <( echo -e "id\tchr\tstart\tend\tSVType" ) - \
  > reg_scoring/$meta.scoring_input.tsv

  # Run scoring model
  sh /opt/CNV_FunctionalAnnotation/bin/annotate-cnv.sh \
     /opt/ \
     reg_scoring/$meta.scoring_input.tsv

  # Prep outputs for analysis
  sed '1d' reg_scoring/$meta.scoring_input.tsv.reg_disruption \
  | awk -v OFS="\t" '{ print $1, $(NF-1), $NF }' \
  | sort -Vk1,1 \
  | join -j 1 -t $'\t' - \
      <( zcat cleaned_cnv/$meta.rCNV.loose_noncoding.bed.gz \
         | sed '1d' | cut -f4-6 | sort -Vk1,1 ) \
  | cat <( echo -e "#id\tpred_exp\treg_dist\tCNV\thpos" ) - \
  | gzip -c \
  > /reg_scoring/$meta.reg_info.all_hpos.tsv.gz

  # Split by CNV type and HPO
  for CNV in DEL DUP; do
    # Controls
    zcat /reg_scoring/$meta.reg_info.all_hpos.tsv.gz \
    | fgrep -w HEALTHY_CONTROL | fgrep -w $CNV \
    | cut -f2-3 \
    | cat <( echo -e "pred_exp\treg_dist" ) - \
    | gzip -c \
    > /reg_scoring/reg_scoring_results/$meta.reg_info.$CNV.controls.tsv.gz

    # Developmental cases
    zcat /reg_scoring/$meta.reg_info.all_hpos.tsv.gz \
    | fgrep -wf refs/rCNV2.hpos_by_severity.developmental.list \
    | fgrep -w $CNV \
    | cut -f2-3 \
    | cat <( echo -e "pred_exp\treg_dist" ) - \
    | gzip -c \
    > /reg_scoring/reg_scoring_results/$meta.reg_info.$CNV.developmental.tsv.gz

    # Adult-onset cases
    zcat /reg_scoring/$meta.reg_info.all_hpos.tsv.gz \
    | fgrep -wvf refs/rCNV2.hpos_by_severity.developmental.list \
    | fgrep -wf refs/rCNV2.hpos_by_severity.adult.list \
    | fgrep -w $CNV \
    | cut -f2-3 \
    | cat <( echo -e "pred_exp\treg_dist" ) - \
    | gzip -c \
    > /reg_scoring/reg_scoring_results/$meta.reg_info.$CNV.adult.tsv.gz
  done

done < <( fgrep -v mega refs/rCNV_metacohort_list.txt )


# Plot distribution of scores
for CNV in DEL DUP; do
  /opt/rCNV2/analysis/noncoding/plot_reg_disruption_dist.R \
    --cnv $CNV \
    /reg_scoring/reg_scoring_results \
    /rCNV2.loose_noncoding_$CNV
done


# Apply regulatory disruption model to all CNVs (including coding CNVs) per cohort
if [ -e /reg_scoring_wCoding ]; then
  rm -rf /reg_scoring_wCoding
fi
mkdir /reg_scoring_wCoding
mkdir /reg_scoring_wCoding/reg_scoring_wCoding_results
while read meta cohorts; do

  echo -e "\n\n\nStarting ${meta}...\n"

  # Reformat CNV file
  zcat cleaned_cnv/$meta.rCNV.bed.gz \
  | fgrep -v "#" \
  | awk -v FS="\t" -v OFS="\t" '{ print $4, $1, $2, $3, $5 }' \
  | cat <( echo -e "id\tchr\tstart\tend\tSVType" ) - \
  > reg_scoring_wCoding/$meta.scoring_input.tsv

  # Run scoring model
  sh /opt/CNV_FunctionalAnnotation/bin/annotate-cnv.sh \
     /opt/ \
     reg_scoring_wCoding/$meta.scoring_input.tsv

  # Prep outputs for analysis
  sed '1d' reg_scoring_wCoding/$meta.scoring_input.tsv.reg_disruption \
  | awk -v OFS="\t" '{ print $1, $(NF-1), $NF }' \
  | sort -Vk1,1 \
  | join -j 1 -t $'\t' - \
      <( zcat cleaned_cnv/$meta.rCNV.bed.gz \
         | sed '1d' | cut -f4-6 | sort -Vk1,1 ) \
  | cat <( echo -e "#id\tpred_exp\treg_dist\tCNV\thpos" ) - \
  | gzip -c \
  > /reg_scoring_wCoding/$meta.reg_info.all_hpos.tsv.gz

  # Split by CNV type and HPO
  for CNV in DEL DUP; do
    # Controls
    zcat /reg_scoring_wCoding/$meta.reg_info.all_hpos.tsv.gz \
    | fgrep -w HEALTHY_CONTROL | fgrep -w $CNV \
    | cut -f2-3 \
    | cat <( echo -e "pred_exp\treg_dist" ) - \
    | gzip -c \
    > /reg_scoring_wCoding/reg_scoring_wCoding_results/$meta.reg_info.$CNV.controls.tsv.gz

    # Developmental cases
    zcat /reg_scoring_wCoding/$meta.reg_info.all_hpos.tsv.gz \
    | fgrep -wf refs/rCNV2.hpos_by_severity.developmental.list \
    | fgrep -w $CNV \
    | cut -f2-3 \
    | cat <( echo -e "pred_exp\treg_dist" ) - \
    | gzip -c \
    > /reg_scoring_wCoding/reg_scoring_wCoding_results/$meta.reg_info.$CNV.developmental.tsv.gz

    # Adult-onset cases
    zcat /reg_scoring_wCoding/$meta.reg_info.all_hpos.tsv.gz \
    | fgrep -wvf refs/rCNV2.hpos_by_severity.developmental.list \
    | fgrep -wf refs/rCNV2.hpos_by_severity.adult.list \
    | fgrep -w $CNV \
    | cut -f2-3 \
    | cat <( echo -e "pred_exp\treg_dist" ) - \
    | gzip -c \
    > /reg_scoring_wCoding/reg_scoring_wCoding_results/$meta.reg_info.$CNV.adult.tsv.gz
  done

done < <( fgrep -v mega refs/rCNV_metacohort_list.txt )


# Plot distribution of scores
for CNV in DEL DUP; do
  /opt/rCNV2/analysis/noncoding/plot_reg_disruption_dist.R \
    --cnv $CNV \
    --prefix All \
    /reg_scoring_wCoding/reg_scoring_wCoding_results \
    /rCNV2.all_$CNV
done

