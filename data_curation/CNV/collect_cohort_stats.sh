#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Collect summary data for each array CNV cohort before and after filtering


# Launch docker image
docker run --rm -it talkowski/rcnv


# Copy all raw and filtered CNV data from Google Bucket (note: requires permissions)
gcloud auth login
mkdir ./raw_cnv
gsutil cp -r gs://rcnv_project/raw_data/cnv/* ./raw_cnv/
mkdir ./cleaned_cnv/
gsutil cp -r gs://rcnv_project/cleaned_data/cnv/* ./cleaned_cnv/


# Collect raw CNV data
while read cohort N_total N_case N_ctrl; do
  for wrapper in 1; do
    echo "$cohort"

    echo "$N_case"
    n_case_cnv=$( zcat raw_cnv/$cohort.raw.bed.gz \
                  | fgrep -v "#" | fgrep -v CTRL | wc -l )
    if [ -s raw_cnv/$cohort.raw.bed.gz ] && [ $n_case_cnv -gt 0 ]; then
      echo "$n_case_cnv"
      zcat raw_cnv/$cohort.raw.bed.gz \
      | fgrep -v "#" \
      | awk '{ print $3-$2 }' \
      | sort -nk1,1 \
      | median
      n_case_del=$( zcat raw_cnv/$cohort.raw.bed.gz | fgrep -v "#" \
                    | fgrep -v CTRL | fgrep -w DEL | wc -l )
      n_case_dup=$( zcat raw_cnv/$cohort.raw.bed.gz | fgrep -v "#" \
                    | fgrep -v CTRL | fgrep -w DUP | wc -l )

    else
      echo -e "0\t-\t-\t-"
    fi

    # echo "$N_ctrl"
    # n_ctrl_cnv=$( zcat raw_cnv/$cohort.raw.bed.gz \
    #               | fgrep -v "#" | fgrep -w CTRL | wc -l )
    # if [ -s raw_cnv/$cohort.raw.bed.gz ] && [ $n_ctrl_cnv -gt 0 ]; then
    #   echo YES
    # else
    #   echo -e "0\t-\t-\t-"
    # fi
  done | paste -s
done < <( fgrep -v "#" /opt/rCNV2/refs/rCNV_sample_counts.txt )
