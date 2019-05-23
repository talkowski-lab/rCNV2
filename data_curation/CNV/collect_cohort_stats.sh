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
gsutil cp gs://rcnv_project/analysis/analysis_refs/rCNV_metacohort_sample_counts.txt ./


# Add helper alias for formatting long integers
alias addcom="sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta'"


# Master function to collect CNV data
collect_stats () {
  infile=$1

  while read cohort N_total N_case N_ctrl bed; do
    for wrapper in 1; do
      echo "$cohort"

      if [ -s $bed ]; then
        n_case_cnv=$( zcat $bed \
                      | fgrep -v "#" | fgrep -v HEALTHY_CONTROL | wc -l )
        n_ctrl_cnv=$( zcat $bed \
                      | fgrep -v "#" | fgrep -w HEALTHY_CONTROL | wc -l )
      else
        n_case_cnv=0
        n_ctrl_cnv=0
      fi

      # Case stats
      echo "$N_case" | addcom
      if [ $n_case_cnv -gt 0 ]; then
        # Count of CNVs
        echo "$n_case_cnv" | addcom
        # Count per sample
        echo "" | awk -v n_case=$N_case -v n_cnv=$n_case_cnv \
                  '{ printf "%0.2f\n", n_cnv/n_case }'
        # Median size
        zcat $bed \
        | fgrep -v "#" \
        | fgrep -v HEALTHY_CONTROL \
        | awk '{ print $3-$2 }' \
        | sort -nk1,1 \
        | median \
        | awk '{ printf "%0.1f kb\n", $1/1000 }'
        # DEL : DUP ratio
        n_case_del=$( zcat $bed | fgrep -v "#" \
                      | fgrep -v HEALTHY_CONTROL | fgrep -w DEL | wc -l )
        n_case_dup=$( zcat $bed | fgrep -v "#" \
                      | fgrep -v HEALTHY_CONTROL | fgrep -w DUP | wc -l )
        if [ $n_case_del -ge $n_case_dup ]; then
          echo "" | awk -v del=$n_case_del -v dup=$n_case_dup \
                    '{ printf "%0.2f:1\n", del/dup }'
        else
          echo "" | awk -v del=$n_case_del -v dup=$n_case_dup \
                    '{ printf "1:%0.2f\n", dup/del }'
        fi
      else
        echo -e "0\t-\t-\t-"
      fi

      # Control stats
      echo "$N_ctrl" | addcom
      if [ $n_ctrl_cnv -gt 0 ]; then
        # Count of CNVs
        echo "$n_ctrl_cnv" | addcom
        # Count per sample
        echo "" | awk -v n_ctrl=$N_ctrl -v n_cnv=$n_ctrl_cnv \
                  '{ printf "%0.2f\n", n_cnv/n_ctrl }'
        # Median size
        zcat $bed \
        | fgrep -v "#" \
        | fgrep -w HEALTHY_CONTROL \
        | awk '{ print $3-$2 }' \
        | sort -nk1,1 \
        | median \
        | awk '{ printf "%0.1f kb\n", $1/1000 }'
        # DEL : DUP ratio
        n_ctrl_del=$( zcat $bed | fgrep -v "#" \
                      | fgrep -w HEALTHY_CONTROL | fgrep -w DEL | wc -l )
        n_ctrl_dup=$( zcat $bed | fgrep -v "#" \
                      | fgrep -w HEALTHY_CONTROL | fgrep -w DUP | wc -l )
        if [ $n_ctrl_del -ge $n_ctrl_dup ]; then
          echo "" | awk -v del=$n_ctrl_del -v dup=$n_ctrl_dup \
                    '{ printf "%0.2f:1\n", del/dup }'
        else
          echo "" | awk -v del=$n_ctrl_del -v dup=$n_ctrl_dup \
                    '{ printf "1:%0.2f\n", dup/del }'
        fi
      else
        echo -e "0\t-\t-\t-"
      fi
    done | paste -s
  done < <( fgrep -v "#" $infile ) \
  | sed -e 's/^/\|\ /g' -e 's/$/\ \|/g' -e 's/\t/\ \|\ /g'
}


# Collect raw CNV data
awk -v OFS="\t" '{ print $0, "/raw_cnv/"$1".raw.bed.gz" }' \
  /opt/rCNV2/refs/rCNV_sample_counts.txt \
| fgrep -v "#" \
> raw_cnv.input.txt
collect_stats raw_cnv.input.txt


# Collect filtered rare CNV data
awk -v OFS="\t" '{ print $0, "/cleaned_cnv/"$1".rCNV.bed.gz" }' \
  /opt/rCNV2/refs/rCNV_sample_counts.txt \
| fgrep -v "#" \
> rare_cnv.input.txt
collect_stats rare_cnv.input.txt


# Collect filtered ultra-rare CNV data
awk -v OFS="\t" '{ print $0, "/cleaned_cnv/"$1".uCNV.bed.gz" }' \
  /opt/rCNV2/refs/rCNV_sample_counts.txt \
| fgrep -v "#" \
> ultrarare_cnv.input.txt
collect_stats ultrarare_cnv.input.txt


# Collect filtered rare CNV data per metacohort
awk -v OFS="\t" '{ print $0, "/cleaned_cnv/"$1".rCNV.bed.gz" }' \
  rCNV_metacohort_sample_counts.txt \
| fgrep -v "#" \
> rare_cnv.metacohorts.input.txt
collect_stats rare_cnv.metacohorts.input.txt

# Collect filtered ultra-rare CNV data per metacohort
awk -v OFS="\t" '{ print $0, "/cleaned_cnv/"$1".uCNV.bed.gz" }' \
  rCNV_metacohort_sample_counts.txt \
| fgrep -v "#" \
> ultrarare_cnv.metacohorts.input.txt
collect_stats ultrarare_cnv.metacohorts.input.txt


