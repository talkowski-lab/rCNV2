#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control CNV burdens in sliding windows, genome-wide


# Launch docker image
docker run --rm -it talkowski/rcnv


# Copy all filtered CNV data, sliding windows, and other references 
# from the project Google Bucket (note: requires permissions)
gcloud auth login
mkdir cleaned_cnv/
gsutil cp -r gs://rcnv_project/cleaned_data/cnv/* cleaned_cnv/
mkdir windows/
gsutil cp -r gs://rcnv_project/cleaned_data/binned_genome/* windows/
mkdir refs/
gsutil cp gs://rcnv_project/analysis/analysis_refs/* refs/


# Test/dev parameters (neurodevelopmental abnormality)
hpo="HP:0012759"
prefix="HP0012759"
meta="meta1"
freq_code="rCNV"
metacohort_list="refs/rCNV_metacohort_list.txt"
metacohort_sample_table="refs/HPOs_by_metacohort.table.tsv"
binned_genome="windows/GRCh37.100kb_bins_10kb_steps.raw.bed.gz"
p_cutoff=0.00001


# Count CNVs in cases and controls per phenotype, split by metacohort and CNV type
# Iterate over phenotypes
while read prefix hpo; do
  descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" '{ print $2 }' )
  # Iterate over metacohorts
  while read meta cohorts; do
    cnv_bed="cleaned_cnv/$meta.$freq_code.bed.gz"
    meta_idx=$( head -n1 "${metacohort_sample_table}" \
                | sed 's/\t/\n/g' \
                | awk -v meta="$meta" '{ if ($1==meta) print NR }' )
    ncase=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" -v meta_idx="$meta_idx" '{ print $meta_idx }' \
             | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
    nctrl=$( fgrep -w "HEALTHY_CONTROL" "${metacohort_sample_table}" \
             | awk -v FS="\t" -v meta_idx="$meta_idx" '{ print $meta_idx }' \
             | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
    title="$descrip (${hpo})\n$ncase cases vs $nctrl controls in '${meta}' cohort"

    # Iterate over CNV types
    for CNV in CNV DEL DUP; do
      # Set CNV-specific parameters
      case "$CNV" in
        DEL)
          highlight_bed=/opt/rCNV2/refs/UKBB_GD.Owen_2018.DEL.bed.gz
          ;;
        DUP)
          highlight_bed=/opt/rCNV2/refs/UKBB_GD.Owen_2018.DUP.bed.gz
          ;;
        *)
          highlight_bed=/opt/rCNV2/refs/UKBB_GD.Owen_2018.bed.gz
          ;;
      esac

      # Count CNVs
      /opt/rCNV2/analysis/sliding_windows/count_cnvs_per_window.py \
        -t $CNV \
        --hpo ${hpo} \
        -z \
        -o "$meta.${prefix}.${freq_code}.$CNV.sliding_window.counts.bed.gz" \
        ${cnv_bed} \
        ${binned_genome}

      # Perform burden test
      /opt/rCNV2/analysis/sliding_windows/window_burden_test.R \
        --pheno-table refs/HPOs_by_metacohort.table.tsv \
        --cohort-name $meta \
        --case-hpo ${hpo} \
        --bgzip \
        "$meta.${prefix}.${freq_code}.$CNV.sliding_window.counts.bed.gz" \
        "$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz"

      # Generate Manhattan & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --p-col-name "fisher_phred_p" \
        --p-is-phred \
        --cutoff ${p_cutoff} \
        --highlight-bed "$highlight_bed" \
        --highlight-name "Known GDs (Owen 2018)" \
        --title "$title" \
        "$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz" \
        "$meta.${prefix}.${freq_code}.$CNV.sliding_window"
    done

    # Generate Miami & QQ plots
    ncase=
    nctrl=
    title="$descrip (${hpo})\nCohort '${meta}': $ncase cases vs $nctrl controls"
    /opt/rCNV2/utils/plot_manhattan_qq.R \
      --miami \
      --p-col-name "fisher_phred_p" \
      --p-is-phred \
      --cutoff ${p_cutoff} \
      --highlight-bed /opt/rCNV2/refs/UKBB_GD.Owen_2018.DUP.bed.gz \
      --highlight-name "Known DUP GDs (Owen 2018)" \
      --label-prefix "DUP" \
      --highlight-bed-2 /opt/rCNV2/refs/UKBB_GD.Owen_2018.DEL.bed.gz \
      --highlight-name-2 "Known DEL GDs (Owen 2018)" \
      --label-prefix-2 "DEL" \
      --title "$title" \
      "$meta.${prefix}.${freq_code}.DUP.sliding_window.stats.bed.gz" \
      "$meta.${prefix}.${freq_code}.DEL.sliding_window.stats.bed.gz" \
      "$meta.${prefix}.${freq_code}.sliding_window"
  done < ${metacohort_list}
done < refs/test_phenotypes.list


# Run burden tests - example
/opt/rCNV2/analysis/sliding_windows/window_burden_test.R \
  --bgzip \
  --case-column Coe.rCNV.CASE \
  --case-n 29083 \
  --control-column UKBB.rCNV.CTRL \
  --control-n 480501 \
  GRCh37.100kb_bins_10kb_steps.raw.rCNV_counts.DEL.bed.gz \
  GRCh37.100kb_bins_10kb_steps.DEL.test_stats.bed.gz


