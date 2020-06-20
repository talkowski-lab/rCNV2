#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control noncoding CNV burdens for cis-regulatory blocks (CRBs)


# Launch docker image
docker run --rm -it talkowski/rcnv
gcloud auth login


# Copy all filtered CNV data, gene coordinates, and other references 
# from the project Google Bucket (note: requires permissions)
mkdir cleaned_cnv/
gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/noncoding/* cleaned_cnv/
gsutil -m cp -r gs://rcnv_project/cleaned_data/genome_annotations/*bed.gz* ./
mkdir refs/
gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/


# Test/dev parameters
# Params for NDDs
hpo="HP:0012759"
prefix="HP0012759"
meta="meta1"
# General params
freq_code="rCNV"
noncoding_filter="loose"
CNV="DEL"
phenotype_list="refs/test_phenotypes.list"
metacohort_list="refs/rCNV_metacohort_list.txt"
metacohort_sample_table="refs/HPOs_by_metacohort.table.tsv"
crbs="rCNV.crbs.bed.gz"
crb_elements="rCNV.crb_elements.bed.gz"
rCNV_bucket="gs://rcnv_project"
pad_controls=0
min_element_ovr=1.0
min_frac_all_elements=0.05
p_cutoff=0.000003741955
max_manhattan_phred_p=30
n_pheno_perms=50
meta_model_prefix="fe"
i=1


# Count CNVs in cases and controls per phenotype, split by metacohort and CNV type
# Iterate over phenotypes
while read pheno hpo; do
  # Set HPO-specific parameters
  descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" '{ print $2 }' )

  # Iterate over metacohorts
  while read meta cohorts; do
    echo $meta

    # Set metacohort-specific parameters
    cnv_bed="cleaned_cnv/$meta.${freq_code}.${noncoding_filter}_noncoding.bed.gz"
    meta_idx=$( head -n1 "${metacohort_sample_table}" \
                | sed 's/\t/\n/g' \
                | awk -v meta="$meta" '{ if ($1==meta) print NR }' )
    ncase=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" -v meta_idx="$meta_idx" '{ print $meta_idx }' \
             | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
    nctrl=$( fgrep -w "HEALTHY_CONTROL" "${metacohort_sample_table}" \
             | awk -v FS="\t" -v meta_idx="$meta_idx" '{ print $meta_idx }' \
             | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
    title="$descrip (${hpo})\n$ncase cases vs $nctrl controls in '$meta' cohort"

    # Iterate over CNV types
    for CNV in DEL DUP; do
      echo $CNV

      # Count CNVs
      /opt/rCNV2/analysis/noncoding/count_cnvs_per_crb.py \
        --cnvs $cnv_bed \
        --crbs ${crbs} \
        --elements ${crb_elements} \
        --pad-controls ${pad_controls} \
        --min-element-ovr ${min_element_ovr} \
        --min-frac-all-elements ${min_frac_all_elements} \
        -t $CNV \
        --hpo ${hpo} \
        -z \
        -o "$meta.${prefix}.${freq_code}.$CNV.noncoding_association.counts.bed.gz" 
      tabix -f "$meta.${prefix}.${freq_code}.$CNV.noncoding_association.counts.bed.gz"

#       # Perform burden test
#       /opt/rCNV2/analysis/genes/gene_burden_test.R \
#         --pheno-table ${metacohort_sample_table} \
#         --cohort-name $meta \
#         --cnv $CNV \
#         --case-hpo ${hpo} \
#         --bgzip \
#         "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.bed.gz" \
#         "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz"
#       tabix -f "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz"

#       # Generate Manhattan & QQ plots
#       /opt/rCNV2/utils/plot_manhattan_qq.R \
#         --p-col-name "fisher_phred_p" \
#         --p-is-phred \
#         --max-phred-p ${max_manhattan_phred_p} \
#         --cutoff ${p_cutoff} \
#         --highlight-bed "${prefix}.highlight_regions.bed" \
#         --highlight-name "Constrained genes associated with this phenotype" \
#         --label-prefix "$CNV" \
#         --title "$title" \
#         "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz" \
#         "$meta.${prefix}.${freq_code}.$CNV.gene_burden"
#     done

#     # Generate Miami & QQ plots
#     /opt/rCNV2/utils/plot_manhattan_qq.R \
#       --miami \
#       --p-col-name "fisher_phred_p" \
#       --p-is-phred \
#       --max-phred-p ${max_manhattan_phred_p} \
#       --cutoff ${p_cutoff} \
#       --highlight-bed "${prefix}.highlight_regions.bed" \
#       --highlight-name "Constrained genes associated with this phenotype" \
#       --label-prefix "DUP" \
#       --highlight-bed-2 "${prefix}.highlight_regions.bed" \
#       --highlight-name-2 "Constrained genes associated with this phenotype" \
#       --label-prefix-2 "DEL" \
#       --title "$title" \
#       "$meta.${prefix}.${freq_code}.DUP.gene_burden.stats.bed.gz" \
#       "$meta.${prefix}.${freq_code}.DEL.gene_burden.stats.bed.gz" \
#       "$meta.${prefix}.${freq_code}.gene_burden"
#   done < ${metacohort_list}
# done < ${phenotype_list}
