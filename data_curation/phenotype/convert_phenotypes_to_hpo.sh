#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Code to convert sample-level phenotype keywords to standardized HPO hierarchy


# Launch docker image
docker run --rm -it talkowski/rcnv


# Add helper alias for formatting long integers
alias addcom="sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta'"


# Copy all raw phenotype data from the project Google Bucket (note: requires permissions)
gcloud auth login
mkdir raw_phenos/
gsutil cp -r gs://rcnv_project/raw_data/phenotypes/* raw_phenos/
mkdir cleaned_phenos/
mkdir cleaned_phenos/all/
mkdir cleaned_phenos/filtered/


# Download current version of HPO and ICD10
wget http://purl.obolibrary.org/obo/hp.obo
gsutil cp gs://rcnv_project/refs/UKBB_ICD10_manifest.tsv ./


# Make tsv of HPO obo file (for convenience)
/opt/rCNV2/data_curation/phenotype/HPO_obo_to_tsv.py \
  --obo hp.obo \
  --outfile HPO_dict.tsv
gzip -f HPO_dict.tsv
gsutil cp HPO_dict.tsv.gz gs://rcnv_project/refs/


# Convert BCH, GDX, Coe, and SSC phenotypes
for cohort in BCH GDX Coe SSC; do
  /opt/rCNV2/data_curation/phenotype/indication_to_HPO.py \
    --obo hp.obo \
    -s /opt/rCNV2/refs/hpo/supplementary_hpo_mappings.tsv \
    -x /opt/rCNV2/refs/hpo/break_hpo_mappings.tsv \
    --no-match-default "HP:0000001;HP:0000118;UNKNOWN" \
    -o cleaned_phenos/all/${cohort}.cleaned_phenos.txt \
    raw_phenos/${cohort}.raw_phenos.txt
done


# Convert CHOP phenotypes
/opt/rCNV2/data_curation/phenotype/indication_to_HPO.py \
  --obo hp.obo \
  -s /opt/rCNV2/refs/hpo/supplementary_hpo_mappings.tsv \
  -x /opt/rCNV2/refs/hpo/break_hpo_mappings.tsv \
  --no-match-default "HP:0000001;HP:0000118;UNKNOWN" \
  -o cleaned_phenos/all/CHOP.cleaned_phenos.preQC.txt \
  raw_phenos/CHOP.raw_phenos.txt
cut -f1 raw_phenos/CHOP.QC_pass_samples.list \
| fgrep -wf - cleaned_phenos/all/CHOP.cleaned_phenos.preQC.txt \
> cleaned_phenos/all/CHOP.cleaned_phenos.txt


# Prep conversion table for cohorts with uniform phenotypes (TSAICG, PGC)
cut -f2 raw_phenos/CHOP.raw_phenos.txt \
| sort \
| uniq \
| awk -v OFS="\t" '{ print $1, $1 }' \
> CHOP.raw_phenos.conversion_input.txt
/opt/rCNV2/data_curation/phenotype/indication_to_HPO.py \
  --obo hp.obo \
  -s /opt/rCNV2/refs/hpo/supplementary_hpo_mappings.tsv \
  -x /opt/rCNV2/refs/hpo/break_hpo_mappings.tsv \
  --no-match-default "HP:0000001;HP:0000118;UNKNOWN" \
  -o CHOP.raw_phenos.conversion_table.txt \
  CHOP.raw_phenos.conversion_input.txt


# Make dummy phenotype files for TSAICG
yes $( fgrep "tourette" CHOP.raw_phenos.conversion_table.txt | cut -f2 ) \
| head -n $( fgrep TSAICG /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f3 ) \
| awk -v OFS="\t" '{ print "TSAICG_CASE_"NR, $1 }' \
> cleaned_phenos/all/TSAICG.cleaned_phenos.txt


# Make dummy phenotype files for PGC
yes $( fgrep "schizophrenia" CHOP.raw_phenos.conversion_table.txt | cut -f2 ) \
| head -n $( fgrep PGC /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f3 ) \
| awk -v OFS="\t" '{ print "PGC_CASE_"NR, $1 }' \
> cleaned_phenos/all/PGC.cleaned_phenos.txt


# Tabulate UKBB ICD-10 dictionary of all terms with at least 48 samples for manual review
# 48 samples = 1:1000 phenotype incidence in UKBB
/opt/rCNV2/data_curation/phenotype/tabulate_ukbb_icd10_dictionary.py \
  --minsamp 48 \
  --outfile UKBB_ICD10_wSampleCounts.for_review.txt \
  UKBB_ICD10_manifest.tsv \
  raw_phenos/UKBB.sample_IDs_w_ICD10.txt


# Convert UKBB ICD-10 codes to indications
/opt/rCNV2/data_curation/phenotype/icd10_to_indication.py \
  --whitelist-terms /opt/rCNV2/refs/icd10/UKBB_ICD_term_whitelist.txt \
  --blacklist-terms /opt/rCNV2/refs/icd10/UKBB_ICD_term_blacklist.txt \
  --blacklist-samples /opt/rCNV2/refs/icd10/UKBB_ICD_sample_blacklist.txt \
  --blacklisted-samples-outfile UKBB.phenotype_blacklisted_samples.txt \
  --default healthy \
  --report-fails \
  --outfile raw_phenos/UKBB.raw_phenos.txt \
  UKBB_ICD10_manifest.tsv \
  raw_phenos/UKBB.sample_IDs_w_ICD10.txt


# Copy list of UKBB samples with phenotype QC failures to Google bucket (note: requires permissions)
gsutil cp UKBB.phenotype_blacklisted_samples.txt \
  gs://rcnv_project/


# Convert UKBB indications
/opt/rCNV2/data_curation/phenotype/indication_to_HPO.py \
  --obo hp.obo \
  -s /opt/rCNV2/refs/hpo/supplementary_hpo_mappings.tsv \
  -x /opt/rCNV2/refs/hpo/break_hpo_mappings.tsv \
  --no-match-default "HP:0000001;HP:0000118;UNKNOWN" \
  -o cleaned_phenos/all/UKBB.cleaned_phenos.preQC.txt \
  raw_phenos/UKBB.raw_phenos.txt


# Pool all phenotypes across cohorts
while read cohort; do
  if [ -e cleaned_phenos/all/${cohort}.cleaned_phenos.txt ]; then
    if [ ${cohort} != "CHOP" ]; then
      cat cleaned_phenos/all/${cohort}.cleaned_phenos.txt
    else
      fgrep -wf raw_phenos/CHOP.QC_pass_samples.list \
        cleaned_phenos/all/${cohort}.cleaned_phenos.txt
    fi
  fi
done < <( cut -f1 /opt/rCNV2/refs/rCNV_sample_counts.txt ) \
> all_phenos.merged.txt


# Determine minimum HPO tree to use
/opt/rCNV2/data_curation/phenotype/collapse_HPO_tree.py \
  --ignore "HP:0000001" \
  --obo hp.obo \
  --raw-counts samples_per_HPO.txt \
  --filter-log HPO_tree_filter.log \
  --outfile phenotype_groups.HPO_metadata.txt \
  all_phenos.merged.txt


# Make simple file mapping HPO codes to directory prefixes (no colons)
paste <( cut -f1 phenotype_groups.HPO_metadata.txt | sed 's/\://g' ) \
      <( cut -f1 phenotype_groups.HPO_metadata.txt ) \
| fgrep -v "#" \
| fgrep -v "HEALTHY_CONTROL" \
> test_phenotypes.list


# Restrict all cohort phenotypes to terms in minimal HPO tree
while read cohort; do
  if [ -e cleaned_phenos/all/${cohort}.cleaned_phenos.txt ]; then
    /opt/rCNV2/data_curation/phenotype/filter_HPO_per_sample.py \
      -o cleaned_phenos/filtered/${cohort}.cleaned_phenos.txt \
      cleaned_phenos/all/${cohort}.cleaned_phenos.txt \
      phenotype_groups.HPO_metadata.txt
  fi
done < <( cut -f1 /opt/rCNV2/refs/rCNV_sample_counts.txt )
/opt/rCNV2/data_curation/phenotype/filter_HPO_per_sample.py \
  -o cleaned_phenos/filtered/CHOP.cleaned_phenos.preQC.txt \
  cleaned_phenos/all/CHOP.cleaned_phenos.preQC.txt \
  phenotype_groups.HPO_metadata.txt


# Copy all final data to Google bucket (requires permissions)
gsutil cp -r cleaned_phenos/* gs://rcnv_project/cleaned_data/phenotypes/
gsutil cp samples_per_HPO.txt \
  gs://rcnv_project/cleaned_data/phenotypes/hpo_logs_metadata/
gsutil cp HPO_tree_filter.log \
  gs://rcnv_project/cleaned_data/phenotypes/hpo_logs_metadata/
gsutil cp phenotype_groups.HPO_metadata.txt \
  gs://rcnv_project/cleaned_data/phenotypes/hpo_logs_metadata/
gsutil cp test_phenotypes.list \
  gs://rcnv_project/analysis/analysis_refs/


# Print HTML table of HPO metadata for README
for wrapper in 1; do
  echo -e "| HPO Term | Description | Samples | HPO Tier | Parent Terms | Child Terms |  "
  echo -e "| :--- | :--- | ---: | ---: | :--- | :--- |  "
  paste \
    <( fgrep -v "#" phenotype_groups.HPO_metadata.txt | cut -f1-2 ) \
    <( fgrep -v "#" phenotype_groups.HPO_metadata.txt | cut -f3 | addcom ) \
    <( fgrep -v "#" phenotype_groups.HPO_metadata.txt | cut -f4- ) \
  | sed -e 's/\t/\ \|\ /g' -e 's/^/\|\ /g' -e 's/$/\ \|\ \ /g' -e 's/\;/,\ /g' \
  | fgrep -v "HEALTHY_CONTROL"
done


# Get summary table of HPO counts per cohort & metacohort
gsutil cp gs://rcnv_project/analysis/analysis_refs/rCNV_metacohort_list.txt ./
/opt/rCNV2/data_curation/phenotype/gather_hpo_per_cohort_table.py \
  --outfile HPOs_by_cohort.table.tsv \
  --meta-cohorts rCNV_metacohort_list.txt \
  --meta-out HPOs_by_metacohort.table.tsv \
  phenotype_groups.HPO_metadata.txt \
  /opt/rCNV2/refs/rCNV_sample_counts.txt \
  cleaned_phenos/filtered/


# Get simplified table of metacohort combined case & control counts
echo -e "#metacohort\tsamples\tCASE\tCTRL" \
> rCNV_metacohort_sample_counts.txt
while read name cohorts; do
  col=$( head -n1 HPOs_by_metacohort.table.tsv \
         | sed 's/\t/\n/g' \
         | awk -v name=${name} '{ if ($1==name) print NR }' )
  ncase=$( fgrep -wv HEALTHY_CONTROL HPOs_by_metacohort.table.tsv \
           | fgrep -v "#" \
           | awk -v FS="\t" -v col=${col} '{ print $col }' \
           | sort -nrk1,1 \
           | head -n1 )
  nctrl=$( fgrep -w HEALTHY_CONTROL HPOs_by_metacohort.table.tsv \
           | awk -v FS="\t" -v col=${col} '{ print $col }' )
  ntotal=$(( ${ncase} + ${nctrl} ))
  echo -e "${name}\t${ntotal}\t${ncase}\t${nctrl}"
done < rCNV_metacohort_list.txt \
>> rCNV_metacohort_sample_counts.txt


# Copy sample counts per HPO term per cohort to Google bucket (requires permissions)
gsutil cp HPOs_by_cohort.table.tsv \
  gs://rcnv_project/analysis/analysis_refs/
gsutil cp HPOs_by_metacohort.table.tsv \
  gs://rcnv_project/analysis/analysis_refs/
gsutil cp rCNV_metacohort_sample_counts.txt \
  gs://rcnv_project/analysis/analysis_refs/


# Print HTML tables of HPO counts per cohort & metacohort
/opt/rCNV2/data_curation/phenotype/gather_hpo_per_cohort_table.py \
  --outfile HPOs_by_cohort.table.html.tsv \
  --meta-cohorts rCNV_metacohort_list.txt \
  --meta-out HPOs_by_metacohort.table.html.tsv \
  --html \
  phenotype_groups.HPO_metadata.txt \
  /opt/rCNV2/refs/rCNV_sample_counts.txt \
  cleaned_phenos/filtered/

