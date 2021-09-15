#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2019-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Code to convert sample-level phenotype keywords to standardized HPO hierarchy


# Launch docker image
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv


# Add helper alias for formatting long integers
alias addcom="sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta'"


# Copy all raw phenotype data from the project Google Bucket (note: requires permissions)
gcloud auth login
mkdir raw_phenos/
gsutil cp -r gs://rcnv_project/raw_data/phenotypes/* raw_phenos/
mkdir cleaned_phenos/
mkdir cleaned_phenos/all/
mkdir cleaned_phenos/intermediate/
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


# Pool old and new GeneDx phenotypes into single file
cat raw_phenos/GDX_2013.raw_phenos.txt \
    raw_phenos/GDX_2021.raw_phenos.txt \
| sort -Vk1,1 -k2,2V \
> raw_phenos/GDX.raw_phenos.txt


# Convert BCH, GDX, Coe, SSC, and IU phenotypes (no post-hoc QC necessary)
for cohort in BCH GDX Coe SSC IU; do
  echo -e "Starting $cohort"
  time /opt/rCNV2/data_curation/phenotype/indication_to_HPO.py \
    --obo hp.obo \
    -s /opt/rCNV2/refs/hpo/supplementary_hpo_mappings.tsv \
    -x /opt/rCNV2/refs/hpo/break_hpo_mappings.tsv \
    --no-match-default "HP:0000001;HP:0000118;UNKNOWN" \
    -o cleaned_phenos/all/${cohort}.cleaned_phenos.txt \
    raw_phenos/${cohort}.raw_phenos.txt
done
# Add dummy control lines to Coe phenotype file
yes "HEALTHY_CONTROL" \
| head -n $( fgrep Coe /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f4 ) \
| awk -v OFS="\t" '{ print "Coe_CONTROL_"NR, $1 }' \
>> cleaned_phenos/all/Coe.cleaned_phenos.txt


# Convert SickKids and Epi25k (post-hoc QC required)
for cohort in SickKids Epi25k CHOP; do
  echo -e "Starting $cohort"
  time /opt/rCNV2/data_curation/phenotype/indication_to_HPO.py \
    --obo hp.obo \
    -s /opt/rCNV2/refs/hpo/supplementary_hpo_mappings.tsv \
    -x /opt/rCNV2/refs/hpo/break_hpo_mappings.tsv \
    --no-match-default "HP:0000001;HP:0000118;UNKNOWN" \
    -o cleaned_phenos/all/${cohort}.cleaned_phenos.preQC.txt \
    raw_phenos/${cohort}.raw_phenos.txt
  fgrep -wf raw_phenos/${cohort}.QC_pass_samples.list \
    cleaned_phenos/all/${cohort}.cleaned_phenos.preQC.txt \
  > cleaned_phenos/all/${cohort}.cleaned_phenos.txt
done


# Prep conversion table for cohorts with uniform phenotypes (TSAICG, PGC)
cut -f2 raw_phenos/CHOP.raw_phenos.txt | sort | uniq \
| awk -v OFS="\t" '{ if ($1=="") $1="UNKNOWN"; print $1, $1 }' \
> CHOP.raw_phenos.conversion_input.txt
time /opt/rCNV2/data_curation/phenotype/indication_to_HPO.py \
  --obo hp.obo \
  -s /opt/rCNV2/refs/hpo/supplementary_hpo_mappings.tsv \
  -x /opt/rCNV2/refs/hpo/break_hpo_mappings.tsv \
  --no-match-default "HP:0000001;HP:0000118;UNKNOWN" \
  -o CHOP.raw_phenos.conversion_table.txt \
  CHOP.raw_phenos.conversion_input.txt


# Convert CHOP phenotypes
awk -v FS="\t" -v OFS="\t" '{ if ($2=="") $2="UNKNOWN"; print $1, $2 }' \
  raw_phenos/CHOP.raw_phenos.txt \
| sort -Vk2,2 \
| join -1 2 -2 1 -t $'\t' - <( sort -Vk1,1 CHOP.raw_phenos.conversion_table.txt ) \
| cut -f2-3 | sort -Vk1,1 \
> cleaned_phenos/all/CHOP.cleaned_phenos.preQC.txt
fgrep -wf raw_phenos/CHOP.QC_pass_samples.list \
  cleaned_phenos/all/CHOP.cleaned_phenos.preQC.txt \
> cleaned_phenos/all/CHOP.cleaned_phenos.txt


# Make dummy phenotype files for TSAICG
yes $( fgrep "tourette" CHOP.raw_phenos.conversion_table.txt | cut -f2 ) \
| head -n $( fgrep TSAICG /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f3 ) \
| awk -v OFS="\t" '{ print "TSAICG_CASE_"NR, $1 }' \
> cleaned_phenos/all/TSAICG.cleaned_phenos.txt
yes "HEALTHY_CONTROL" \
| head -n $( fgrep TSAICG /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f4 ) \
| awk -v OFS="\t" '{ print "TSAICG_CONTROL_"NR, $1 }' \
>> cleaned_phenos/all/TSAICG.cleaned_phenos.txt


# Make dummy phenotype files for PGC
yes $( fgrep "schizophrenia" CHOP.raw_phenos.conversion_table.txt | cut -f2 ) \
| head -n $( fgrep PGC /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f3 ) \
| awk -v OFS="\t" '{ print "PGC_CASE_"NR, $1 }' \
> cleaned_phenos/all/PGC.cleaned_phenos.txt
yes "HEALTHY_CONTROL" \
| head -n $( fgrep PGC /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f4 ) \
| awk -v OFS="\t" '{ print "PGC_CONTROL_"NR, $1 }' \
>> cleaned_phenos/all/PGC.cleaned_phenos.txt


# Make dummy phenotype files for control-only cohorts (Cooper, GSD, Ontario, TCGA)
for cohort in Cooper GSD Ontario TCGA; do
  yes "HEALTHY_CONTROL" \
  | head -n $( fgrep ${cohort} /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f4 ) \
  | awk -v OFS="\t" -v cohort=${cohort} '{ print cohort"_CONTROL_"NR, $1 }' \
  > cleaned_phenos/all/${cohort}.cleaned_phenos.txt
done


# Make dummy phenotype file for RUMC ID cases
echo -e "dummy\tintellectual_disability" > id.dummy.tsv
time /opt/rCNV2/data_curation/phenotype/indication_to_HPO.py \
  --obo hp.obo \
  -s /opt/rCNV2/refs/hpo/supplementary_hpo_mappings.tsv \
  -x /opt/rCNV2/refs/hpo/break_hpo_mappings.tsv \
  --no-match-default "HP:0000001;HP:0000118;UNKNOWN" \
  -o id.conversion_table.txt \
  id.dummy.tsv
yes $( fgrep "dummy" id.conversion_table.txt | cut -f2 ) \
| head -n $( fgrep RUMC /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f3 ) \
| awk -v OFS="\t" '{ print "RUMC_CASE_"NR, $1 }' \
> cleaned_phenos/all/RUMC.cleaned_phenos.txt



# Tabulate UKBB ICD-10 dictionary of all terms with at least 48 samples for manual review
# 48 samples = 1:1000 phenotype incidence in UKBB
/opt/rCNV2/data_curation/phenotype/tabulate_ukbb_icd10_dictionary.py \
  --minsamp 48 \
  --outfile UKBB_ICD10_wSampleCounts.for_review.txt \
  UKBB_ICD10_manifest.tsv \
  raw_phenos/UKBB.sample_IDs_w_ICD10.txt


# Convert UKBB ICD-10 codes to indications
time /opt/rCNV2/data_curation/phenotype/icd10_to_indication.py \
  --whitelist-terms /opt/rCNV2/refs/icd10/UKBB_ICD_term_whitelist.txt \
  --blacklist-terms /opt/rCNV2/refs/icd10/UKBB_ICD_term_blacklist.txt \
  --blacklist-samples /opt/rCNV2/refs/icd10/UKBB_ICD_sample_blacklist.txt \
  --blacklisted-samples-outfile UKBB.phenotype_blacklisted_samples.txt \
  --default healthy \
  --report-fails \
  --outfile raw_phenos/UKBB.raw_phenos.txt \
  UKBB_ICD10_manifest.tsv \
  raw_phenos/UKBB.sample_IDs_w_ICD10.txt


# Copy UKBB sample IDs that failed phenotype QC to Google bucket (note: requires permissions)
gsutil cp UKBB.phenotype_blacklisted_samples.txt \
  gs://rcnv_project/raw_data/phenotypes/


# Convert UKBB indications
cut -f2 raw_phenos/UKBB.raw_phenos.txt | sort | uniq \
| awk -v OFS="\t" '{ print $1, $1 }' \
> UKBB.raw_phenos.conversion_input.txt
time /opt/rCNV2/data_curation/phenotype/indication_to_HPO.py \
  --obo hp.obo \
  -s /opt/rCNV2/refs/hpo/supplementary_hpo_mappings.tsv \
  -x /opt/rCNV2/refs/hpo/break_hpo_mappings.tsv \
  --no-match-default "HP:0000001;HP:0000118;UNKNOWN" \
  -o UKBB.raw_phenos.conversion_table.txt \
  UKBB.raw_phenos.conversion_input.txt
sort -Vk2,2 raw_phenos/UKBB.raw_phenos.txt \
| join -1 2 -2 1 -t $'\t' - <( sort -Vk1,1 UKBB.raw_phenos.conversion_table.txt ) \
| cut -f2-3 | sort -Vk1,1 \
> cleaned_phenos/all/UKBB.cleaned_phenos.preQC.txt
fgrep -wf raw_phenos/UKBB.QC_pass_samples.list \
  cleaned_phenos/all/UKBB.cleaned_phenos.preQC.txt \
> cleaned_phenos/all/UKBB.cleaned_phenos.txt


# Pool all phenotypes across cohorts
while read cohort; do
  pheno_file="cleaned_phenos/all/${cohort}.cleaned_phenos.txt"
  if [ -e $pheno_file ]; then
    cat $pheno_file
  fi
done < <( cut -f1 /opt/rCNV2/refs/rCNV_sample_counts.txt | fgrep -v "#" ) \
> all_phenos.merged.txt


# Make list of counts per HPO pair for EstBB & BioVU
echo -e "EstBB\traw_phenos/EstBB.HPO_terms_full_cooccurrence_table.tsv.gz" > hpo_pair_cohorts.inputs.tsv
echo -e "BioVU\traw_phenos/BioVU.hpo_coocurrence_table.tsv.gz" >> hpo_pair_cohorts.inputs.tsv


# Determine minimum HPO tree to use
/opt/rCNV2/data_curation/phenotype/collapse_HPO_tree.py \
  --hpo-pair-cohorts hpo_pair_cohorts.inputs.tsv \
  --ignore "HP:0000001" \
  --ignore "HP:0031796" \
  --ignore "HP:0031797" \
  --ignore "HP:0011008" \
  --ignore "HP:0025303" \
  --ignore "HP:0025142" \
  --ignore "HP:0012823" \
  --obo hp.obo \
  --raw-counts samples_per_HPO.txt \
  --filter-log HPO_tree_filter.log \
  --min-samples 3000 \
  --min-diff 3000 \
  --max-jac 0.8 \
  --outfile phenotype_groups.HPO_metadata.intermediate.txt \
  all_phenos.merged.txt


# Restrict all cohort phenotypes to terms in minimal HPO tree
while read cohort; do
  if [ -e cleaned_phenos/all/${cohort}.cleaned_phenos.txt ]; then
    /opt/rCNV2/data_curation/phenotype/filter_HPO_per_sample.py \
      -o cleaned_phenos/intermediate/${cohort}.cleaned_phenos.txt \
      cleaned_phenos/all/${cohort}.cleaned_phenos.txt \
      phenotype_groups.HPO_metadata.intermediate.txt
  fi
done < <( cut -f1 /opt/rCNV2/refs/rCNV_sample_counts.txt )
for cohort in UKBB CHOP Epi25k; do
  /opt/rCNV2/data_curation/phenotype/filter_HPO_per_sample.py \
    -o cleaned_phenos/intermediate/${cohort}.cleaned_phenos.preQC.txt \
    cleaned_phenos/all/${cohort}.cleaned_phenos.preQC.txt \
    phenotype_groups.HPO_metadata.intermediate.txt
done


# Get summary table of HPO counts per cohort & metacohort
# Only keep HPO terms with at least 300 cases from three or more metacohorts
gsutil cp gs://rcnv_project/analysis/analysis_refs/rCNV_metacohort_list.txt ./
/opt/rCNV2/data_curation/phenotype/gather_hpo_per_cohort_table.py \
  --outfile HPOs_by_cohort.table.tsv \
  --meta-cohorts rCNV_metacohort_list.txt \
  --meta-out HPOs_by_metacohort.table.tsv \
  --min-metacohorts 3 \
  --min-per-metacohort 300 \
  --hpo-pair-cohorts hpo_pair_cohorts.inputs.tsv \
  --pairwise-curation /opt/rCNV2/refs/hpo/manual_hpo_pair_similarity_drops.tsv \
  phenotype_groups.HPO_metadata.intermediate.txt \
  /opt/rCNV2/refs/rCNV_sample_counts.txt \
  cleaned_phenos/intermediate/


# Clean up output from original HPO tree consolidation to reflect metacohort filtering
fgrep -v "#" HPOs_by_cohort.table.tsv \
| cut -f1 \
| sort -Vk1,1 \
> final_HPOs.txt
while read hpo; do
  awk -v hpo=${hpo} '{ if ($1==hpo) print $0 }' \
  phenotype_groups.HPO_metadata.intermediate.txt
done < final_HPOs.txt \
> phenotype_groups.HPO_metadata.tmp
while IFS=$'\t' read hpo descrip n tier oparents ochildren; do
  parents=$( echo "${oparents}" | sed -e 's/;/\n/g' \
             | fgrep -wf final_HPOs.txt | paste -s -d\; )
  if [ -z ${parents} ] || [ ${parents} == "" ]; then
    parents="NA"
  fi
  children=$( echo "${ochildren}" | sed -e 's/;/\n/g' \
             | fgrep -wf final_HPOs.txt | paste -s -d\; )
  if [ -z ${children} ] || [ ${children} == "" ]; then
    children="NA"
  fi
  echo -e "${hpo}\t${descrip}\t${n}\t${tier}\t${parents}\t${children}"
done < phenotype_groups.HPO_metadata.tmp \
| sort -t$'\t' -nrk3,3 \
| cat <( grep -e '^#' phenotype_groups.HPO_metadata.intermediate.txt ) - \
> phenotype_groups.HPO_metadata.txt
rm phenotype_groups.HPO_metadata.tmp


# Print HTML table of HPO metadata for README
for wrapper in 1; do
  echo -e "| HPO Term | Description | Samples | HPO Tier | Parent Terms | Child Terms |  "
  echo -e "| :--- | :--- | ---: | ---: | :--- | :--- |  "
  paste \
    <( fgrep -v "#" phenotype_groups.HPO_metadata.txt | cut -f1-2 ) \
    <( fgrep -v "#" phenotype_groups.HPO_metadata.txt | cut -f3 | addcom ) \
    <( fgrep -v "#" phenotype_groups.HPO_metadata.txt | cut -f4- ) \
  | sed -e 's/\t/\ \|\ /g' -e 's/^/\|\ /g' -e 's/$/\ \|\ \ /g' -e 's/\;/,\ /g' 
done


# Barplot of samples per metacohort per HPO
/opt/rCNV2/data_curation/phenotype/plot_hpo_per_cohort.R \
  HPOs_by_metacohort.table.tsv \
  HPOs_by_metacohort.barplot.jpg
gsutil cp HPOs_by_metacohort.barplot.jpg \
  gs://rcnv_project/public/
gsutil acl ch -u AllUsers:R gs://rcnv_project/public/HPOs_by_metacohort.barplot.jpg


# Make simple file mapping HPO codes to directory prefixes (no colons)
paste <( cut -f1 phenotype_groups.HPO_metadata.txt | sed 's/\://g' ) \
      <( cut -f1 phenotype_groups.HPO_metadata.txt ) \
| fgrep -v "#" \
| fgrep -v "HEALTHY_CONTROL" \
> test_phenotypes.list


# Restrict all cohort phenotypes to final analysis terms
while read cohort; do
  if [ -e cleaned_phenos/all/${cohort}.cleaned_phenos.txt ]; then
    /opt/rCNV2/data_curation/phenotype/filter_HPO_per_sample.py \
      -o cleaned_phenos/filtered/${cohort}.cleaned_phenos.txt \
      cleaned_phenos/all/${cohort}.cleaned_phenos.txt \
      phenotype_groups.HPO_metadata.txt
  fi
done < <( cut -f1 /opt/rCNV2/refs/rCNV_sample_counts.txt )
for cohort in UKBB CHOP; do
  /opt/rCNV2/data_curation/phenotype/filter_HPO_per_sample.py \
    -o cleaned_phenos/filtered/${cohort}.cleaned_phenos.preQC.txt \
    cleaned_phenos/all/${cohort}.cleaned_phenos.preQC.txt \
    phenotype_groups.HPO_metadata.txt
done


# Filter HPO coocurrence tables for EstBB & BioVU
/opt/rCNV2/data_curation/phenotype/filter_coocurrence_table.py \
  raw_phenos/EstBB.HPO_terms_full_cooccurrence_table.tsv.gz \
  <( fgrep -v "#" phenotype_groups.HPO_metadata.txt | cut -f1 ) \
  cleaned_phenos/filtered/EstBB.final_cooccurrence_table.tsv
gzip -f cleaned_phenos/filtered/EstBB.final_cooccurrence_table.tsv
/opt/rCNV2/data_curation/phenotype/filter_coocurrence_table.py \
  raw_phenos/BioVU.hpo_coocurrence_table.tsv.gz \
  <( fgrep -v "#" phenotype_groups.HPO_metadata.txt | cut -f1 ) \
  cleaned_phenos/filtered/BioVU.final_cooccurrence_table.tsv
gzip -f cleaned_phenos/filtered/BioVU.final_cooccurrence_table.tsv


# Merge final phenotype lists per metacohort
while read name cohorts; do
  for cohort in $( echo "$cohorts" | sed 's/;/\n/g' ); do
    pheno_file="cleaned_phenos/filtered/${cohort}.cleaned_phenos.txt"
    if [ -e $pheno_file ]; then
      cat $pheno_file
    fi
  done \
  > cleaned_phenos/filtered/${name}.cleaned_phenos.txt
done < rCNV_metacohort_list.txt


# Copy all final data to Google bucket (requires permissions)
gsutil -m cp -r cleaned_phenos/* gs://rcnv_project/cleaned_data/phenotypes/
gsutil cp \
  samples_per_HPO.txt \
  HPO_tree_filter.log \
  phenotype_groups.HPO_metadata.txt \
  gs://rcnv_project/cleaned_data/phenotypes/hpo_logs_metadata/
gsutil cp test_phenotypes.list \
  gs://rcnv_project/analysis/analysis_refs/


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
gsutil cp \
  HPOs_by_cohort.table.tsv \
  HPOs_by_metacohort.table.tsv \
  rCNV_metacohort_sample_counts.txt \
  gs://rcnv_project/analysis/analysis_refs/


# Print HTML tables of HPO counts per cohort & metacohort
/opt/rCNV2/data_curation/phenotype/gather_hpo_per_cohort_table.py \
  --outfile HPOs_by_cohort.table.html.tsv \
  --meta-cohorts rCNV_metacohort_list.txt \
  --meta-out HPOs_by_metacohort.table.html.tsv \
  --hpo-pair-cohorts hpo_pair_cohorts.inputs.tsv \
  --html \
  phenotype_groups.HPO_metadata.txt \
  /opt/rCNV2/refs/rCNV_sample_counts.txt \
  cleaned_phenos/filtered/

