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


# Copy all raw phenotype data from the project Google Bucket (note: requires permissions)
gcloud auth login
mkdir raw_phenos/
gsutil cp -r gs://rcnv_project/raw_data/phenotypes/* raw_phenos/
mkdir cleaned_phenos/


# Convert BCH & GDX phenotypes
/opt/rCNV2/data_curation/phenotype/indication_to_HPO.py \
  --obo http://purl.obolibrary.org/obo/hp.obo \
  -s /opt/rCNV2/refs/hpo/supplementary_hpo_mappings.tsv \
  -x /opt/rCNV2/refs/hpo/break_hpo_mappings.tsv \
  --no-match-default "UNKNOWN" \
  -o cleaned_phenos/BCH_GDX.cleaned_phenos.txt \
  raw_phenos/BCH_GDX.raw_phenos.txt

# Convert Coe phenotypes
/opt/rCNV2/data_curation/phenotype/indication_to_HPO.py \
  --obo http://purl.obolibrary.org/obo/hp.obo \
  -s /opt/rCNV2/refs/hpo/supplementary_hpo_mappings.tsv \
  -x /opt/rCNV2/refs/hpo/break_hpo_mappings.tsv \
  --no-match-default "UNKNOWN" \
  -o cleaned_phenos/Coe.cleaned_phenos.txt \
  raw_phenos/Coe.raw_phenos.txt

# Convert SSC phenotypes
/opt/rCNV2/data_curation/phenotype/indication_to_HPO.py \
  --obo http://purl.obolibrary.org/obo/hp.obo \
  -s /opt/rCNV2/refs/hpo/supplementary_hpo_mappings.tsv \
  -x /opt/rCNV2/refs/hpo/break_hpo_mappings.tsv \
  --no-match-default "UNKNOWN" \
  -o cleaned_phenos/SSC.cleaned_phenos.txt \
  raw_phenos/SSC.raw_phenos.txt

# Prep conversion table for cohorts with uniform phenotypes (CHOP, TSAICG, PGC)
cut -f2 raw_phenos/CHOP.raw_phenos.txt \
| sort \
| uniq \
| awk -v OFS="\t" '{ print $'

/opt/rCNV2/data_curation/phenotype/indication_to_HPO.py \
  --obo http://purl.obolibrary.org/obo/hp.obo \
  -s /opt/rCNV2/refs/hpo/supplementary_hpo_mappings.tsv \
  -x /opt/rCNV2/refs/hpo/break_hpo_mappings.tsv \
  -e /opt/rCNV2/refs/hpo/eligible_hpo_terms.txt \
  --no-match-default "UNKNOWN" \
  -o cleaned_phenos/CHOP.cleaned_phenos.txt \
  raw_phenos/CHOP.raw_phenos.txt








time /opt/rCNV2/data_curation/phenotype/indication_to_HPO.py \
  --obo http://purl.obolibrary.org/obo/hp.obo \
  -s /opt/rCNV2/refs/hpo/supplementary_hpo_mappings.tsv \
  -x /opt/rCNV2/refs/hpo/break_hpo_mappings.tsv \
  -e /opt/rCNV2/refs/hpo/eligible_hpo_terms.txt \
  --no-match-default "UNKNOWN" \
  -o test_phenotypes_reclassified.txt \
  test_phenotypes.txt

time /opt/rCNV2/data_curation/phenotype/indication_to_HPO.py \
  --obo http://purl.obolibrary.org/obo/hp.obo \
  -s /opt/rCNV2/refs/hpo/supplementary_hpo_mappings.tsv \
  -x /opt/rCNV2/refs/hpo/break_hpo_mappings.tsv \
  --no-match-default "UNKNOWN" \
  -o test_phenotypes_reclassified.txt \
  test_phenotypes.txt

