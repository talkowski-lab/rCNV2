#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master code block for secondary analyses of gene dosage sensitivity scores


# Launch docker image & authenticate GCP credentials
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d1"


# Download necessary data (note: requires permissions)
gsutil -m cp -r \
  ${rCNV_bucket}/results/gene_scoring/rCNV.gene_scores.tsv.gz \
  ${rCNV_bucket}/analysis/gene_scoring/all_models \
  ${rCNV_bucket}/analysis/gene_scoring/gene_lists \
  ./
mkdir refs/
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz \
  refs/


# Plot performance of various ML models
ls -l all_models/rCNV.DEL*tsv | awk -v FS="." '{ print $(NF-1) }' | sort | uniq \
> models.list
for CNV in DEL DUP; do
  while read model; do
    echo -e "${model}\tall_models/rCNV.${CNV}.gene_scores.${model}.tsv"
  done < models.list \
  > evaluation.${CNV}.input.tsv
  case $CNV in
    DEL)
      true_genes=gene_lists/gold_standard.haploinsufficient.genes.list
      false_genes=gene_lists/gold_standard.haplosufficient.genes.list
      ;;
    DUP)
      true_genes=gene_lists/gold_standard.triplosensitive.genes.list
      false_genes=gene_lists/gold_standard.triploinsensitive.genes.list
      ;;
  esac
  /opt/rCNV2/analysis/paper/plot/gene_scores/plot_ml_performance.R \
    --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
    scores.input.tsv \
    "$true_genes" \
    "$false_genes" \
    ${prefix}.${CNV}
done


# Plot simple scatterplot distributions of scores
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_scores_scatter.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.gene_scores.tsv.gz \
  ${prefix}


# Perform feature regressions between subgroups of genes




# Plot distributions of features between subgroups of genes


# Copy all plots to final gs:// directory
gsutil -m cp \
  ${prefix}.model_comparison*pdf \
  ${prefix}.gene_scores_scatterplot*pdf \
  ${rCNV_bucket}/analysis/paper/plots/gene_scores/
