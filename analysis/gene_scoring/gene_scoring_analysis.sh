#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Prediction of dosage sensitivity scores for all genes


# Launch docker image
docker run --rm -it talkowski/rcnv
gcloud auth login


# Test/dev parameters
hpo="HP:0000118"
prefix="HP0000118"
freq_code="rCNV"
CNV="DEL"
phenotype_list="test_phenotypes.list"
# metacohort_list="refs/rCNV_metacohort_list.txt"
metacohort_sample_table="refs/HPOs_by_metacohort.table.tsv"
# gtf="genes/gencode.v19.canonical.gtf.gz"
rCNV_bucket="gs://rcnv_project"
theta0_del=1.167
theta0_dup=0.778
theta1=2.373
var=1.467
prior=0.115
gene_features="gencode.v19.canonical.pext_filtered.all_features.eigenfeatures.bed.gz"
raw_gene_features="gencode.v19.canonical.pext_filtered.all_features.bed.gz"
elnet_alpha=0.1
elnet_l1_l2_mix=1


# Copy all gene lists, metadata, and rCNV association stats
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/genes/gene_lists ./
mkdir gene_metadata && \
  gsutil -m cp -r ${rCNV_bucket}/cleaned_data/genes/metadata ./gene_metadata
mkdir stats && \
  gsutil -m cp -r \
    ${rCNV_bucket}/analysis/gene_burden/${prefix}/${freq_code}/stats/**meta_analysis**bed.gz \
    ./stats/


# Define list of high-confidence dosage-sensitive genes
fgrep -wf gene_lists/DDG2P.hc_lof.genes.list \
  gene_lists/ClinGen.hc_haploinsufficient.genes.list \
> gold_standard.ad_disease.genes.list
fgrep -wf gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
  gold_standard.ad_disease.genes.list \
> gold_standard.haploinsufficient.genes.list


# Define list of high-confidence dosage-insensitive genes
cat gene_lists/HP0000118.HPOdb.genes.list \
  gene_lists/DDG2P*.genes.list \
  gene_lists/ClinGen*.genes.list \
| fgrep -wvf - gene_lists/gencode.v19.canonical.pext_filtered.genes.list \
> gold_standard.no_disease_assoc.genes.list
fgrep -wf gene_lists/gnomad.v2.1.1.mutation_tolerant.genes.list \
  gold_standard.no_disease_assoc.genes.list \
> gold_standard.haplosufficient.genes.list


# Score all genes
for CNV in DEL DUP; do
  # Set CNV-specific variables
  case $CNV in
    "DEL")
      theta0=${theta0_del}
      ;;
    "DUP")
      theta0=${theta0_dup}
      ;;
  esac

  # Compute BF & BFDR for all genes
  /opt/rCNV2/analysis/gene_scoring/calc_gene_bfs.py \
    --theta0 $theta0 \
    --theta1 ${theta1} \
    --var0 ${var} \
    --prior ${prior} \
    --outfile ${freq_code}.$CNV.gene_abfs.tsv \
    stats/${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed.gz

  # Score all genes
  /opt/rCNV2/analysis/gene_scoring/score_genes.py \
    --regularization-alpha ${elnet_alpha} \
    --regularization-l1-l2-mix ${elnet_l1_l2_mix} \
    --outfile ${freq_code}.$CNV.gene_scores.tsv \
    ${freq_code}.$CNV.gene_abfs.tsv \
    ${gene_features}
done


# Merge scores
/opt/rCNV2/analysis/gene_scoring/merge_del_dup_scores.R \
  ${freq_code}.DEL.gene_scores.tsv \
  ${freq_code}.DUP.gene_scores.tsv \
  ${freq_code}.gene_scores.tsv
gzip -f ${freq_code}.gene_scores.tsv


#### Quality assessment of gene scores ####

# Plot correlations of raw features vs scores
/opt/rCNV2/analysis/gene_scoring/plot_score_feature_cors.R \
  ${freq_code}.gene_scores.tsv.gz \
  ${raw_gene_features} \
  ${freq_code}.gene_scores.raw_feature_cors

# Make all gene truth sets
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/genes/gene_lists ./

# Make CNV type-dependent truth sets
for CNV in DEL DUP; do
  case $CNV in
    "DEL")
      # Union (ClinGen HI + DDG2P dominant LoF)
      cat gene_lists/ClinGen.hmc_haploinsufficient.genes.list \
          gene_lists/DDG2P.hmc_lof.genes.list \
      | sort -Vk1,1 | uniq \
      > union_truth_set.lof.tsv

      # Intersection (ClinGen HI + DDG2P dominant LoF)
      fgrep -wf \
        gene_lists/ClinGen.hmc_haploinsufficient.genes.list \
        gene_lists/DDG2P.hmc_lof.genes.list \
      | sort -Vk1,1 | uniq \
      > intersection_truth_set.lof.tsv

      # ClinGen HI alone
      cat gene_lists/ClinGen.all_haploinsufficient.genes.list \
      | sort -Vk1,1 | uniq \
      > clingen_truth_set.lof.tsv

      # DDG2P lof + other alone
      cat gene_lists/DDG2P.all_lof.genes.list \
          gene_lists/DDG2P.all_other.genes.list \
      | sort -Vk1,1 | uniq \
      > ddg2p_truth_set.lof.tsv

      # Combine all truth sets into master union
      cat union_truth_set.lof.tsv \
          clingen_truth_set.lof.tsv \
          ddg2p_truth_set.lof.tsv \
      | sort -Vk1,1 | uniq \
      > master_union_truth_set.lof.tsv

      # Write truth set input tsv
      for wrapper in 1; do
        echo -e "Union truth set\tmaster_union_truth_set.lof.tsv\tgrey25"
        echo -e "ClinGen & DECIPHER (union)\tunion_truth_set.lof.tsv\t#9F2B1C"
        echo -e "ClinGen & DECIPHER (int.)\tintersection_truth_set.lof.tsv\t#D43925"
        echo -e "ClinGen dom. HI\tclingen_truth_set.lof.tsv\t#DD6151"
        echo -e "DECIPHER dom. LoF/unk.\tddg2p_truth_set.lof.tsv\t#E5887C"
      done > DEL.roc_truth_sets.tsv
      ;;

    "DUP")
      # Union (ClinGen HI + DDG2P dominant CG)
      cat gene_lists/ClinGen.hmc_triplosensitive.genes.list \
          gene_lists/DDG2P.hmc_gof.genes.list \
      | sort -Vk1,1 | uniq \
      > union_truth_set.gof.tsv

      # ClinGen triplo alone
      cat gene_lists/ClinGen.all_triplosensitive.genes.list \
      | sort -Vk1,1 | uniq \
      > clingen_truth_set.triplo.tsv

      # DDG2P gof + other alone
      cat gene_lists/DDG2P.all_gof.genes.list \
          gene_lists/DDG2P.all_other.genes.list \
      | sort -Vk1,1 | uniq \
      > ddg2p_truth_set.gof.tsv

      # Combine all truth sets into master union
      cat union_truth_set.gof.tsv \
          clingen_truth_set.triplo.tsv \
          ddg2p_truth_set.gof.tsv \
      | sort -Vk1,1 | uniq \
      > master_union_truth_set.gof.tsv

      # Write truth set input tsv
      for wrapper in 1; do
        echo -e "Union truth set\tmaster_union_truth_set.gof.tsv\tgrey25"
        echo -e "ClinGen & DECIPHER GoF (union)\tunion_truth_set.gof.tsv\t#1A5985"
        echo -e "ClinGen dom. TS\tclingen_truth_set.triplo.tsv\t#2376B2"
        echo -e "DECIPHER dom. GoF/unk.\tddg2p_truth_set.gof.tsv\t#4F91C1"
      done > DUP.roc_truth_sets.tsv
      ;;
  esac    
done

# Plot ROC, PRC, and enrichments
mkdir ${freq_code}_gene_scoring_QC_plots/
/opt/rCNV2/analysis/gene_scoring/plot_gene_score_qc.R \
  ${freq_code}.gene_scores.tsv.gz \
  DEL.finemap_roc_truth_sets.tsv \
  DUP.finemap_roc_truth_sets.tsv \
  ${freq_code}_gene_scoring_QC_plots/${freq_code}_gene_score_qc

