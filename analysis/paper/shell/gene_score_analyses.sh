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
  ${rCNV_bucket}/analysis/gene_scoring/data/rCNV.*.gene_abfs.tsv \
  ${rCNV_bucket}/analysis/gene_scoring/data/rCNV2_analysis_d1.rCNV.*.gene_burden.meta_analysis.stats.bed.gz \
  ${rCNV_bucket}/analysis/gene_scoring/data/rCNV2_analysis_d1.rCNV.*.gene_burden.underpowered_genes.bed.gz \
  ${rCNV_bucket}/analysis/gene_scoring/refs/rCNV.gene_scoring.training_gene_blacklist.bed.gz \
  ${rCNV_bucket}/results/gene_scoring/rCNV.gene_scores.tsv.gz \
  ${rCNV_bucket}/analysis/gene_scoring/all_models \
  ${rCNV_bucket}/analysis/gene_scoring/gene_lists \
  ./
mkdir refs/
gsutil -m cp -r \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.constraint_features.bed.gz \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.variation_features.bed.gz \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.genomic_features.eigenfeatures.bed.gz \
  ${rCNV_bucket}/analysis/paper/data/misc/asc_spark_* \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists \
  ${rCNV_bucket}/analysis/analysis_refs/gene_feature_transformations.tsv \
  ${rCNV_bucket}/analysis/paper/data/misc/gene_feature_metadata.tsv \
  ${rCNV_bucket}/analysis/paper/data/large_segments/${prefix}.master_segments.bed.gz \
  refs/


# Plot UpSet comparisons of gold-standard haploinsufficient genes
echo -e "LoF constrained\trefs/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list" > hi.gs.upset.input.tsv
cat \
  refs/gene_lists/DDG2P.hc_lof.genes.list \
  refs/gene_lists/ClinGen.hc_haploinsufficient.genes.list \
| sort | uniq > hc_lof_disease.genes.list
echo -e "Dominant LoF disease\thc_lof_disease.genes.list" >> hi.gs.upset.input.tsv
echo -e "Low-expressor intolerant\trefs/gene_lists/gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.low_expression_invariant.genes.list" >> hi.gs.upset.input.tsv
/opt/rCNV2/analysis/paper/plot/misc/plot_upset.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --cnv-coloring "DEL" \
  hi.gs.upset.input.tsv \
  ${prefix}.gold_standard_genes.haploinsufficient.upset.pdf


# Plot UpSet comparisons of gold-standard haplosufficient genes
echo -e "Mutationally tolerant\trefs/gene_lists/gnomad.v2.1.1.mutation_tolerant.genes.list" > hs.gs.upset.input.tsv
cat \
  refs/gene_lists/HP0000118.HPOdb.genes.list \
  refs/gene_lists/DDG2P*.genes.list \
  refs/gene_lists/ClinGen*.genes.list \
| fgrep -wvf - refs/gene_lists/gencode.v19.canonical.pext_filtered.genes.list \
| sort | uniq > no_disease_assoc.genes.list
echo -e "No disease assoc.\tno_disease_assoc.genes.list" >> hs.gs.upset.input.tsv
echo -e "Low-expressor tolerant\trefs/gene_lists/gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.low_expression_variable.genes.list" >> hs.gs.upset.input.tsv
/opt/rCNV2/analysis/paper/plot/misc/plot_upset.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --cnv-coloring "DEL" \
  hs.gs.upset.input.tsv \
  ${prefix}.gold_standard_genes.haplosufficient.upset.pdf


# Plot UpSet comparisons of gold-standard triplosensitive genes
echo -e "Missense constrained\trefs/gene_lists/gnomad.v2.1.1.mis_constrained.genes.list" > ts.gs.upset.input.tsv
cat \
  refs/gene_lists/DDG2P.hc_gof.genes.list \
  refs/gene_lists/DDG2P.hc_other.genes.list \
  refs/gene_lists/ClinGen.all_triplosensitive.genes.list \
| sort | uniq > hc_notlof_disease.genes.list
echo -e "Dominant GoF disease\thc_lof_disease.genes.list" >> ts.gs.upset.input.tsv
echo -e "High-expressor intolerant\trefs/gene_lists/gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.high_expression_invariant.genes.list" >> ts.gs.upset.input.tsv
/opt/rCNV2/analysis/paper/plot/misc/plot_upset.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --cnv-coloring "DUP" \
  ts.gs.upset.input.tsv \
  ${prefix}.gold_standard_genes.triplosensitive.upset.pdf


# Plot UpSet comparisons of gold-standard triploinsensitive genes
echo -e "Mutationally tolerant\trefs/gene_lists/gnomad.v2.1.1.mutation_tolerant.genes.list" > ti.gs.upset.input.tsv
cat \
  refs/gene_lists/HP0000118.HPOdb.genes.list \
  refs/gene_lists/DDG2P*.genes.list \
  refs/gene_lists/ClinGen*.genes.list \
| fgrep -wvf - refs/gene_lists/gencode.v19.canonical.pext_filtered.genes.list \
| sort | uniq > no_disease_assoc.genes.list
echo -e "No disease assoc.\tno_disease_assoc.genes.list" >> ti.gs.upset.input.tsv
echo -e "High-expressor tolerant\trefs/gene_lists/gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.high_expression_variable.genes.list" >> ti.gs.upset.input.tsv
/opt/rCNV2/analysis/paper/plot/misc/plot_upset.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --cnv-coloring "DUP" \
  ti.gs.upset.input.tsv \
  ${prefix}.gold_standard_genes.triploinsensitive.upset.pdf


# Plot distribution of training effect sizes and BFDPs
for CNV in DEL DUP; do
  # Combine training blacklist
  zcat *.$CNV.gene_burden.underpowered_genes.bed.gz \
    rCNV.gene_scoring.training_gene_blacklist.bed.gz \
  | fgrep -v "#" | cut -f4 | sort -V | uniq \
  > rCNV.$CNV.training_blacklist.genes.list
  # Set CNV-specific parameters
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
  # Plot
  /opt/rCNV2/analysis/paper/plot/gene_scores/plot_ds_model_training_distribs.R \
    --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
    rCNV2_analysis_d1.rCNV.$CNV.gene_burden.meta_analysis.stats.bed.gz \
    rCNV.$CNV.gene_abfs.tsv \
    $true_genes \
    $false_genes \
    rCNV.$CNV.training_blacklist.genes.list \
    $CNV \
    ${prefix}.$CNV
done


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
    evaluation.${CNV}.input.tsv \
    "$true_genes" \
    "$false_genes" \
    ${prefix}.${CNV}
done


# Compare performance of final scores on HI- or TS-only genes (but not both)
fgrep -wvf \
  gene_lists/gold_standard.triplosensitive.genes.list \
  gene_lists/gold_standard.haploinsufficient.genes.list \
> hi_only.genes.list
fgrep -wvf \
  gene_lists/gold_standard.haploinsufficient.genes.list \
  gene_lists/gold_standard.triplosensitive.genes.list \
> ts_only.genes.list
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_phi_vs_pts_trainingsets.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.gene_scores.tsv.gz \
  hi_only.genes.list \
  gene_lists/gold_standard.haplosufficient.genes.list \
  ${prefix}.HI_only
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_phi_vs_pts_trainingsets.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.gene_scores.tsv.gz \
  ts_only.genes.list \
  gene_lists/gold_standard.triploinsensitive.genes.list \
  ${prefix}.TS_only


# Plot simple scatterplot distributions of scores
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_scores_scatter.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.gene_scores.tsv.gz \
  ${prefix}


# Compare distribution of pHI for selected gene sets
echo -e "ClinGen dom. HI\trefs/gene_lists/ClinGen.hc_haploinsufficient.genes.list\tTRUE" > phi_vs_gene_sets.input.tsv
echo -e "DECIPHER dom. HI\trefs/gene_lists/DDG2P.hc_lof.genes.list\tTRUE" >> phi_vs_gene_sets.input.tsv
echo -e "Mouse het. lethal\trefs/gene_lists/mouse_het_lethal.genes.list\tFALSE" >> phi_vs_gene_sets.input.tsv
echo -e "Cell essential\trefs/gene_lists/cell_essential.genes.list\tFALSE" >> phi_vs_gene_sets.input.tsv
# echo -e "Mouse dispensable\trefs/gene_lists/mouse_dispensable.genes.list\tFALSE" >> phi_vs_gene_sets.input.tsv
echo -e "Cell non-essential\trefs/gene_lists/cell_nonessential.genes.list\tFALSE" >> phi_vs_gene_sets.input.tsv
echo -e "Olfactory receptors\trefs/gene_lists/olfactory_receptors.genes.list\tFALSE" >> phi_vs_gene_sets.input.tsv
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_gene_set_enrichments.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --height 1.7 \
  rCNV.gene_scores.tsv.gz \
  phi_vs_gene_sets.input.tsv \
  pHI \
  ${prefix}


# Compare distribution of pTS for selected gene sets
echo -e "ClinGen dom. TS\trefs/gene_lists/ClinGen.all_triplosensitive.genes.list\tTRUE" > pts_vs_gene_sets.input.tsv
echo -e "DECIPHER dom. GoF\trefs/gene_lists/DDG2P.all_gof.genes.list\tTRUE" >> pts_vs_gene_sets.input.tsv
echo -e "Proto-oncogenes\trefs/gene_lists/COSMIC.hc_oncogenes.genes.list\tFALSE" >> pts_vs_gene_sets.input.tsv
echo -e "Olfactory receptors\trefs/gene_lists/olfactory_receptors.genes.list\tFALSE" >> pts_vs_gene_sets.input.tsv
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_gene_set_enrichments.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  --height 1.35 \
  rCNV.gene_scores.tsv.gz \
  pts_vs_gene_sets.input.tsv \
  pTS \
  ${prefix}


# Compare de novo CNVs in ASD vs gene scores
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_asd_denovo_cnv_analysis.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.gene_scores.tsv.gz \
  refs/gencode.v19.canonical.pext_filtered.constraint_features.bed.gz \
  refs/asc_spark_denovo_cnvs.cleaned.b37.annotated.bed.gz \
  refs/asc_spark_child_phenotypes.list \
  ${prefix}


# Compare rates of LoF deletions & CG duplications in gnomAD-SV vs rCNV gene scores
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_gnomad-sv_comparisons.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.gene_scores.tsv.gz \
  refs/gencode.v19.canonical.pext_filtered.variation_features.bed.gz \
  refs/gencode.v19.canonical.pext_filtered.genomic_features.eigenfeatures.bed.gz \
  refs/gene_lists/gnomad.v2.1.1.likely_unconstrained.genes.list \
  ${prefix}


# Compute enrichment of de novo SNVs from ASC & DDD vs rCNV gene scores
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_dnm_enrichments.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.gene_scores.tsv.gz \
  refs/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
  refs/gene_lists/gnomad.v2.1.1.likely_unconstrained.genes.list \
  ${prefix}


# Perform feature regressions between subgroups of genes
athena transform \
  --transformations-tsv refs/gene_feature_transformations.tsv \
  --ignore-columns 4 \
  --bgzip \
  refs/gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz \
  refs/gencode.v19.canonical.pext_filtered.all_features.no_variation.transformed.bed.gz
/opt/rCNV2/analysis/paper/plot/gene_scores/ds_feature_regressions.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.gene_scores.tsv.gz \
  refs/gencode.v19.canonical.pext_filtered.all_features.no_variation.transformed.bed.gz \
  refs/gene_feature_metadata.tsv \
  ${prefix}


# Plot distributions of features between subgroups of genes
mkdir feature_distribs_by_ds_group/
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_ds_feature_distribs.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.gene_scores.tsv.gz \
  refs/gencode.v19.canonical.pext_filtered.all_features.no_variation.transformed.bed.gz \
  refs/gene_feature_metadata.tsv \
  feature_distribs_by_ds_group/${prefix}


# Plot distributions of features between subgroups of genes
/opt/rCNV2/analysis/paper/plot/gene_scores/driver_gene_prediction.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.gene_scores.tsv.gz \
  refs/${prefix}.master_segments.bed.gz \
  ${prefix}


# Copy all plots to final gs:// directory
gsutil -m cp -r \
  ${prefix}.gold_standard_genes.*.upset.pdf \
  ${prefix}.*.gene_scoring_training_distribs.*pdf \
  ${prefix}.*.model_eval*pdf \
  ${prefix}.*pHI_vs_pTS.*pdf \
  ${prefix}.gene_scores_scatterplot*pdf \
  ${prefix}.*.geneset_enrichments.*pdf \
  ${prefix}.asc_spark_denovo_cnvs*pdf \
  ${prefix}.scores_vs_gnomAD-SV*pdf \
  ${prefix}*dnm_enrichments.*pdf \
  ${prefix}.gradient_regression*pdf \
  feature_distribs_by_ds_group \
  ${prefix}.gd_driver_genes* \
  ${rCNV_bucket}/analysis/paper/plots/gene_scores/
