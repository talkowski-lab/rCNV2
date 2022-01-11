#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master code block for secondary analyses of gene dosage sensitivity scores


# Launch docker image & authenticate GCP credentials
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d2"


# Download necessary data (note: requires permissions)
gsutil -m cp -r \
  ${rCNV_bucket}/analysis/gene_scoring/data/rCNV.*.gene_abfs.tsv \
  ${rCNV_bucket}/analysis/gene_scoring/data/${prefix}.rCNV.*.gene_burden.meta_analysis.stats.bed.gz \
  ${rCNV_bucket}/analysis/gene_scoring/data/${prefix}.rCNV.*.gene_burden.underpowered_genes.bed.gz \
  ${rCNV_bucket}/analysis/gene_scoring/refs/rCNV.gene_scoring.training_gene_excludelist.bed.gz \
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
if ! [ -e training_data ]; then
  mkdir training_data
fi
echo -e "LoF constrained\trefs/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list" > hi.gs.upset.input.tsv
cat \
  refs/gene_lists/DDG2P.hc_lof.genes.list \
  refs/gene_lists/ClinGen.hc_haploinsufficient.genes.list \
| sort | uniq > hc_lof_disease.genes.list
echo -e "Dominant LoF disease\thc_lof_disease.genes.list" >> hi.gs.upset.input.tsv
echo -e "Low-expressor intolerant\trefs/gene_lists/gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.low_expression_invariant.genes.list" >> hi.gs.upset.input.tsv
echo -e "No LoF SVs in gnomAD\trefs/gene_lists/gnomad_sv.v2.1.nonneuro.no_lof_dels.genes.list" >> hi.gs.upset.input.tsv
/opt/rCNV2/analysis/paper/plot/misc/plot_upset.R \
  --min-highlight 3 \
  --cnv-coloring "DEL" \
  hi.gs.upset.input.tsv \
  training_data/${prefix}.gold_standard_genes.haploinsufficient.upset.pdf


# Plot UpSet comparisons of gold-standard haplosufficient genes
if ! [ -e training_data ]; then
  mkdir training_data
fi
echo -e "Mutationally tolerant\trefs/gene_lists/gnomad.v2.1.1.mutation_tolerant.genes.list" > hs.gs.upset.input.tsv
cat \
  refs/gene_lists/HP0000118.HPOdb.genes.list \
  refs/gene_lists/DDG2P*.genes.list \
  refs/gene_lists/ClinGen*.genes.list \
| fgrep -wvf - refs/gene_lists/gencode.v19.canonical.pext_filtered.genes.list \
| sort | uniq > no_disease_assoc.genes.list
echo -e "No disease assoc.\tno_disease_assoc.genes.list" >> hs.gs.upset.input.tsv
echo -e "Low-expressor tolerant\trefs/gene_lists/gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.low_expression_variable.genes.list" >> hs.gs.upset.input.tsv
echo -e "Has LoF SVs in gnomAD\trefs/gene_lists/gnomad_sv.v2.1.nonneuro.has_lof_dels.genes.list" >> hs.gs.upset.input.tsv
/opt/rCNV2/analysis/paper/plot/misc/plot_upset.R \
  --min-highlight 4 \
  --cnv-coloring "DEL" \
  hs.gs.upset.input.tsv \
  training_data/${prefix}.gold_standard_genes.haplosufficient.upset.pdf


# Plot UpSet comparisons of gold-standard triplosensitive genes
if ! [ -e training_data ]; then
  mkdir training_data
fi
echo -e "Missense constrained\trefs/gene_lists/gnomad.v2.1.1.mis_constrained.genes.list" > ts.gs.upset.input.tsv
cat \
  refs/gene_lists/DDG2P.hc_gof.genes.list \
  refs/gene_lists/DDG2P.hc_other.genes.list \
  refs/gene_lists/ClinGen.all_triplosensitive.genes.list \
| sort | uniq > hc_notlof_disease.genes.list
echo -e "Dominant GoF disease\thc_lof_disease.genes.list" >> ts.gs.upset.input.tsv
echo -e "High-expressor intolerant\trefs/gene_lists/gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.high_expression_invariant.genes.list" >> ts.gs.upset.input.tsv
echo -e "No CG DUPs in gnomAD\trefs/gene_lists/gnomad_sv.v2.1.nonneuro.no_cg_dups.genes.list" >> ts.gs.upset.input.tsv
/opt/rCNV2/analysis/paper/plot/misc/plot_upset.R \
  --min-highlight 3 \
  --cnv-coloring "DUP" \
  ts.gs.upset.input.tsv \
  training_data/${prefix}.gold_standard_genes.triplosensitive.upset.pdf


# Plot UpSet comparisons of gold-standard triploinsensitive genes
if ! [ -e training_data ]; then
  mkdir training_data
fi
echo -e "Mutationally tolerant\trefs/gene_lists/gnomad.v2.1.1.mutation_tolerant.genes.list" > ti.gs.upset.input.tsv
cat \
  refs/gene_lists/HP0000118.HPOdb.genes.list \
  refs/gene_lists/DDG2P*.genes.list \
  refs/gene_lists/ClinGen*.genes.list \
| fgrep -wvf - refs/gene_lists/gencode.v19.canonical.pext_filtered.genes.list \
| sort | uniq > no_disease_assoc.genes.list
echo -e "No disease assoc.\tno_disease_assoc.genes.list" >> ti.gs.upset.input.tsv
echo -e "High-expressor tolerant\trefs/gene_lists/gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.high_expression_variable.genes.list" >> ti.gs.upset.input.tsv
echo -e "Has CG DUPs in gnomAD\trefs/gene_lists/gnomad_sv.v2.1.nonneuro.has_cg_dups.genes.list" >> ti.gs.upset.input.tsv
/opt/rCNV2/analysis/paper/plot/misc/plot_upset.R \
  --min-highlight 4 \
  --cnv-coloring "DUP" \
  ti.gs.upset.input.tsv \
  training_data/${prefix}.gold_standard_genes.triploinsensitive.upset.pdf


# Plot distribution of training effect sizes and BFDPs
if ! [ -e training_data ]; then
  mkdir training_data
fi
for CNV in DEL DUP; do
  # Combine training excludelist
  zcat *.$CNV.gene_burden.underpowered_genes.bed.gz \
    rCNV.gene_scoring.training_gene_excludelist.bed.gz \
  | fgrep -v "#" | cut -f4 | sort -V | uniq \
  > rCNV.$CNV.training_excludelist.genes.list
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
    ${prefix}.rCNV.$CNV.gene_burden.meta_analysis.stats.bed.gz \
    rCNV.$CNV.gene_abfs.tsv \
    $true_genes \
    $false_genes \
    rCNV.$CNV.training_excludelist.genes.list \
    $CNV \
    training_data/${prefix}.$CNV
done


# Plot performance of various ML models
if ! [ -e model_comparisons ]; then
  mkdir model_comparisons
fi
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
    evaluation.${CNV}.input.tsv \
    "$true_genes" \
    "$false_genes" \
    model_comparisons/${prefix}.${CNV}
done


# Compare performance of final scores on HI- or TS-only genes (but not both)
if ! [ -e model_comparisons ]; then
  mkdir model_comparisons
fi
fgrep -wvf \
  gene_lists/gold_standard.triplosensitive.genes.list \
  gene_lists/gold_standard.haploinsufficient.genes.list \
> hi_only.genes.list
fgrep -wvf \
  gene_lists/gold_standard.haploinsufficient.genes.list \
  gene_lists/gold_standard.triplosensitive.genes.list \
> ts_only.genes.list
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_phaplo_vs_ptriplo_trainingsets.R \
  rCNV.gene_scores.tsv.gz \
  hi_only.genes.list \
  gene_lists/gold_standard.haplosufficient.genes.list \
  model_comparisons/${prefix}.HI_only
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_phaplo_vs_ptriplo_trainingsets.R \
  rCNV.gene_scores.tsv.gz \
  ts_only.genes.list \
  gene_lists/gold_standard.triploinsensitive.genes.list \
  model_comparisons/${prefix}.TS_only


# Plot basic distributions of scores
if ! [ -e basic_distribs ]; then
  mkdir basic_distribs
fi
/opt/rCNV2/analysis/paper/plot/gene_scores/empirical_score_cutoffs.R \
  rCNV.gene_scores.tsv.gz \
  ${prefix}.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz \
  ${prefix}.rCNV.DUP.gene_burden.meta_analysis.stats.bed.gz \
  refs/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
  rCNV.gene_scoring.training_gene_excludelist.bed.gz \
  basic_distribs/${prefix}
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_scores_scatter.R \
  rCNV.gene_scores.tsv.gz \
  basic_distribs/${prefix}
# TODO: remove mean-assigned genes with missing scores from pHaplo/pTriplo score correlations
/opt/rCNV2/analysis/paper/plot/gene_scores/score_vs_score_correlations.R \
  rCNV.gene_scores.tsv.gz \
  refs/gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz \
  refs/gene_feature_metadata.tsv \
  basic_distribs/${prefix}


# Compare distribution of pHaplo for selected gene sets
if ! [ -e enrichments ]; then
  mkdir enrichments
fi
echo -e "ClinGen dom. HI\trefs/gene_lists/ClinGen.hc_haploinsufficient.genes.list\tTRUE" > phaplo_vs_gene_sets.input.tsv
echo -e "DECIPHER dom. HI\trefs/gene_lists/DDG2P.hc_lof.genes.list\tTRUE" >> phaplo_vs_gene_sets.input.tsv
echo -e "Mouse het. lethal\trefs/gene_lists/mouse_het_lethal.genes.list\tFALSE" >> phaplo_vs_gene_sets.input.tsv
echo -e "Cell essential\trefs/gene_lists/cell_essential.genes.list\tFALSE" >> phaplo_vs_gene_sets.input.tsv
echo -e "Cell non-essential\trefs/gene_lists/cell_nonessential.genes.list\tFALSE" >> phaplo_vs_gene_sets.input.tsv
echo -e "Olfactory receptors\trefs/gene_lists/olfactory_receptors.genes.list\tFALSE" >> phaplo_vs_gene_sets.input.tsv
# TODO: ADD DDD+ASC TO THIS ANALYSIS
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_gene_set_enrichments.R \
  --height 1.7 \
  rCNV.gene_scores.tsv.gz \
  phaplo_vs_gene_sets.input.tsv \
  pHaplo \
  enrichments/${prefix}


# Compare distribution of pTriplo for selected gene sets
if ! [ -e enrichments ]; then
  mkdir enrichments
fi
echo -e "ClinGen dom. TS\trefs/gene_lists/ClinGen.all_triplosensitive.genes.list\tTRUE" > ptriplo_vs_gene_sets.input.tsv
echo -e "DECIPHER dom. GoF\trefs/gene_lists/DDG2P.all_gof.genes.list\tTRUE" >> ptriplo_vs_gene_sets.input.tsv
echo -e "Proto-oncogenes\trefs/gene_lists/COSMIC.hc_oncogenes.genes.list\tFALSE" >> ptriplo_vs_gene_sets.input.tsv
echo -e "Olfactory receptors\trefs/gene_lists/olfactory_receptors.genes.list\tFALSE" >> ptriplo_vs_gene_sets.input.tsv
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_gene_set_enrichments.R \
  --height 1.35 \
  rCNV.gene_scores.tsv.gz \
  ptriplo_vs_gene_sets.input.tsv \
  pTriplo \
  enrichments/${prefix}


# Compare de novo CNVs in ASD vs gene scores
if ! [ -e enrichments ]; then
  mkdir enrichments
fi
### TODO: COULD CONSIDER ADDING OTHER NEW SCORES TO THESE COMPARISONS
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_asd_denovo_cnv_analysis.R \
  rCNV.gene_scores.tsv.gz \
  refs/gencode.v19.canonical.pext_filtered.constraint_features.bed.gz \
  refs/asc_spark_denovo_cnvs.cleaned.b37.annotated.bed.gz \
  refs/asc_spark_child_phenotypes.list \
  enrichments/${prefix}


# Compare rates of LoF deletions & CG duplications in gnomAD-SV vs rCNV gene scores
if ! [ -e enrichments ]; then
  mkdir enrichments
fi
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_gnomad-sv_comparisons.R \
  rCNV.gene_scores.tsv.gz \
  refs/gencode.v19.canonical.pext_filtered.variation_features.bed.gz \
  refs/gencode.v19.canonical.pext_filtered.genomic_features.eigenfeatures.bed.gz \
  refs/gene_lists/gnomad.v2.1.1.likely_unconstrained.genes.list \
  enrichments/${prefix}


# Compute enrichment of de novo SNVs from ASC & DDD vs rCNV gene scores
if ! [ -e enrichments ]; then
  mkdir enrichments
fi
### TODO: DEBUG THIS. SOMETHING DOESNT LOOK RIGHT WITH ASC & ASC UNAFFECTED
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_dnm_enrichments.R \
  rCNV.gene_scores.tsv.gz \
  refs/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
  refs/gene_lists/gnomad.v2.1.1.likely_unconstrained.genes.list \
  enrichments/${prefix}


# Transform features prior to feature regressions
# Note: the code chunk below must be run in a docker with athena installed, such as:
# us.gcr.io/broad-dsmap/athena-cloud
# The output from the below athena transform call must be copied back into this docker
docker run --rm -it us.gcr.io/broad-dsmap/athena-cloud
gcloud auth login
export rCNV_bucket="gs://rcnv_project"
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz \
  ${rCNV_bucket}/analysis/analysis_refs/gene_feature_transformations.tsv \
  ./
athena transform \
  --transformations-tsv gene_feature_transformations.tsv \
  --ignore-columns 4 \
  --bgzip \
  gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz \
  gencode.v19.canonical.pext_filtered.all_features.no_variation.transformed.bed.gz


# Perform feature regressions between subgroups of genes
if ! [ -e feature_regressions ]; then
  mkdir feature_regressions
fi
/opt/rCNV2/analysis/paper/plot/gene_scores/ds_feature_regressions.R \
  rCNV.gene_scores.tsv.gz \
  gencode.v19.canonical.pext_filtered.all_features.no_variation.transformed.bed.gz \
  refs/gene_feature_metadata.tsv \
  feature_regressions/${prefix}


# Plot distributions of features between subgroups of genes
if ! [ -e feature_distribs_by_ds_group ]; then
  mkdir feature_distribs_by_ds_group
fi
/opt/rCNV2/analysis/paper/plot/gene_scores/plot_ds_feature_distribs.R \
  rCNV.gene_scores.tsv.gz \
  gencode.v19.canonical.pext_filtered.all_features.no_variation.transformed.bed.gz \
  refs/gene_feature_metadata.tsv \
  feature_distribs_by_ds_group/${prefix}


# Plot predicted driver genes for GD segments
if ! [ -e basic_distribs ]; then
  mkdir basic_distribs
fi
# TODO: FIX SEGMENT LABELS IN PLOTS
/opt/rCNV2/analysis/paper/plot/gene_scores/driver_gene_prediction.R \
  rCNV.gene_scores.tsv.gz \
  refs/${prefix}.master_segments.bed.gz \
  basic_distribs/${prefix}


# Copy all plots to final gs:// directory
gsutil -m cp -r \
  training_data \
  model_comparisons \
  basic_distribs \
  enrichments \
  feature_regressions \
  feature_distribs_by_ds_group \
  ${rCNV_bucket}/analysis/paper/plots/gene_scores/
