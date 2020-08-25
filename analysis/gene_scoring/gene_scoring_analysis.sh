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


# Set parameters
export rCNV_bucket="gs://rcnv_project"


# Copy all filtered CNV data, gene coordinates, and other references 
# from the project Google Bucket (note: requires permissions)
mkdir cleaned_cnv/ phenos/
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/cnv/* cleaned_cnv/
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/phenotypes/filtered/meta*.cleaned_phenos.txt phenos/
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/genes ./
mkdir refs/
gsutil -m cp ${rCNV_bucket}/refs/GRCh37.*.bed.gz refs/
gsutil -m cp ${rCNV_bucket}/analysis/analysis_refs/* refs/
gsutil -m cp ${rCNV_bucket}/analysis/paper/data/large_segments/lit_GDs.*.bed.gz refs/


# Test/dev parameters
prefix="HP0000118"
freq_code="rCNV"
CNV="DEL"
meta="meta1"
phenotype_list="refs/test_phenotypes.list"
metacohort_list="refs/rCNV_metacohort_list.txt"
metacohort_sample_table="refs/HPOs_by_metacohort.table.tsv"
gtf="genes/gencode.v19.canonical.gtf.gz"
gtf_index="genes/gencode.v19.canonical.gtf.gz.tbi"
contig=18
pad_controls=0
max_cnv_size=300000000
weight_mode="weak"
min_cds_ovr_del=0.8
min_cds_ovr_dup=0.8
max_genes_per_cnv=24
meta_model_prefix="fe"
min_cnvs_per_gene_training=3
cen_tel_dist=1000000
prior_frac=0.12


# Create training blacklist: Remove all genes within ±1Mb of a telomere/centromere, 
# or those within known genomic disorder regions
zcat refs/GRCh37.centromeres_telomeres.bed.gz \
| awk -v FS="\t" -v OFS="\t" -v d=${cen_tel_dist} \
  '{ if ($2-d<0) print $1, "0", $3+d; else print $1, $2-d, $3+d }' \
| cat - <( zcat refs/lit_GDs.*.bed.gz | cut -f1-3 | fgrep -v "#" ) \
| sort -Vk1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
> ${freq_code}.gene_scoring.training_mask.bed.gz
gsutil -m cp \
  ${rCNV_bucket}/analysis/gene_burden/HP0000118/rCNV/stats/HP0000118.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz \
  ./
bedtools intersect -u \
  -a HP0000118.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz \
  -b ${freq_code}.gene_scoring.training_mask.bed.gz \
| cut -f1-4 \
| cat <( echo -e "#chr\tstart\tend\tgene" ) - \
| bgzip -c \
> ${freq_code}.gene_scoring.training_gene_blacklist.bed.gz

# Copy training blacklist to project bucket (note: requires permissions)
gsutil -m cp \
  ${freq_code}.gene_scoring.training_gene_blacklist.bed.gz \
  ${rCNV_bucket}/analysis/gene_scoring/refs/


# Create training gene lists
# Define list of high-confidence haploinsufficient genes
cat \
  genes/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
  <( cat genes/gene_lists/DDG2P.hc_lof.genes.list genes/gene_lists/ClinGen.hc_haploinsufficient.genes.list | sort | uniq ) \
  genes/gene_lists/gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.low_expression_invariant.genes.list \
| sort | uniq -c \
| awk '{ if ($1>=2) print $2 }' \
| fgrep -wf - genes/gene_lists/gencode.v19.canonical.pext_filtered.genes.list \
> gold_standard.haploinsufficient.genes.list

# Define list of high-confidence haplosufficient genes
cat \
  genes/gene_lists/gnomad.v2.1.1.mutation_tolerant.genes.list \
  <( cat genes/gene_lists/HP0000118.HPOdb.genes.list genes/gene_lists/DDG2P*.genes.list genes/gene_lists/ClinGen*.genes.list \
     | fgrep -wvf - genes/gene_lists/gencode.v19.canonical.pext_filtered.genes.list | sort | uniq ) \
  genes/gene_lists/gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.low_expression_variable.genes.list \
| sort | uniq -c \
| awk '{ if ($1>=2) print $2 }' \
| fgrep -wf - genes/gene_lists/gencode.v19.canonical.pext_filtered.genes.list \
> gold_standard.haplosufficient.genes.list

# Define list of high-confidence triplosensitive genes
cat \
  genes/gene_lists/gnomad.v2.1.1.mis_constrained.genes.list \
  <( cat genes/gene_lists/DDG2P.hc_gof.genes.list genes/gene_lists/DDG2P.hc_other.genes.list genes/gene_lists/ClinGen.all_triplosensitive.genes.list | sort | uniq ) \
  genes/gene_lists/gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.high_expression_invariant.genes.list \
| sort | uniq -c \
| awk '{ if ($1>=2) print $2 }' \
| fgrep -wf - genes/gene_lists/gencode.v19.canonical.pext_filtered.genes.list \
> gold_standard.triplosensitive.genes.list

# Define list of high-confidence haplosufficient genes
cat \
  genes/gene_lists/gnomad.v2.1.1.mutation_tolerant.genes.list \
  <( cat genes/gene_lists/HP0000118.HPOdb.genes.list genes/gene_lists/DDG2P*.genes.list genes/gene_lists/ClinGen*.genes.list \
     | fgrep -wvf - genes/gene_lists/gencode.v19.canonical.pext_filtered.genes.list | sort | uniq ) \
  genes/gene_lists/gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.high_expression_variable.genes.list \
| sort | uniq -c \
| awk '{ if ($1>=2) print $2 }' \
| fgrep -wf - genes/gene_lists/gencode.v19.canonical.pext_filtered.genes.list \
> gold_standard.triploinsensitive.genes.list

# Copy all gold-standard gene lists to output bucket (note: requires permissions)
gsutil -m cp \
  gold_standard.*.genes.list \
  ${rCNV_bucket}/analysis/gene_scoring/gene_lists/


# Note: data for parameter optimization is computed in parallel in FireCloud with
# gene_scoring_analysis.wdl, and must be run before the below code:
gsutil -m cp -r ${rCNV_bucket}/analysis/gene_scoring/optimization_data ./

# Postprocess optimization data
# Note: the --max-genes-for-summary parameter must be determined first by running the
# subsequent code block (number of genes per CNV), then running this block a second time
while read meta cohorts; do
  echo $meta
  # Process deletions
  echo DEL
  /opt/rCNV2/analysis/gene_scoring/distill_optimization_data.py \
    --genes-per-cnv optimization_data/${freq_code}.DEL.$meta.genes_per_cnv.tsv.gz \
    --hpos <( cut -f2 refs/test_phenotypes.list ) \
    --positive-truth-genes gold_standard.haploinsufficient.genes.list \
    --negative-truth-genes gold_standard.haplosufficient.genes.list \
    --exclude-genes <( zcat ${freq_code}.gene_scoring.training_gene_blacklist.bed.gz \
                       | fgrep -v "#" | cut -f4 ) \
    --max-genes-for-summary 24 \
    --summary-counts ${meta}.optimization_data.counts_per_hpo.DEL.tsv \
    --cnv-stats ${meta}.optimization_data.counts_per_cnv.DEL.tsv \
    --gzip
  # Process duplications
  echo DUP
  /opt/rCNV2/analysis/gene_scoring/distill_optimization_data.py \
    --genes-per-cnv optimization_data/${freq_code}.DUP.$meta.genes_per_cnv.tsv.gz \
    --hpos <( cut -f2 refs/test_phenotypes.list ) \
    --positive-truth-genes gold_standard.triplosensitive.genes.list \
    --negative-truth-genes gold_standard.triploinsensitive.genes.list \
    --exclude-genes <( zcat ${freq_code}.gene_scoring.training_gene_blacklist.bed.gz \
                       | fgrep -v "#" | cut -f4 ) \
    --max-genes-for-summary 24 \
    --summary-counts ${meta}.optimization_data.counts_per_hpo.DUP.tsv \
    --cnv-stats ${meta}.optimization_data.counts_per_cnv.DUP.tsv \
    --gzip
done < <( fgrep -v mega refs/rCNV_metacohort_list.txt )

# Parameter optimization: number of genes per CNV
for CNV in DEL DUP; do
  while read meta cohorts; do
    echo -e "${meta}\t${meta}.optimization_data.counts_per_cnv.${CNV}.tsv.gz"
  done < <( fgrep -v mega refs/rCNV_metacohort_list.txt ) \
  > optimize_genes_per_cnv.${CNV}.input.tsv
done
/opt/rCNV2/analysis/gene_scoring/optimize_genes_per_cnv.R \
  optimize_genes_per_cnv.DEL.input.tsv \
  optimize_genes_per_cnv.DUP.input.tsv \
  optimize_genes_per_cnv.results.pdf

# Parameter optimization: phenotype selection
for CNV in DEL DUP; do
  while read meta cohorts; do
    echo -e "${meta}\t${meta}.optimization_data.counts_per_hpo.${CNV}.tsv.gz"
  done < <( fgrep -v mega refs/rCNV_metacohort_list.txt ) \
  > select_hpos.${CNV}.input.tsv
done
/opt/rCNV2/analysis/gene_scoring/select_hpos.R \
  select_hpos.DEL.input.tsv \
  select_hpos.DUP.input.tsv \
  gene_scoring.hpos_to_keep.list \
  select_hpos.results.pdf
while read meta cohorts; do
  /opt/rCNV2/analysis/gene_scoring/determine_effective_case_n.py \
    phenos/$meta.cleaned_phenos.txt \
    gene_scoring.hpos_to_keep.list \
    | paste <( echo $meta ) -
done < <( fgrep -v mega refs/rCNV_metacohort_list.txt ) \
> gene_scoring.effective_case_sample_sizes.tsv
gsutil -m cp \
  gene_scoring.hpos_to_keep.list \
  gene_scoring.effective_case_sample_sizes.tsv \
  ${rCNV_bucket}/analysis/gene_scoring/refs/


# Parameter optimization: minimum number of CNVs per gene for training
# NOTE: requires running meta-analysis and get_underpowered_genes tasks first in WDL
gsutil -m cp ${rCNV_bucket}/analysis/gene_scoring/data/** ./
for CNV in DEL DUP; do
  /opt/rCNV2/analysis/gene_scoring/variance_vs_counts.R \
    rCNV2_analysis_d1.rCNV.$CNV.gene_burden.meta_analysis.stats.bed.gz \
    rCNV2_analysis_d1.rCNV.$CNV.counts_per_gene.tsv \
    variance_vs_cnvs.$CNV.pdf
done


# Recompute association stats per cohort
export training_hpo_list=gene_scoring.hpos_to_keep.list
export effective_case_sample_sizes=gene_scoring.effective_case_sample_sizes.tsv
for contig in $( seq 1 22 ); do

  # Extract contig of interest from GTF
  tabix ${gtf} ${contig} | bgzip -c > ${contig}.gtf.gz

  # Create string of HPOs to keep
  keep_hpos=$( cat ${training_hpo_list} | paste -s -d\; )

  # Iterate over metacohorts to compute single-cohort stats
  while read meta cohorts; do
    echo $meta

    # Set metacohort-specific parameters
    cnv_bed="cleaned_cnv/$meta.${freq_code}.bed.gz"
    effective_case_n=$( fgrep -w $meta ${effective_case_sample_sizes} | cut -f2 )

    # Iterate over CNV types
    for CNV in DEL DUP; do
      echo $CNV

      # Set CNV-specific parameters
      case "$CNV" in
        DEL)
          min_cds_ovr=${min_cds_ovr_del}
          ;;
        DUP)
          min_cds_ovr=${min_cds_ovr_dup}
          ;;
      esac

      # Count CNVs
      /opt/rCNV2/analysis/genes/count_cnvs_per_gene.py \
        --max-cnv-size ${max_cnv_size} \
        --weight-mode ${weight_mode} \
        --min-cds-ovr $min_cds_ovr \
        --max-genes ${max_genes_per_cnv} \
        -t $CNV \
        --hpo "$keep_hpos" \
        --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
        --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
        --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
        -z \
        --verbose \
        -o "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.${contig}.bed.gz" \
        "$cnv_bed" \
        ${contig}.gtf.gz
      tabix -f "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.${contig}.bed.gz"

      # Perform burden test
      /opt/rCNV2/analysis/genes/gene_burden_test.R \
        --pheno-table ${metacohort_sample_table} \
        --cohort-name $meta \
        --cnv $CNV \
        --effective-case-n $effective_case_n \
        --bgzip \
        "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.${contig}.bed.gz" \
        "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.${contig}.bed.gz"
      tabix -f "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.${contig}.bed.gz"
    done
  done < <( fgrep -v "mega" ${metacohort_list} )

done

# Recompute meta-analysis statistics
for CNV in DEL DUP; do
  echo $CNV

  # Make list of stats files to be considered
  find / -name "*.${prefix}.${freq_code}.${CNV}.gene_burden.stats.*.bed.gz" \
  > stats.paths.list

  # Merge burden stats per cohort
  zcat $( sed -n '1p' stats.paths.list ) | sed -n '1p' > header.tsv
  while read meta cohort; do
    zcat $( fgrep $meta stats.paths.list ) \
    | grep -ve '^#' \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | cat header.tsv - \
    | bgzip -c \
    > "$meta.${prefix}.${freq_code}.${CNV}.gene_burden.stats.bed.gz"
    tabix -f "$meta.${prefix}.${freq_code}.${CNV}.gene_burden.stats.bed.gz"
  done < <( fgrep -v mega ${metacohort_list} )

  # Make input for meta-analysis
  while read meta cohorts; do
    echo -e "$meta\t$meta.${prefix}.${freq_code}.${CNV}.gene_burden.stats.bed.gz"
  done < <( fgrep -v mega ${metacohort_list} ) \
  > ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.input.txt
  
  # Run meta-analysis
  /opt/rCNV2/analysis/genes/gene_meta_analysis.R \
    --model ${meta_model_prefix} \
    --p-is-phred \
    --spa \
    ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.input.txt \
    ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed
  bgzip -f ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed
  tabix -p bed -f ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz

  # Extract list of genes with < min_cnvs_per_gene_training (to be used as blacklist later for training)
  # Note: this cutoff must be pre-determined with eval_or_vs_cnv_counts.R, but is not automated here
  /opt/rCNV2/analysis/gene_scoring/get_underpowered_genes.R \
    --min-cnvs ${min_cnvs_per_gene_training} \
    --gene-counts-out ${prefix}.${freq_code}.${CNV}.counts_per_gene.tsv \
    ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.input.txt \
    ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed
  awk -v OFS="\t" '{ print $1, $2, $3, $4 }' \
    ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed \
  | bgzip -c \
  > ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed.gz
  tabix -p bed -f ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed.gz
done


# Localize meta-analysis data (only necessary for local development)
gsutil -m cp \
  ${rCNV_bucket}/analysis/gene_scoring/data/**.gene_burden.meta_analysis.stats.bed.gz \
  ${rCNV_bucket}/analysis/gene_scoring/data/*.gene_burden.underpowered_genes.bed.gz \
  ./
del_meta_stats="rCNV2_analysis_d1.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz"
dup_meta_stats="rCNV2_analysis_d1.rCNV.DUP.gene_burden.meta_analysis.stats.bed.gz"


# Create CNV-type-specific gene blacklists
for CNV in DEL DUP; do
  zcat *.$CNV.gene_burden.underpowered_genes.bed.gz \
    ${freq_code}.gene_scoring.training_gene_blacklist.bed.gz \
  | fgrep -v "#" | cut -f4 | sort -V | uniq \
  > ${freq_code}.$CNV.training_blacklist.genes.list
done

# Compute prior effect sizes
prior_lnor_thresholding_pct=1
/opt/rCNV2/analysis/gene_scoring/estimate_prior_effect_sizes.R \
  --pct ${prior_lnor_thresholding_pct} \
  ${del_meta_stats} \
  ${dup_meta_stats} \
  ${freq_code}.DEL.training_blacklist.genes.list \
  ${freq_code}.DUP.training_blacklist.genes.list \
  gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
  gold_standard.haploinsufficient.genes.list \
  gold_standard.haplosufficient.genes.list \
  gold_standard.triplosensitive.genes.list \
  gold_standard.triploinsensitive.genes.list \
  ${freq_code}.prior_estimation
awk -v FS="\t" '{ if ($1=="theta0" && $2=="DEL") print $3 }' \
  ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
> theta0_del.tsv
awk -v FS="\t" '{ if ($1=="theta0" && $2=="DUP") print $3 }' \
  ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
> theta0_dup.tsv
awk -v FS="\t" '{ if ($1=="theta1" && $2=="DEL") print $3 }' \
  ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
> theta1_del.tsv
awk -v FS="\t" '{ if ($1=="theta1" && $2=="DUP") print $3 }' \
  ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
> theta1_dup.tsv
awk -v FS="\t" '{ if ($1=="var0" && $2=="DEL") print $3 }' \
  ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
> var0_del.tsv
awk -v FS="\t" '{ if ($1=="var0" && $2=="DUP") print $3 }' \
  ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
> var0_dup.tsv
awk -v FS="\t" '{ if ($1=="var1" && $2=="DEL") print $3 }' \
  ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
> var1_del.tsv
awk -v FS="\t" '{ if ($1=="var1" && $2=="DUP") print $3 }' \
  ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
> var1_dup.tsv

# Compute BFDP per gene
for CNV in DEL DUP; do
  # Set CNV-specific variables
  case $CNV in
    "DEL")
      theta0=$( cat theta0_del.tsv )
      theta1=$( cat theta1_del.tsv )
      statsbed=${del_meta_stats}
      ;;
    "DUP")
      theta0=$( cat theta0_dup.tsv )
      theta1=$( cat theta1_dup.tsv )
      statsbed=${dup_meta_stats}
      ;;
  esac

  # Compute BF & BFDR for all genes
  /opt/rCNV2/analysis/gene_scoring/calc_gene_bfs.py \
    --theta0 $theta0 \
    --theta1 $theta1 \
    --var0 1 \
    --prior ${prior_frac} \
    --blacklist ${freq_code}.$CNV.training_blacklist.genes.list \
    --outfile ${freq_code}.$CNV.gene_abfs.tsv \
    "$statsbed"
done




# Score genes for a single model & CNV type
# Dev/test parameters
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.no_variation.* \
  ./
freq_code="rCNV"
rCNV_bucket="${rCNV_bucket}"
CNV="DEL"
BFDP_stats="${freq_code}.${CNV}.gene_abfs.tsv"
blacklist="${freq_code}.gene_scoring.training_gene_blacklist.bed.gz"
underpowered_genes="rCNV2_analysis_d1.rCNV.DEL.gene_burden.underpowered_genes.bed.gz"
gene_features="gencode.v19.canonical.pext_filtered.all_features.no_variation.eigenfeatures.bed.gz"
raw_gene_features="gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz"
max_true_bfdp=0.5
min_false_bfdp=0.5
model="logit"
elnet_alpha=0.1
elnet_l1_l2_mix=1

# Copy necessary references
gsutil -m cp \
  ${rCNV_bucket}/refs/GRCh37.centromeres_telomeres.bed.gz \
  ${rCNV_bucket}/analysis/gene_scoring/gene_lists/*genes.list \
  ./

# Merge blacklist and list of underpowered genes
zcat ${blacklist} ${underpowered_genes} \
| grep -ve '^#' \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
| uniq \
| cat <( zcat ${blacklist} | grep -e '^#' ) - \
| bgzip -c  \
> blacklist_plus_underpowered.bed.gz

# Set CNV type-specific parameters
case ${CNV} in
  "DEL")
    true_pos="gold_standard.haploinsufficient.genes.list"
    true_neg="gold_standard.haplosufficient.genes.list"
    ;;
  "DUP")
    true_pos="gold_standard.triplosensitive.genes.list"
    true_neg="gold_standard.triploinsensitive.genes.list"
    ;;
esac

# Score genes
/opt/rCNV2/analysis/gene_scoring/score_genes.py \
  --centromeres GRCh37.centromeres_telomeres.bed.gz \
  --blacklist blacklist_plus_underpowered.bed.gz \
  --true-positives $true_pos \
  --true-negatives $true_neg \
  --model ${model} \
  --max-true-bfdp ${max_true_bfdp} \
  --min-false-bfdp ${min_false_bfdp} \
  --regularization-alpha ${elnet_alpha} \
  --regularization-l1-l2-mix ${elnet_l1_l2_mix} \
  --outfile ${freq_code}.${CNV}.gene_scores.${model}.tsv \
  ${BFDP_stats} \
  ${gene_features}


# Ensemble classifier - requires running all individual models first
for CNV in DEL DUP; do
  # Make inputs for ensemble classifier
  find ./ -name "${freq_code}.${CNV}.gene_scores.*.tsv" \
  > ensemble_input.tsv
  case ${CNV} in
    DEL)
      pos_genes="gold_standard.haploinsufficient.genes.list"
      neg_genes="gold_standard.haplosufficient.genes.list"
      ;;
    DUP)
      pos_genes="gold_standard.triplosensitive.genes.list"
      neg_genes="gold_standard.triploinsensitive.genes.list"
      ;;
  esac

  # Run ensemble classifier
  /opt/rCNV2/analysis/gene_scoring/ensemble_classifier.R \
    ensemble_input.tsv \
    $pos_genes \
    $neg_genes \
    ${freq_code}.${CNV}.gene_scores.ensemble.tsv
done


# Copy gene lists
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/genes/gene_lists ./
gsutil -m cp ${rCNV_bucket}/analysis/gene_scoring/gene_lists/* ./gene_lists/

# Evaluate every model to determine best overall predictor
compdir=${freq_code}_gene_scoring_model_comparisons
if ! [ -e $compdir ]; then
  mkdir $compdir
fi
for CNV in DEL DUP; do
  for wrapper in 1; do
    for model in logit svm randomforest lda naivebayes sgd neuralnet ensemble; do
      echo $model
      echo ${freq_code}.$CNV.gene_scores.$model.tsv
    done | paste - -
  done > ${freq_code}.$CNV.model_evaluation.input.tsv
done
/opt/rCNV2/analysis/gene_scoring/compare_models.R \
  ${freq_code}.DEL.model_evaluation.input.tsv \
  ${freq_code}.DUP.model_evaluation.input.tsv \
  gene_lists/gold_standard.haploinsufficient.genes.list \
  gene_lists/gold_standard.triplosensitive.genes.list \
  gene_lists/gold_standard.haplosufficient.genes.list \
  gene_lists/gold_standard.triploinsensitive.genes.list \
  $compdir/${freq_code}_gene_scoring_model_comparisons

# Merge scores from best model (highest harmonic mean AUC)
sed -n '2p' $compdir/${freq_code}_gene_scoring_model_comparisons.summary_table.tsv \
| cut -f1 > best_model.tsv
best_model=$( cat best_model.tsv )
/opt/rCNV2/analysis/gene_scoring/merge_del_dup_scores.R \
  ${freq_code}.DEL.gene_scores.$best_model.tsv \
  ${freq_code}.DUP.gene_scores.$best_model.tsv \
  ${freq_code}.gene_scores.tsv
gzip -f ${freq_code}.gene_scores.tsv

# Plot correlations of raw features vs scores
mkdir gene_score_corplots
/opt/rCNV2/analysis/gene_scoring/plot_score_feature_cors.R \
  ${freq_code}.gene_scores.tsv.gz \
  ${raw_gene_features} \
  gene_score_corplots/${freq_code}.gene_scores.raw_feature_cors

# Make CNV type-dependent truth sets
for CNV in DEL DUP; do
  case $CNV in
    "DEL")
      cat gene_lists/ClinGen.hc_haploinsufficient.genes.list \
          gene_lists/DDG2P.hc_lof.genes.list \
          gene_lists/cell_essential.genes.list \
          gene_lists/mouse_het_lethal.genes.list \
      | sort -Vk1,1 | uniq \
      > master_union_truth_set.lof.tsv

      # Write truth set input tsv
      for wrapper in 1; do
        echo -e "Union truth set\tmaster_union_truth_set.lof.tsv\tgrey25"
        echo -e "ClinGen dom. HI\tgene_lists/ClinGen.hc_haploinsufficient.genes.list\t#9F2B1C"
        echo -e "DECIPHER dom. LoF\tgene_lists/DDG2P.hc_lof.genes.list\t#D43925"
        echo -e "Cell essential\tgene_lists/cell_essential.genes.list\t#DD6151"
        echo -e "Mouse het. lethal\tgene_lists/mouse_het_lethal.genes.list\t#E5887C"
      done > DEL.roc_truth_sets.tsv
      ;;

    "DUP")
      cat gene_lists/ClinGen.all_triplosensitive.genes.list \
          gene_lists/DDG2P.hc_gof.genes.list \
          gene_lists/COSMIC.hc_oncogenes.genes.list \
      | sort -Vk1,1 | uniq \
      > master_union_truth_set.gof.tsv

      # Write truth set input tsv
      for wrapper in 1; do
        echo -e "Union truth set\tmaster_union_truth_set.gof.tsv\tgrey25"
        echo -e "ClinGen dom. TS\tgene_lists/ClinGen.all_triplosensitive.genes.list\t#1A5985"
        echo -e "DECIPHER dom. GoF\tgene_lists/DDG2P.hc_gof.genes.list\t#2376B2"
        echo -e "COSMIC dom. oncogenes\tgene_lists/COSMIC.hc_oncogenes.genes.list\t#4F91C1"
      done > DUP.roc_truth_sets.tsv
      ;;
  esac    
done

# Plot ROC, PRC, and enrichments
mkdir ${freq_code}_gene_scoring_QC_plots/
/opt/rCNV2/analysis/gene_scoring/plot_gene_score_qc.R \
  ${freq_code}.gene_scores.tsv.gz \
  DEL.roc_truth_sets.tsv \
  DUP.roc_truth_sets.tsv \
  gene_lists/gold_standard.haplosufficient.genes.list \
  gene_lists/gold_standard.triploinsensitive.genes.list \
  ${freq_code}_gene_scoring_QC_plots/${freq_code}_gene_score_qc

# Compute enrichments versus population variation data
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.variation_features.bed.gz \
  ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.genomic_features.eigenfeatures.bed.gz \
  ./
/opt/rCNV2/analysis/gene_scoring/scores_vs_gnomad-SV.R \
  rCNV.gene_scores.tsv.gz \
  gencode.v19.canonical.pext_filtered.variation_features.bed.gz \
  gencode.v19.canonical.pext_filtered.genomic_features.eigenfeatures.bed.gz \
  rCNV_scores_vs_gnomad_sv.pdf

