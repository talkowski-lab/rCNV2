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
mkdir cleaned_cnv/
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/cnv/* cleaned_cnv/
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/genes ./
mkdir refs/
gsutil -m cp ${rCNV_bucket}/refs/GRCh37.*.bed.gz refs/
gsutil -m cp ${rCNV_bucket}/analysis/analysis_refs/* refs/
gsutil -m cp ${rCNV_bucket}/analysis/paper/data/large_segments/lit_GDs.*.bed.gz refs/


# Test/dev parameters
hpo="HP:0000118"
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
max_cnv_size=5000000
weight_mode="weak"
min_cds_ovr_del=1.0
min_cds_ovr_dup=1.0
max_genes_per_cnv=10
meta_model_prefix="fe"
min_cnvs_per_gene_training=5


for contig in $( seq 1 22 ); do

  # Extract contig of interest from GTF
  tabix ${gtf} ${contig} | bgzip -c > ${contig}.gtf.gz

  # Iterate over metacohorts to compute single-cohort stats
  while read meta cohorts; do
    echo $meta

    # Set metacohort-specific parameters
    cnv_bed="cleaned_cnv/$meta.${freq_code}.bed.gz"

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
        --hpo ${hpo} \
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
        --case-hpo ${hpo} \
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
  /opt/rCNV2/analysis/gene_scoring/get_underpowered_genes.R \
    --min-cnvs ${min_cnvs_per_gene_training} \
    ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.input.txt \
    ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed
  bgzip -f ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed
  tabix -p bed -f ${prefix}.${freq_code}.${CNV}.gene_burden.underpowered_genes.bed.gz
done


# The following code chunk performs three tasks in serial:
# 1. Determine genes to blacklist during training
# 2. Estimate effect size priors
# 3. Compute BFDP per gene
freq_code="rCNV"
rCNV_bucket="${rCNV_bucket}"
del_meta_stats="HP0000118.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz"
dup_meta_stats="HP0000118.rCNV.DUP.gene_burden.meta_analysis.stats.bed.gz"
cen_tel_dist=1000000
prior_frac=0.115


# Create blacklist: Remove all genes within ±1Mb of a telomere/centromere, 
# or those within known genomic disorder regions
gsutil -m cp ${rCNV_bucket}/refs/GRCh37.centromeres_telomeres.bed.gz ./
zcat GRCh37.centromeres_telomeres.bed.gz \
| awk -v FS="\t" -v OFS="\t" -v d=${cen_tel_dist} \
  '{ if ($2-d<0) print $1, "0", $3+d; else print $1, $2-d, $3+d }' \
| cat - <( zcat refs/lit_GDs.*.bed.gz | cut -f1-3 | fgrep -v "#" ) \
| sort -Vk1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
> ${freq_code}.gene_scoring.training_mask.bed.gz
bedtools intersect -u \
  -a ${del_meta_stats} \
  -b ${freq_code}.gene_scoring.training_mask.bed.gz \
| cut -f1-4 \
| cat <( echo -e "#chr\tstart\tend\tgene" ) - \
| bgzip -c \
> ${freq_code}.gene_scoring.training_gene_blacklist.bed.gz

# Copy gene lists
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/genes/gene_lists ./

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

# Compute prior effect sizes
/opt/rCNV2/analysis/gene_scoring/estimate_prior_effect_sizes.R \
  ${del_meta_stats} \
  ${dup_meta_stats} \
  ${freq_code}.gene_scoring.training_gene_blacklist.bed.gz \
  gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
  gold_standard.haploinsufficient.genes.list \
  gold_standard.haplosufficient.genes.list \
  ${freq_code}.prior_estimation
awk -v FS="\t" '{ if ($1=="theta0" && $2=="DEL") print $3 }' \
  ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
> theta0_del.tsv
awk -v FS="\t" '{ if ($1=="theta0" && $2=="DUP") print $3 }' \
  ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
> theta0_dup.tsv
awk -v FS="\t" '{ if ($1=="theta1" && $2=="DEL") print $3 }' \
  ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
> theta1.tsv
awk -v FS="\t" '{ if ($1=="var1" && $2=="DEL") print $3 }' \
  ${freq_code}.prior_estimation.empirical_prior_estimates.tsv \
> var1.tsv

# Compute BFDP per gene
for CNV in DEL DUP; do
  # Set CNV-specific variables
  case $CNV in
    "DEL")
      theta0=$( cat theta0_del.tsv )
      statsbed=${del_meta_stats}
      ;;
    "DUP")
      theta0=$( cat theta0_dup.tsv )
      statsbed=${dup_meta_stats}
      ;;
  esac

  # Compute BF & BFDR for all genes
  /opt/rCNV2/analysis/gene_scoring/calc_gene_bfs.py \
    --theta0 $theta0 \
    --theta1 $( cat theta1.tsv ) \
    --var0 $( cat var1.tsv ) \
    --prior ${prior_frac} \
    --blacklist ${freq_code}.gene_scoring.training_gene_blacklist.bed.gz \
    --outfile ${freq_code}.$CNV.gene_abfs.tsv \
    "$statsbed"
done




# Score genes for a single model & CNV type
# Dev/test parameters
freq_code="rCNV"
rCNV_bucket="${rCNV_bucket}"
CNV="DEL"
BFDP_stats="${freq_code}.${CNV}.gene_abfs.tsv"
blacklist="${freq_code}.gene_scoring.training_gene_blacklist.bed.gz"
underpowered_genes="underpowered_genes.test.bed.gz"
gene_features="gencode.v19.canonical.pext_filtered.all_features.eigenfeatures.bed.gz"
raw_gene_features="gencode.v19.canonical.pext_filtered.all_features.bed.gz"
max_true_bfdp=0.2
min_false_bfdp=0.8
model="logit"
elnet_alpha=0.1
elnet_l1_l2_mix=1

# Copy centromeres bed
gsutil -m cp ${rCNV_bucket}/refs/GRCh37.centromeres_telomeres.bed.gz ./

# Merge blacklist and list of underpowered genes
zcat ${blacklist} ${underpowered_genes} \
| grep -ve '^#' \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
| uniq \
| cat <( zcat ${blacklist} | grep -e '^#' ) - \
| bgzip -c  \
> blacklist_plus_underpowered.bed.gz

# Score genes
/opt/rCNV2/analysis/gene_scoring/score_genes.py \
  --centromeres GRCh37.centromeres_telomeres.bed.gz \
  --blacklist blacklist_plus_underpowered.bed.gz \
  --model ${model} \
  --max-true-bfdp ${max_true_bfdp} \
  --min-false-bfdp ${min_false_bfdp} \
  --regularization-alpha ${elnet_alpha} \
  --regularization-l1-l2-mix ${elnet_l1_l2_mix} \
  --outfile ${freq_code}.${CNV}.gene_scores.${model}.tsv \
  ${BFDP_stats} \
  ${gene_features}




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
    for model in logit svm randomforest lda naivebayes sgd neuralnet; do
      echo $model
      echo ${freq_code}.$CNV.gene_scores.$model.tsv
    done | paste - -
  done > ${freq_code}.$CNV.model_evaluation.input.tsv
done
/opt/rCNV2/analysis/gene_scoring/compare_models.R \
  ${freq_code}.DEL.model_evaluation.input.tsv \
  ${freq_code}.DUP.model_evaluation.input.tsv \
  gene_lists/gold_standard.haploinsufficient.genes.list \
  gene_lists/gold_standard.haplosufficient.genes.list \
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
  DEL.roc_truth_sets.tsv \
  DUP.roc_truth_sets.tsv \
  gene_lists/gold_standard.haplosufficient.genes.list \
  ${freq_code}_gene_scoring_QC_plots/${freq_code}_gene_score_qc

