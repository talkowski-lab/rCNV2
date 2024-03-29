#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Prediction of dosage sensitivity scores for all genes


# Launch docker image
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv
gcloud auth login


# Set parameters
export rCNV_bucket="gs://rcnv_project"


# Copy all filtered CNV data, gene coordinates, and other references 
# from the project Google Bucket (note: requires permissions)
mkdir cleaned_cnv/ phenos/ refs/
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/cnv/* cleaned_cnv/
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/phenotypes/filtered/meta*.cleaned_phenos.txt phenos/
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/genes ./
gsutil -m cp \
  ${rCNV_bucket}/refs/GRCh37.*.bed.gz \
  ${rCNV_bucket}/analysis/analysis_refs/* \
  ${rCNV_bucket}/analysis/paper/data/large_segments/lit_GDs.*.bed.gz \
  ${rCNV_bucket}/analysis/paper/data/large_segments/loose_unclustered_nahr_regions.bed.gz \
  refs/


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
min_cds_ovr_del=0.84
min_cds_ovr_dup=0.84
max_genes_per_cnv=37
meta_model_prefix="fe"
max_standard_error=1
exclusion_bed=refs/gencode.v19.canonical.pext_filtered.cohort_exclusion.bed.gz
winsorize_meta_z=1.0
meta_min_cases=300
prior_frac=0.137


# Create training excludelist: Remove all genes within known NAHR-mediated genomic disorder regions
# (Note: cen/tel exclusion no longer applied — no evidence these genes will bias results)
zcat refs/lit_GDs.*.bed.gz | cut -f1-3 | fgrep -v "#" \
| bedtools intersect -r -f 0.25 -u -a - \
  -b refs/loose_unclustered_nahr_regions.bed.gz \
| sort -Vk1,1 -k2,2n -k3,3n \
| bedtools merge -i - \
> ${freq_code}.gene_scoring.training_mask.bed.gz
gsutil -m cp \
  ${rCNV_bucket}/analysis/gene_burden/HP0000118/rCNV/stats/meta1.HP0000118.rCNV.DEL.gene_burden.stats.bed.gz \
  ./
bedtools intersect -u \
  -a meta1.HP0000118.rCNV.DEL.gene_burden.stats.bed.gz \
  -b ${freq_code}.gene_scoring.training_mask.bed.gz \
| cut -f1-4 \
| cat <( echo -e "#chr\tstart\tend\tgene" ) - \
| bgzip -c \
> ${freq_code}.gene_scoring.training_gene_excludelist.bed.gz

# Copy training excludelist to project bucket (note: requires permissions)
gsutil -m cp \
  ${freq_code}.gene_scoring.training_gene_excludelist.bed.gz \
  ${rCNV_bucket}/analysis/gene_scoring/refs/


# Create training gene lists
# Define list of high-confidence haploinsufficient genes
cat \
  genes/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
  <( cat genes/gene_lists/DDG2P.hc_lof.genes.list genes/gene_lists/ClinGen.hc_haploinsufficient.genes.list | sort | uniq ) \
  genes/gene_lists/gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.low_expression_invariant.genes.list \
  genes/gene_lists/gnomad_sv.v2.1.nonneuro.no_lof_dels.genes.list \
| sort | uniq -c \
| awk '{ if ($1>=3) print $2 }' \
| fgrep -wf - genes/gene_lists/gencode.v19.canonical.pext_filtered.genes.list \
> gold_standard.haploinsufficient.genes.list

# Define list of high-confidence haplosufficient genes
cat \
  genes/gene_lists/gnomad.v2.1.1.mutation_tolerant.genes.list \
  <( cat genes/gene_lists/HP0000118.HPOdb.genes.list genes/gene_lists/DDG2P*.genes.list genes/gene_lists/ClinGen*.genes.list \
     | fgrep -wvf - genes/gene_lists/gencode.v19.canonical.pext_filtered.genes.list | sort | uniq ) \
  genes/gene_lists/gnomad_sv.v2.1.nonneuro.has_lof_dels.genes.list \
  genes/gene_lists/gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.low_expression_variable.genes.list \
| sort | uniq -c \
| awk '{ if ($1>=4) print $2 }' \
| fgrep -wf - genes/gene_lists/gencode.v19.canonical.pext_filtered.genes.list \
> gold_standard.haplosufficient.genes.list

# Define list of high-confidence triplosensitive genes
cat \
  genes/gene_lists/gnomad.v2.1.1.mis_constrained.genes.list \
  <( cat genes/gene_lists/DDG2P.hc_gof.genes.list genes/gene_lists/DDG2P.hc_other.genes.list genes/gene_lists/ClinGen.all_triplosensitive.genes.list | sort | uniq ) \
  genes/gene_lists/gnomad_sv.v2.1.nonneuro.no_cg_dups.genes.list \
  genes/gene_lists/gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.high_expression_invariant.genes.list \
| sort | uniq -c \
| awk '{ if ($1>=3) print $2 }' \
| fgrep -wf - genes/gene_lists/gencode.v19.canonical.pext_filtered.genes.list \
> gold_standard.triplosensitive.genes.list

# Define list of high-confidence triploinsensitive genes
cat \
  genes/gene_lists/gnomad.v2.1.1.mutation_tolerant.genes.list \
  <( cat genes/gene_lists/HP0000118.HPOdb.genes.list genes/gene_lists/DDG2P*.genes.list genes/gene_lists/ClinGen*.genes.list \
     | fgrep -wvf - genes/gene_lists/gencode.v19.canonical.pext_filtered.genes.list | sort | uniq ) \
  genes/gene_lists/gnomad_sv.v2.1.nonneuro.has_cg_dups.genes.list \
  genes/gene_lists/gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.high_expression_variable.genes.list \
| sort | uniq -c \
| awk '{ if ($1>=4) print $2 }' \
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
# subsequent code block (number of genes per CNV) with a very large --max-genes-for-summary, 
# then running this block a second time after updating --max-genes-for-summary to
# the optimal value determined by optimize_genes_per_cnv.R
while read meta cohorts; do
  echo $meta
  # Process deletions
  echo DEL
  /opt/rCNV2/analysis/gene_scoring/distill_optimization_data.py \
    --genes-per-cnv optimization_data/${freq_code}.DEL.$meta.genes_per_cnv.tsv.gz \
    --hpos <( cut -f2 refs/test_phenotypes.list ) \
    --positive-truth-genes gold_standard.haploinsufficient.genes.list \
    --negative-truth-genes gold_standard.haplosufficient.genes.list \
    --exclude-genes <( zcat ${freq_code}.gene_scoring.training_gene_excludelist.bed.gz \
                       | fgrep -v "#" | cut -f4 ) \
    --max-genes-for-summary 37 \
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
    --exclude-genes <( zcat ${freq_code}.gene_scoring.training_gene_excludelist.bed.gz \
                       | fgrep -v "#" | cut -f4 ) \
    --max-genes-for-summary 37 \
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
  --max-eval 100 \
  optimize_genes_per_cnv.DEL.input.tsv \
  optimize_genes_per_cnv.DUP.input.tsv \
  optimize_genes_per_cnv.results.pdf


# Recompute association stats per cohort
for contig in $( seq 1 22 ); do

  # Extract contig of interest from GTF
  tabix ${gtf} ${contig} | bgzip -c > ${contig}.gtf.gz

  # Create string of HPOs to keep
  keep_hpos=$( cat refs/rCNV2.hpos_by_severity.developmental.list | paste -s -d\; )

  # Iterate over metacohorts to compute single-cohort stats
  while read meta cohorts; do
    echo $meta

    # Exclude all NAHR CNVs from burden testing
    bedtools intersect -v -r -f 0.5 \
      -a cleaned_cnv/$meta.${freq_code}.bed.gz \
      -b refs/loose_unclustered_nahr_regions.bed.gz \
    | bgzip -c \
    > cleaned_cnv/$meta.${freq_code}.no_NAHR.bed.gz
    tabix -f cleaned_cnv/$meta.${freq_code}.no_NAHR.bed.gz

    # Set metacohort-specific parameters
    cnv_bed="cleaned_cnv/$meta.${freq_code}.no_NAHR.bed.gz"
    effective_case_n=$( fgrep -w $meta refs/rCNV2.hpos_by_severity.developmental.counts.tsv | cut -f2 )

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
      /opt/rCNV2/analysis/generic_scripts/fisher_test_single_cohort.R \
        --pheno-table ${metacohort_sample_table} \
        --cohort-name $meta \
        --effective-case-n $effective_case_n \
        --keep-n-columns 4 \
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
  /opt/rCNV2/analysis/generic_scripts/meta_analysis.R \
    --model ${meta_model_prefix} \
    --conditional-exclusion $exclusion_bed \
    --p-is-neg-log10 \
    --spa \
    --spa-exclude /opt/rCNV2/refs/lit_GDs.all.${CNV}.bed.gz \
    --winsorize ${winsorize_meta_z} \
    --min-cases ${meta_min_cases} \
    --keep-n-columns 4 \
    ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.input.txt \
    ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed
  bgzip -f ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed
  tabix -p bed -f ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz

  # Extract list of genes with standard error < max_standard_error (to be used as excludelist later for training)
  /opt/rCNV2/analysis/gene_scoring/get_underpowered_genes.R \
    --max-se "${max_standard_error}" \
    ${prefix}.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz \
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
del_meta_stats="rCNV2_analysis_d2.rCNV.DEL.gene_burden.meta_analysis.stats.bed.gz"
dup_meta_stats="rCNV2_analysis_d2.rCNV.DUP.gene_burden.meta_analysis.stats.bed.gz"


# Create CNV-type-specific gene excludelists
for CNV in DEL DUP; do
  zcat *.$CNV.gene_burden.underpowered_genes.bed.gz \
    ${freq_code}.gene_scoring.training_gene_excludelist.bed.gz \
  | fgrep -v "#" | cut -f4 | sort -V | uniq \
  > ${freq_code}.$CNV.training_excludelist.genes.list
done

# Prepare files for on-the fly meta-analysis of CNV effect sizes
echo -e "cohort\tn_case\tn_control\tDEL_path\tDUP_path" \
> prior_estimation.meta_inputs.tsv
zcat optimization_data/rCNV.DEL.meta1.genes_per_cnv.tsv.gz | head -n1 \
> optimization_data/opt_data.header.tsv
while read meta cohorts; do
  # Subset CNVs to developmental cases & controls
  for CNV in DEL DUP; do
    zcat optimization_data/rCNV.$CNV.$meta.genes_per_cnv.tsv.gz \
    | fgrep -wf <( cat refs/rCNV2.hpos_by_severity.developmental.list \
                       <( echo "HEALTHY_CONTROL" ) ) \
    | sort -Vk1,1 | cat optimization_data/opt_data.header.tsv - | gzip -c \
    > optimization_data/rCNV.$meta.annotated_developmental_and_control.$CNV.tsv.gz
  done

  # Write info to meta-analysis input
  for dummy in 1; do
    fgrep -w $meta refs/rCNV2.hpos_by_severity.developmental.counts.tsv
    cidx=$( head -n1 refs/HPOs_by_metacohort.table.tsv | sed 's/\t/\n/g' \
            | awk -v OFS="\t" '{ print NR, $0 }' | fgrep -w $meta | cut -f1 )
    fgrep -w "HEALTHY_CONTROL" refs/HPOs_by_metacohort.table.tsv | cut -f$cidx
    for CNV in DEL DUP; do
      echo optimization_data/rCNV.$meta.annotated_developmental_and_control.$CNV.tsv.gz
    done
  done | paste -s >> prior_estimation.meta_inputs.tsv
done < <( fgrep -v "mega" refs/rCNV_metacohort_list.txt )

# Compute prior effect sizes
/opt/rCNV2/analysis/gene_scoring/estimate_prior_effect_sizes.R \
  prior_estimation.meta_inputs.tsv \
  ${freq_code}.gene_scoring.training_gene_excludelist.bed.gz \
  genes/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list \
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
    --blacklist ${freq_code}.$CNV.training_excludelist.genes.list \
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
excludelist="${freq_code}.gene_scoring.training_gene_excludelist.bed.gz"
underpowered_genes="rCNV2_analysis_d2.rCNV.DEL.gene_burden.underpowered_genes.bed.gz"
gene_features="gencode.v19.canonical.pext_filtered.all_features.no_variation.eigenfeatures.bed.gz"
raw_gene_features="gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz"
max_true_bfdp=0.5
min_false_bfdp=0.5
model="logit"

# Copy necessary references
gsutil -m cp \
  ${rCNV_bucket}/refs/GRCh37.centromeres_telomeres.bed.gz \
  ${rCNV_bucket}/analysis/gene_scoring/gene_lists/*genes.list \
  ./

# Merge excludelist and list of underpowered genes
zcat ${excludelist} ${underpowered_genes} \
| grep -ve '^#' \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
| uniq \
| cat <( zcat ${excludelist} | grep -e '^#' ) - \
| bgzip -c  \
> excludelist_plus_underpowered.bed.gz

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
  --blacklist excludelist_plus_underpowered.bed.gz \
  --true-positives $true_pos \
  --true-negatives $true_neg \
  --model ${model} \
  --max-true-bfdp ${max_true_bfdp} \
  --min-false-bfdp ${min_false_bfdp} \
  --chromsplit \
  --no-out-of-sample-prediction \
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
    for model in logit svm randomforest lda naivebayes neuralnet gbdt knn ensemble; do
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

