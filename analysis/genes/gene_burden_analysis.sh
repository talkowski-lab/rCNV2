#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2019-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control CNV burdens for canonical protein-coding genes


# Launch docker image
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv
gcloud auth login


# Test/dev parameters
# Params for NDDs
hpo="HP:0012759"
prefix="HP0012759"
meta="meta1"
# General params
freq_code="rCNV"
CNV="DEL"
phenotype_list="refs/test_phenotypes.list"
metacohort_list="refs/rCNV_metacohort_list.txt"
metacohort_sample_table="refs/HPOs_by_metacohort.table.tsv"
gtf="genes/gencode.v19.canonical.pext_filtered.gtf.gz"
rCNV_bucket="gs://rcnv_project"
pad_controls=0
min_cds_ovr_del=0.02
min_cds_ovr_dup=0.84
max_genes_per_cnv=20000
p_cutoff=0.000002973712
max_manhattan_neg_log10_p=30
n_pheno_perms=50
meta_model_prefix="fe"
i=1
winsorize_meta_z=0.99
meta_min_cases=300


# Copy all filtered CNV data, gene coordinates, and other references 
# from the project Google Bucket (note: requires permissions)
mkdir cleaned_cnv/
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/cnv/* cleaned_cnv/
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/genes ./
mkdir refs/
gsutil -m cp ${rCNV_bucket}/refs/GRCh37.*.bed.gz refs/
gsutil -m cp ${rCNV_bucket}/analysis/analysis_refs/* refs/


# Count CNVs in cases and controls per phenotype, split by metacohort and CNV type
# Iterate over phenotypes
while read pheno hpo; do
  # Set HPO-specific parameters
  descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" '{ print $2 }' )
  zcat refs/gencode.v19.canonical.constrained.bed.gz \
  | fgrep -wf genes/gene_lists/${prefix}.HPOdb.constrained.genes.list \
  | cat <( echo -e "#chr\tstart\tend\tgene" ) - \
  > ${prefix}.highlight_regions.bed

  # Iterate over metacohorts
  while read meta cohorts; do
    echo $meta

    # Set metacohort-specific parameters
    cnv_bed="cleaned_cnv/$meta.${freq_code}.bed.gz"
    meta_idx=$( head -n1 "${metacohort_sample_table}" \
                | sed 's/\t/\n/g' \
                | awk -v meta="$meta" '{ if ($1==meta) print NR }' )
    ncase=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" -v meta_idx="$meta_idx" '{ print $meta_idx }' \
             | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
    nctrl=$( fgrep -w "HEALTHY_CONTROL" "${metacohort_sample_table}" \
             | awk -v FS="\t" -v meta_idx="$meta_idx" '{ print $meta_idx }' \
             | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
    title="$descrip (${hpo})\n$ncase cases vs $nctrl controls in '$meta' cohort"

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
        --pad-controls ${pad_controls} \
        --min-cds-ovr $min_cds_ovr \
        --max-genes ${max_genes_per_cnv} \
        -t $CNV \
        --hpo ${hpo} \
        --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
        --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
        --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
        -z \
        --verbose \
        -o "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.bed.gz" \
        "$cnv_bed" \
        ${gtf}
      tabix -f "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.bed.gz"

      # Perform burden test
      /opt/rCNV2/analysis/generic_scripts/fisher_test_single_cohort.R \
        --pheno-table ${metacohort_sample_table} \
        --cohort-name $meta \
        --case-hpo ${hpo} \
        --keep-n-columns 4 \
        --bgzip \
        "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.bed.gz" \
        "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz"
      tabix -f "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz"

      # Generate Manhattan & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --p-col-name "fisher_neg_log10_p" \
        --p-is-neg-log10 \
        --max-neg-log10-p ${max_manhattan_neg_log10_p} \
        --cutoff ${p_cutoff} \
        --highlight-bed "${prefix}.highlight_regions.bed" \
        --highlight-name "Constrained genes associated with ${hpo}" \
        --label-prefix "$CNV" \
        --title "$title" \
        "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz" \
        "$meta.${prefix}.${freq_code}.$CNV.gene_burden"
    done

    # Generate Miami & QQ plots
    /opt/rCNV2/utils/plot_manhattan_qq.R \
      --miami \
      --p-col-name "fisher_neg_log10_p" \
      --p-is-neg-log10 \
      --max-neg-log10-p ${max_manhattan_neg_log10_p} \
      --cutoff ${p_cutoff} \
      --highlight-bed "${prefix}.highlight_regions.bed" \
      --highlight-name "Constrained genes associated with ${hpo}" \
      --label-prefix "DUP" \
      --highlight-bed-2 "${prefix}.highlight_regions.bed" \
      --highlight-name-2 "Constrained genes associated with ${hpo}" \
      --label-prefix-2 "DEL" \
      --title "$title" \
      "$meta.${prefix}.${freq_code}.DUP.gene_burden.stats.bed.gz" \
      "$meta.${prefix}.${freq_code}.DEL.gene_burden.stats.bed.gz" \
      "$meta.${prefix}.${freq_code}.gene_burden"
  done < ${metacohort_list}
done < ${phenotype_list}






### Determine probe density-based conditional exclusion list
# NOTE: This code must be run using a DIFFERENT DOCKER: us.gcr.io/broad-dsmap/athena-cloud
# Test/dev parameters
genes_bed="$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz" #Note: in the WDL, this file is extracted as the first stats bed from the array (as a test case)
# For convenience, this can be downloaded from a precomputed file here (but will be generated dynamically in the WDL)
gsutil -m cp ${rCNV_bucket}/freeze_pre_rev1/analysis_freeze_pre_rev1/gene_burden/${prefix}/${freq_code}/stats/$genes_bed ./
gtf_prefix="gencode.v19.canonical" #Note: this can be inferred in WDL as basename(gtf, ".gtf.gz")
min_probes_per_gene=10
min_frac_controls_probe_exclusion=0.9
metacohort_list="refs/rCNV_metacohort_list.txt"
rCNV_bucket="gs://rcnv_project"
freq_code="rCNV"

# Download probesets (note: requires permissions)
gsutil -m cp -r \
  ${rCNV_bucket}/cleaned_data/control_probesets \
  ./

# Subset gene coordinates to minimal BED4
zcat ${genes_bed} | cut -f1-4 | bgzip -c > gene_coords.bed.gz

# Clone rCNV2 repo (not present in athena-cloud Docker)
cd opt && \
git clone https://github.com/talkowski-lab/rCNV2.git && \
cd -

# Compute conditional cohort exclusion mask
for file in control_probesets/*bed.gz; do
  echo -e "$file\t$( basename $file | sed 's/\.bed\.gz//g' )"
done > probeset_tracks.tsv
/opt/rCNV2/data_curation/other/probe_based_exclusion.py \
  --outfile ${gtf_prefix}.cohort_exclusion.bed.gz \
  --probecounts-outfile ${gtf_prefix}.probe_counts.bed.gz \
  --control-mean-counts-outfile ${gtf_prefix}.mean_probe_counts_per_cohort.bed.gz \
  --frac-pass-outfile ${gtf_prefix}.frac_passing.bed.gz \
  --min-probes ${min_probes_per_gene} \
  --min-frac-samples ${min_frac_controls_probe_exclusion} \
  --keep-n-columns 4 \
  --bgzip \
  gene_coords.bed.gz \
  probeset_tracks.tsv \
  control_probesets/rCNV.control_counts_by_array.tsv \
  <( fgrep -v mega ${metacohort_list} )

# Estimate number of effective tests while requiring at least two cohorts to have
# adequate probe density for a gene to be evaluated
zcat ${gtf_prefix}.cohort_exclusion.bed.gz | sed 's/;/\t/g' \
| awk -v FS="\t" -v OFS="\t" '{ if (NF<=7) print $1, $2, $3 }' \
| fgrep -v "#" | wc -l




# Permute CNVs in cases and controls per phenotype, split by metacohort and CNV type
# Iterate over phenotypes
while read prefix hpo; do
  for i in $( seq 1 ${n_pheno_perms} ); do

    # Shuffle phenotypes for each metacohort CNV dataset, and restrict CNVs to
    # phenotype of relevance
    if ! [ -e shuffled_cnv/ ]; then
      mkdir shuffled_cnv/
    fi
    yes $i | head -n1000000 > seed_$i.txt
    while read meta cohorts; do

      cnvbed="cleaned_cnv/$meta.${freq_code}.bed.gz"

      # Shuffle phenotypes, matching by CNV type
      for CNV in DEL DUP; do
        zcat $cnvbed | fgrep -w $CNV \
        > cnv_subset.bed
        paste <( cut -f1-5 cnv_subset.bed ) \
              <( cut -f6 cnv_subset.bed | shuf --random-source seed_$i.txt ) \
        | awk -v hpo=${hpo} '{ if ($NF ~ "HEALTHY_CONTROL" || $NF ~ hpo) print $0 }'
        rm cnv_subset.bed
      done \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | cat <( tabix -H $cnvbed ) - \
      | bgzip -c \
      > shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed.gz
      tabix -f shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed.gz

    done < <( fgrep -v "mega" ${metacohort_list} )

    # Iterate over CNV types
    for CNV in DEL DUP; do

      # Set CNV-specific parameters
      case "$CNV" in
        DEL)
          min_cds_ovr=${min_cds_ovr_del}
          ;;
        DUP)
          min_cds_ovr=${min_cds_ovr_dup}
          ;;
      esac

      # Perform association test for each metacohort
      while read meta cohorts; do

        echo -e "[$( date )] Starting permutation $i for $CNV in ${hpo} from $meta...\n"

        cnv_bed="shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed.gz"

        # Count CNVs
        /opt/rCNV2/analysis/genes/count_cnvs_per_gene.py \
          --pad-controls ${pad_controls} \
          --min-cds-ovr $min_cds_ovr \
          --max-genes ${max_genes_per_cnv} \
          -t $CNV \
          --hpo ${hpo} \
          --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
          --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
          --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
          -z \
          --verbose \
          -o "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.bed.gz" \
          $cnv_bed \
          ${gtf}
        tabix -f "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.bed.gz"

        # Perform burden test
        /opt/rCNV2/analysis/generic_scripts/fisher_test_single_cohort.R \
          --pheno-table ${metacohort_sample_table} \
          --cohort-name $meta \
          --case-hpo ${hpo} \
          --keep-n-columns 4 \
          --bgzip \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.bed.gz" \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz"
        tabix -f "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz"

      done < <( fgrep -v "mega" ${metacohort_list} )

      # Perform meta-analysis
      while read meta cohorts; do
        echo -e "$meta\t$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz"
      done < <( fgrep -v mega ${metacohort_list} ) \
      > ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.input.txt
      /opt/rCNV2/analysis/generic_scripts/meta_analysis.R \
        --model ${meta_model_prefix} \
        --conditional-exclusion ${exclusion_bed} \
        --p-is-neg-log10 \
        --spa \
        --spa-exclude /opt/rCNV2/refs/lit_GDs.all.$CNV.bed.gz \
        --winsorize ${winsorize_meta_z} \
        --min-cases ${meta_min_cases} \
        --keep-n-columns 4 \
        ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.input.txt \
        ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.perm_$i.bed
      bgzip -f ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.perm_$i.bed
      tabix -f ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.perm_$i.bed.gz

  done
done < ${phenotype_list}

# Gather all permutation results and compute FDR CDFs
mkdir perm_res/
while read prefix hpo; do
  gsutil -m cp \
    "${rCNV_bucket}/analysis/gene_burden/$prefix/${freq_code}/permutations/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.perm_*.bed.gz" \
    perm_res/
  for i in $( seq 1 ${n_pheno_perms} ); do
    p_idx=$( zcat perm_res/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.perm_$i.bed.gz \
             | sed -n '1p' | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' \
             | fgrep -w ${p_val_column_name} | cut -f2 )
    zcat perm_res/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.perm_$i.bed.gz \
    | grep -ve '^#' \
    | awk -v p_idx=$p_idx '{ print $(p_idx) }' \
    | cat <( echo "$prefix.${CNV}.$i" ) - \
    > perm_res/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.permuted_p_values.$i.txt
  done
  rm perm_res/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.perm_*.bed.gz
done < ${phenotype_list}
paste perm_res/*.gene_burden.meta_analysis.permuted_p_values.*.txt \
| gzip -c \
> ${freq_code}.${CNV}.permuted_pval_matrix.txt.gz

# Analyze p-values and compute FDR
/opt/rCNV2/analysis/sliding_windows/calc_empirical_fdr.R \
  --cnv ${CNV} \
  --fdr-target ${meta_p_cutoff} \
  --linear-fit \
  --flat-ladder \
  --plot gene_burden.${freq_code}.${CNV}.${fdr_table_suffix}_permutation_results.png \
  ${freq_code}.${CNV}.permuted_pval_matrix.txt.gz \
  ${metacohort_sample_table} \
  gene_burden.${freq_code}.${CNV}.${fdr_table_suffix}

# Also produce an optional table of flat Bonferroni P-value cutoffs
awk -v pval=${p_cutoff} -v FS="\t" -v OFS="\t" \
  '{ print $1, pval }' ${phenotype_list} \
> gene_burden.${freq_code}.${CNV}.bonferroni_pval.hpo_cutoffs.tsv





### Run meta-analysis for each phenotype
# # Test/dev parameters (seizures)
# hpo="HP:0001250"
# prefix="HP0001250"
# Test/dev parameters (NDDs)
hpo="HP:0012759"
prefix="HP0012759"
freq_code="rCNV"
cnv="DEL"
phenotype_list="refs/test_phenotypes.list"
metacohort_list="refs/rCNV_metacohort_list.txt"
metacohort_sample_table="refs/HPOs_by_metacohort.table.tsv"
rCNV_bucket="gs://rcnv_project"
p_cutoff=0.000002973712
meta_p_cutoff=0.000002973712
max_manhattan_neg_log10_p=30
meta_model_prefix="fe"
exclusion_bed="gencode.v19.canonical.cohort_exclusion.bed.gz" #Note: this file must be generated above

# Copy necessary data for local testing (without running the above -- this is not in the WDL)
gsutil -m cp \
  ${rCNV_bucket}/analysis/analysis_refs/gene_burden.${freq_code}.*.bonferroni_pval.hpo_cutoffs.tsv \
  ./
gsutil -m cp \
  ${rCNV_bucket}/analysis/gene_burden/${prefix}/${freq_code}/stats/meta**.stats.bed.gz* \
  ./


while read prefix hpo; do

  # Get metadata for meta-analysis (while accounting for cohorts below inclusion criteria)
  last_cohort_col=$( head -n1 "${metacohort_sample_table}" | awk '{ print NF-1 }' )
  keep_cols=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
               | cut -f4-$last_cohort_col \
               | sed 's/\t/\n/g' \
               | awk -v min_n=${meta_min_cases} '{ if ($1>=min_n) print NR+3 }' \
               | paste -s -d, )
  ncase=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
           | cut -f$keep_cols \
           | sed 's/\t/\n/g' \
           | awk '{ sum+=$1 }END{ print sum }' )
  nctrl=$( fgrep -w "HEALTHY_CONTROL" "${metacohort_sample_table}" \
           | cut -f$keep_cols \
           | sed 's/\t/\n/g' \
           | awk '{ sum+=$1 }END{ print sum }' )
  descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" '{ print $2 }' )
  title="$descrip (${hpo})\nMeta-analysis of $ncase cases and $nctrl controls"
  DEL_p_cutoff=$( awk -v hpo=${prefix} '{ if ($1==hpo) print $2 }' \
                  gene_burden.${freq_code}.DEL.bonferroni_pval.hpo_cutoffs.tsv )
  DUP_p_cutoff=$( awk -v hpo=${prefix} '{ if ($1==hpo) print $2 }' \
                    gene_burden.${freq_code}.DUP.bonferroni_pval.hpo_cutoffs.tsv )

  # Set HPO-specific parameters
  descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" '{ print $2 }' )
  zcat refs/gencode.v19.canonical.constrained.bed.gz \
  | fgrep -wf genes/gene_lists/${prefix}.HPOdb.constrained.genes.list \
  | cat <( echo -e "#chr\tstart\tend\tgene" ) - \
  > ${prefix}.highlight_regions.bed


  # Run meta-analysis for each CNV type
  for CNV in DEL DUP; do
    # Set CNV-specific parameters
    case "$CNV" in
      DEL)
        meta_p_cutoff=$DEL_p_cutoff
        ;;
      DUP)
        meta_p_cutoff=$DUP_p_cutoff
        ;;
    esac

    # Perform meta-analysis for unweighted CNVs
    while read meta cohorts; do
      echo -e "$meta\t$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz"
    done < <( fgrep -v mega ${metacohort_list} ) \
    > ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.input.txt
    /opt/rCNV2/analysis/generic_scripts/meta_analysis.R \
      --model ${meta_model_prefix} \
      --conditional-exclusion ${exclusion_bed} \
      --p-is-neg-log10 \
      --spa \
      --spa-exclude /opt/rCNV2/refs/lit_GDs.all.$CNV.bed.gz \
      --winsorize ${winsorize_meta_z} \
      --min-cases ${meta_min_cases} \
      --keep-n-columns 4 \
      ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.input.txt \
      ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed
    bgzip -f ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed
    tabix -f ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed.gz

    # Generate Manhattan & QQ plots
    /opt/rCNV2/utils/plot_manhattan_qq.R \
      --p-col-name "meta_neg_log10_p" \
      --p-is-neg-log10 \
      --max-neg-log10-p ${max_manhattan_neg_log10_p} \
      --cutoff $meta_p_cutoff \
      --highlight-bed "${prefix}.highlight_regions.bed" \
      --highlight-name "Constrained genes associated with ${hpo}" \
      --label-prefix "$CNV" \
      --title "$title" \
      "${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed.gz" \
      "${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis"
  done

  # Generate Miami & QQ plots
  /opt/rCNV2/utils/plot_manhattan_qq.R \
    --miami \
    --p-col-name "meta_neg_log10_p" \
    --p-is-neg-log10 \
    --max-neg-log10-p ${max_manhattan_neg_log10_p} \
    --cutoff $DUP_p_cutoff \
    --highlight-bed "${prefix}.highlight_regions.bed" \
    --highlight-name "Constrained genes associated with ${hpo}" \
    --label-prefix "DUP" \
    --cutoff-2 $DEL_p_cutoff \
    --highlight-bed-2 "${prefix}.highlight_regions.bed" \
    --highlight-name-2 "Constrained genes associated with ${hpo}" \
    --label-prefix-2 "DEL" \
    --title "$title" \
    "${prefix}.${freq_code}.DUP.gene_burden.meta_analysis.stats.bed.gz" \
    "${prefix}.${freq_code}.DEL.gene_burden.meta_analysis.stats.bed.gz" \
    "${prefix}.${freq_code}.gene_burden.meta_analysis"

done < refs/test_phenotypes.list

# Collapse all meta-analysis p-values into single matrix for visualizing calibration
mkdir meta_res/
# Download data
while read prefix hpo; do
  for CNV in DEL DUP; do
    echo -e "${rCNV_bucket}/analysis/gene_burden/${prefix}/${freq_code}/stats/${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed.gz"
  done
done < ${phenotype_list} \
| gsutil -m cp -I meta_res/
# Primary p-values
p_val_column_name="meta_neg_log10_p"
while read prefix hpo; do
  echo -e "$prefix\n\n"
  for CNV in DEL DUP; do
      stats=meta_res/$prefix.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed.gz
      if [ -e $stats ]; then
        p_idx=$( zcat $stats | sed -n '1p' | sed 's/\t/\n/g' \
                 | awk -v OFS="\t" '{ print $1, NR }' \
                 | fgrep -w ${p_val_column_name} | cut -f2 )
        zcat $stats | grep -ve '^#' \
        | awk -v p_idx=$p_idx '{ print $(p_idx) }' \
        | cat <( echo "$prefix.$CNV" ) - \
        > meta_res/$prefix.${freq_code}.$CNV.gene_burden.meta_analysis.p_values.txt
      fi
  done
done < ${phenotype_list}
paste meta_res/*.${freq_code}.$CNV.gene_burden.meta_analysis.p_values.txt \
| gzip -c \
> ${freq_code}.observed_pval_matrix.txt.gz
# Secondary p-values
p_val_column_name="meta_neg_log10_p_secondary"
while read prefix hpo; do
  echo -e "$prefix\n\n"
  for CNV in DEL DUP; do
      stats=meta_res/$prefix.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed.gz
      if [ -e $stats ]; then
        p_idx=$( zcat $stats | sed -n '1p' | sed 's/\t/\n/g' \
                 | awk -v OFS="\t" '{ print $1, NR }' \
                 | fgrep -w ${p_val_column_name} | cut -f2 )
        zcat $stats | grep -ve '^#' \
        | awk -v p_idx=$p_idx '{ print $(p_idx) }' \
        | cat <( echo "$prefix.$CNV" ) - \
        > meta_res/$prefix.${freq_code}.$CNV.gene_burden.meta_analysis.secondary_p_values.txt
      fi
  done
done < ${phenotype_list}
paste meta_res/*.${freq_code}.$CNV.gene_burden.meta_analysis.secondary_p_values.txt \
| gzip -c \
> ${freq_code}.observed_pval_matrix_secondary.txt.gz
# DEV NOTE: these p-values can be visualized with plot_gene_burden_meta_analysis_p_values.R,
#           which is currently just a code snippet referencing local filepaths




# Fine-map significant gene blocks
# Test/dev parameters
freq_code="rCNV"
CNV="DEL"
gtf="genes/gencode.v19.canonical.pext_filtered.gtf.gz"
rCNV_bucket="gs://rcnv_project"
phenotype_list="refs/test_phenotypes.list"
metacohort_sample_table="refs/HPOs_by_metacohort.table.tsv"
meta_p_cutoffs_tsv="refs/gene_burden.rCNV.DEL.bonferroni_pval.hpo_cutoffs.tsv"
meta_secondary_p_cutoff=0.05
meta_nominal_cohorts_cutoff=2
FDR_cutoff=0.01
gene_features="gencode.v19.canonical.pext_filtered.all_features.no_variation.eigenfeatures.bed.gz"
finemap_output_label="all_features"
finemap_elnet_alpha=0.05
finemap_elnet_l1_l2_mix=0.75
finemap_cluster_distance=1000000
finemap_nonsig_distance=2000000
finemap_min_cnv_covariance=0.2
finemap_conf_pip=0.15
finemap_vconf_pip=0.85

# Copy association stats & other reference data from the project Google Bucket (note: requires permissions)
mkdir stats
gsutil -m cp \
  ${rCNV_bucket}/analysis/gene_burden/**.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz \
  stats/
gsutil -m cp -r \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists \
  ${rCNV_bucket}/analysis/paper/data/large_segments/clustered_nahr_regions.bed.gz \
  ${rCNV_bucket}/analysis/analysis_refs/rCNV2.hpos_by_severity.*list \
  ./
gsutil -m cp \
  ${rCNV_bucket}/analysis/gene_scoring/gene_lists/* \
  ./gene_lists/
gsutil -m cp -r \
  ${rCNV_bucket}/cleaned_data/cnv/* \
  cleaned_cnv/


# Make list of genes from predicted NAHR-mediated CNV regions for training exclusion
zcat clustered_nahr_regions.bed.gz | fgrep -v "#" \
| awk -v FS="\t" '{ if ($5>0) print $NF }' \
| sed 's/;/\n/g' | sort | uniq > nahr.genes.list

# Write tsv inputs
while read prefix hpo; do
  for wrapper in 1; do
    echo "$hpo"
    echo "stats/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz"
    awk -v x=$prefix -v FS="\t" '{ if ($1==x) print $2 }' ${meta_p_cutoffs_tsv}
  done | paste -s
done < ${phenotype_list} \
> ${freq_code}.${CNV}.gene_fine_mapping.stats_input.tsv

# Make lists of a priori "true causal" genes for null variance estimation
case ${CNV} in
  "DEL")
    echo "gene_lists/gold_standard.haploinsufficient.genes.list" > known_causal_gene_lists.tsv
    ;;
  "DUP")
    echo "gene_lists/gold_standard.triplosensitive.genes.list" > known_causal_gene_lists.tsv
    ;;
esac

# Compute gene-gene covariance
case "$CNV" in
  DEL)
    min_cds_ovr=${min_cds_ovr_del}
    ;;
  DUP)
    min_cds_ovr=${min_cds_ovr_dup}
    ;;
esac

# Parallelize per contig
for contig in $( zcat cleaned_cnv/mega.${freq_code}.bed.gz | fgrep -v "#" | cut -f1 | sort -V | uniq ); do
  # Subset files to contig of interest
  tabix -h ${gtf} ${contig} | bgzip -c > ${contig}.gtf.gz
  tabix -h cleaned_cnv/mega.${freq_code}.bed.gz ${contig} | bgzip -c > ${contig}_cnvs.bed.gz

  # Compute covariance
  /opt/rCNV2/analysis/genes/count_cnvs_per_gene.py \
    --min-cds-ovr $min_cds_ovr \
    --max-genes 20000 \
    -t ${CNV} \
    --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
    --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
    --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
    --covariance ${prefix}.${CNV}.gene_covariance.${contig}.tsv.gz \
    --verbose \
    -o /dev/null \
    ${contig}_cnvs.bed.gz \
    ${contig}.gtf.gz
done

# Find all covariance shards
find / -name "${prefix}.$CNV.gene_covariance.*.tsv.gz" \
> covariance_shards.$CNV.txt

# Write covariance header
zcat $( head -n1 covariance_shards.$CNV.txt ) | head -n1 \
> ${prefix}.$CNV.gene_covariance.tsv

# Merge & sort covariance shards
cat covariance_shards.$CNV.txt | zcat - | fgrep -v "#" | sort -Vk1,1 -k2,2V \
>> ${prefix}.$CNV.gene_covariance.tsv
gzip -f ${prefix}.$CNV.gene_covariance.tsv

# Run functional fine-mapping procedure
/opt/rCNV2/analysis/genes/finemap_genes.py \
  --cnv ${CNV} \
  --secondary-p-cutoff ${meta_secondary_p_cutoff} \
  --min-nominal ${meta_nominal_cohorts_cutoff} \
  --secondary-or-nominal \
  --fdr-q-cutoff ${FDR_cutoff} \
  --secondary-for-fdr \
  --regularization-alpha ${finemap_elnet_alpha} \
  --regularization-l1-l2-mix ${finemap_elnet_l1_l2_mix} \
  --distance ${finemap_cluster_distance} \
  --nonsig-distance ${finemap_nonsig_distance} \
  --covariance-tsv ${freq_code}.gene_cnv_covariance.tsv.gz \
  --min-covariance ${finemap_min_cnv_covariance} \
  --training-exclusion nahr.genes.list \
  --use-max-pip-per-gene \
  --confident-pip ${finemap_conf_pip} \
  --very-confident-pip ${finemap_vconf_pip} \
  --known-causal-gene-lists known_causal_gene_lists.tsv \
  --outfile ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.${finemap_output_label}.tsv \
  --all-genes-outfile ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.${finemap_output_label}.tsv \
  --naive-outfile ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.naive_priors.tsv \
  --genetic-outfile ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genetics_only.tsv \
  --coeffs-out ${freq_code}.${CNV}.gene_fine_mapping.logit_coeffs.${finemap_output_label}.tsv \
  --sig-genes-bed ${freq_code}.${CNV}.gene_fine_mapping.genes.${finemap_output_label}.bed \
  --sig-assoc-bed ${freq_code}.${CNV}.gene_fine_mapping.associations.${finemap_output_label}.bed \
  --hpo-credsets-bed ${freq_code}.${CNV}.gene_fine_mapping.credible_sets_per_hpo.${finemap_output_label}.bed \
  --joint-credsets-bed ${freq_code}.${CNV}.gene_fine_mapping.joint_credible_sets.${finemap_output_label}.bed \
  ${freq_code}.${CNV}.gene_fine_mapping.stats_input.tsv \
  ${gene_features} \
  ${metacohort_sample_table}

# # Merge & sort outputs for developmental & adult-onset subgroups
# mkdir finemapping_merged_outputs/
# cat <( fgrep -v "#" ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.${finemap_output_label}.developmental.tsv ) \
#     <( fgrep -v "#" ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.${finemap_output_label}.adult.tsv ) \
# | sort -nrk4,4 \
# | cat <( grep -e '^#' ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.${finemap_output_label}.developmental.tsv ) - \
# > finemapping_merged_outputs/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.${finemap_output_label}.tsv
# cat <( fgrep -v "#" ${freq_code}.${CNV}.final_genes.associations.${finemap_output_label}.developmental.bed ) \
#     <( fgrep -v "#" ${freq_code}.${CNV}.final_genes.associations.${finemap_output_label}.adult.bed ) \
# | sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V -k6,6V \
# | cat <( grep -e '^#' ${freq_code}.${CNV}.final_genes.associations.${finemap_output_label}.developmental.bed ) - \
# > finemapping_merged_outputs/${freq_code}.${CNV}.final_genes.associations.${finemap_output_label}.bed
# cat <( fgrep -v "#" ${freq_code}.${CNV}.final_genes.credible_sets.${finemap_output_label}.developmental.bed ) \
#     <( fgrep -v "#" ${freq_code}.${CNV}.final_genes.credible_sets.${finemap_output_label}.adult.bed ) \
# | sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V -k6,6V \
# | cat <( grep -e '^#' ${freq_code}.${CNV}.final_genes.credible_sets.${finemap_output_label}.developmental.bed ) - \
# > finemapping_merged_outputs/${freq_code}.${CNV}.final_genes.credible_sets.${finemap_output_label}.bed
# cat <( fgrep -v "#" ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.${finemap_output_label}.developmental.tsv ) \
#     <( fgrep -v "#" ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.${finemap_output_label}.adult.tsv ) \
# | sort -nrk4,4 \
# | cat <( grep -e '^#' ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.${finemap_output_label}.developmental.tsv ) - \
# > finemapping_merged_outputs/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.${finemap_output_label}.tsv
# cat <( fgrep -v "#" ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.naive_priors.developmental.tsv ) \
#     <( fgrep -v "#" ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.naive_priors.adult.tsv ) \
# | sort -nrk4,4 \
# | cat <( grep -e '^#' ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.naive_priors.developmental.tsv ) - \
# > finemapping_merged_outputs/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.naive_priors.tsv
# cat <( fgrep -v "#" ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genetics_only.developmental.tsv ) \
#     <( fgrep -v "#" ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genetics_only.adult.tsv ) \
# | sort -nrk4,4 \
# | cat <( grep -e '^#' ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genetics_only.developmental.tsv ) - \
# > finemapping_merged_outputs/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genetics_only.tsv
# cat <( fgrep -v "#" ${freq_code}.${CNV}.final_genes.genes.${finemap_output_label}.developmental.bed \
#        | awk -v OFS="\t" '{ print $0, "developmental" }' ) \
#     <( fgrep -v "#" ${freq_code}.${CNV}.final_genes.genes.${finemap_output_label}.adult.bed \
#        | awk -v OFS="\t" '{ print $0, "adult" }' ) \
# | sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V -k6,6V \
# | cat <( grep -e "^#" ${freq_code}.${CNV}.final_genes.genes.${finemap_output_label}.developmental.bed \
#        | awk -v OFS="\t" '{ print $0, "finemapping_model" }' ) \
# > finemapping_merged_outputs/${freq_code}.${CNV}.final_genes.genes.${finemap_output_label}.bed


# Plot results of fine-mapping
# Test/dev parameters
freq_code="rCNV"
CNV="DEL"
phenotype_list="refs/test_phenotypes.list"
raw_features_genomic="gencode.v19.canonical.pext_filtered.genomic_features.bed.gz"
raw_features_expression="gencode.v19.canonical.pext_filtered.expression_features.bed.gz"
raw_features_chromatin="gencode.v19.canonical.pext_filtered.chromatin_features.bed.gz"
raw_features_protein="gencode.v19.canonical.pext_filtered.protein_features.bed.gz"
raw_features_constraint="gencode.v19.canonical.pext_filtered.constraint_features.bed.gz"
raw_features_variation="gencode.v19.canonical.pext_filtered.variation_features.bed.gz"
raw_features_merged_no_variation="gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz"
raw_features_merged="gencode.v19.canonical.pext_filtered.all_features.bed.gz"
rCNV_bucket="gs://rcnv_project"

# Copy all fine-mapped gene lists
mkdir finemap_stats/
gsutil -m cp \
  ${rCNV_bucket}/analysis/gene_burden/fine_mapping/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.*.tsv \
  finemap_stats/

# Make input tsvs
for wrapper in 1; do
  echo -e "Prior\tgrey70\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.naive_priors.tsv"
  echo -e "Posterior\t'#264653'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genetics_only.tsv"
  echo -e "Genomic features\t'#490C65'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genomic_features.tsv"
  echo -e "Gene expression\t'#BA7FD0'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.expression_features.tsv"
  echo -e "Chromatin\t'#001588'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.chromatin_features.tsv"
  echo -e "Protein\t'#0180C9'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.protein_features.tsv"
  echo -e "Gene constraint\t'#F6313E'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.constraint_features.tsv"
  echo -e "Variation\t'#FFA300'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.variation_features.tsv"
  echo -e "Full (no var.)\t'#46A040'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.merged_no_variation_features.tsv"
  echo -e "Full model\t'#00441B'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.merged_features.tsv"
done > finemap_roc_input.tsv
for wrapper in 1; do
  echo -e "Genomic features\t'#490C65'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.genomic_features.tsv\t${raw_features_genomic}"
  echo -e "Gene expression\t'#BA7FD0'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.expression_features.tsv\t${raw_features_expression}"
  echo -e "Chromatin\t'#001588'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.chromatin_features.tsv\t${raw_features_chromatin}"
  echo -e "Protein\t'#0180C9'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.protein_features.tsv\t${raw_features_protein}"
  echo -e "Gene constraint\t'#F6313E'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.constraint_features.tsv\t${raw_features_constraint}"
  echo -e "Variation\t'#FFA300'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.variation_features.tsv\t${raw_features_variation}"
  echo -e "Full (no var.)\t'#46A040'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.merged_no_variation_features.tsv\t${raw_features_merged_no_variation}"
  echo -e "Full model\t'#00441B'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.all_genes_from_blocks.merged_features.tsv\t${raw_features_merged}"
done > finemap_feature_cor_input.tsv

# Make all gene truth sets
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/genes/gene_lists ./

# HPO-associated
while read prefix hpo; do
  awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }' \
    gene_lists/$prefix.HPOdb.genes.list
done < ${phenotype_list} \
| sort -Vk1,1 -k2,2V | uniq \
| cat <( echo -e "#HPO\tgene" ) - \
> hpo_truth_set.tsv

# Make CNV type-dependent truth sets
case ${CNV} in
  "DEL")
    # Union (ClinGen HI + DDG2P dominant LoF)
    while read prefix hpo; do
      cat gene_lists/ClinGen.hmc_haploinsufficient.genes.list \
          gene_lists/DDG2P.hmc_lof.genes.list \
      | awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }'
    done < ${phenotype_list} \
    | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > union_truth_set.lof.tsv

    # Intersection (ClinGen HI + DDG2P dominant LoF)
    while read prefix hpo; do
      fgrep -wf \
        gene_lists/ClinGen.hmc_haploinsufficient.genes.list \
        gene_lists/DDG2P.hmc_lof.genes.list \
      | awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }'
    done < ${phenotype_list} \
    | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > intersection_truth_set.lof.tsv

    # ClinGen HI alone
    while read prefix hpo; do
      awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }' \
        gene_lists/ClinGen.all_haploinsufficient.genes.list
    done < ${phenotype_list} \
    | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > clingen_truth_set.lof.tsv

    # DDG2P lof + other alone
    while read prefix hpo; do
      cat gene_lists/DDG2P.all_lof.genes.list \
          gene_lists/DDG2P.all_other.genes.list \
      | awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }'
    done < ${phenotype_list} \
    | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > ddg2p_truth_set.lof.tsv

    # Combine all truth sets into master union
    cat union_truth_set.lof.tsv \
        clingen_truth_set.lof.tsv \
        ddg2p_truth_set.lof.tsv \
        hpo_truth_set.tsv \
    | fgrep -v "#" | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > master_union_truth_set.lof.tsv

    # Write truth set input tsv
    for wrapper in 1; do
      echo -e "Union of all truth sets\tmaster_union_truth_set.lof.tsv"
      echo -e "ClinGen HI & DECIPHER LoF (union)\tunion_truth_set.lof.tsv"
      echo -e "ClinGen HI & DECIPHER LoF (intersection)\tintersection_truth_set.lof.tsv"
      echo -e "ClinGen HI (any confidence)\tclingen_truth_set.lof.tsv"
      echo -e "DECIPHER dominant LoF/unk. (any confidence)\tddg2p_truth_set.lof.tsv"
      echo -e "HPO-matched disease genes\thpo_truth_set.tsv"
    done > finemap_roc_truth_sets.tsv
    ;;

  "DUP")
    # Union (ClinGen HI + DDG2P dominant CG)
    while read prefix hpo; do
      cat gene_lists/ClinGen.hmc_triplosensitive.genes.list \
          gene_lists/DDG2P.hmc_gof.genes.list \
      | awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }'
    done < ${phenotype_list} \
    | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > union_truth_set.gof.tsv

    # ClinGen triplo alone
    while read prefix hpo; do
      awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }' \
        gene_lists/ClinGen.all_triplosensitive.genes.list
    done < ${phenotype_list} \
    | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > clingen_truth_set.triplo.tsv

    # DDG2P gof + other alone
    while read prefix hpo; do
      cat gene_lists/DDG2P.all_gof.genes.list \
          gene_lists/DDG2P.all_other.genes.list \
      | awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }'
    done < ${phenotype_list} \
    | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > ddg2p_truth_set.gof.tsv

    # Combine all truth sets into master union
    cat union_truth_set.gof.tsv \
        clingen_truth_set.gof.tsv \
        ddg2p_truth_set.gof.tsv \
        hpo_truth_set.tsv \
    | fgrep -v "#" | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > master_union_truth_set.gof.tsv

    # Write truth set input tsv
    for wrapper in 1; do
      echo -e "Union of all truth sets\tmaster_union_truth_set.gof.tsv"
      echo -e "ClinGen TS & DECIPHER GoF (union)\tunion_truth_set.gof.tsv"
      echo -e "ClinGen TS (any confidence)\tclingen_truth_set.triplo.tsv"
      echo -e "DECIPHER dominant GoF/unk. (any confidence)\tddg2p_truth_set.gof.tsv"
      echo -e "HPO-matched disease genes\thpo_truth_set.tsv"
    done > finemap_roc_truth_sets.tsv
    ;;

esac

fgrep -wvf \
  gene_lists/HP0000118.HPOdb.genes.list \
  gene_lists/gnomad.v2.1.1.mutation_tolerant.genes.list \
> true_negatives.genes.tsv

# Plot finemapping QC & feature correlations
for subgroup in developmental adult; do
  mkdir ${freq_code}_${CNV}_finemap_plots_$subgroup/
  /opt/rCNV2/analysis/genes/plot_finemap_results.R \
    finemap_roc_input.tsv \
    finemap_roc_truth_sets.tsv \
    true_negatives.genes.tsv \
    $subgroup \
    ${freq_code}_${CNV}_finemap_plots_$subgroup/${freq_code}.${CNV}.$subgroup.finemap_results
  /opt/rCNV2/analysis/genes/plot_finemap_coefficients.R \
    finemap_feature_cor_input.tsv \
    $subgroup \
    ${freq_code}_${CNV}_finemap_plots_$subgroup/${freq_code}.${CNV}.$subgroup.finemap_feature_cors

  # Compress results
  tar -czvf \
    ${freq_code}_${CNV}_finemap_plots_$subgroup.tgz \
    ${freq_code}_${CNV}_finemap_plots_$subgroup
done
