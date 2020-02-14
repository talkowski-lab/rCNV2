#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control CNV burdens in sliding windows, genome-wide


# Launch docker image
docker run --rm -it talkowski/rcnv


# Copy all filtered CNV data, sliding windows, and other references 
# from the project Google Bucket (note: requires permissions)
gcloud auth login
mkdir cleaned_cnv/
gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/* cleaned_cnv/
mkdir windows/
gsutil -m cp -r gs://rcnv_project/cleaned_data/binned_genome/* windows/
mkdir refs/
gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/


# Test/dev parameters (seizures)
hpo="HP:0001250"
prefix="HP0001250"
meta="meta1"
freq_code="rCNV"
phenotype_list="refs/test_phenotypes.list"
metacohort_list="refs/rCNV_metacohort_list.txt"
metacohort_sample_table="refs/HPOs_by_metacohort.table.tsv"
binned_genome="windows/GRCh37.200kb_bins_10kb_steps.raw.bed.gz"
rCNV_bucket="gs://rcnv_project"
p_cutoff=0.00000385862
meta_p_cutoff=0.00000385862
meta_model_prefix="re"
bin_overlap=0.5
pad_controls=50000




# Count CNVs in cases and controls per phenotype, split by metacohort and CNV type
# Iterate over phenotypes
while read prefix hpo; do
  descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" '{ print $2 }' )
  # Iterate over metacohorts
  while read meta cohorts; do

    # Set metacohort parameters
    cnv_bed="cleaned_cnv/$meta.${freq_code}.bed.gz"
    descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
               | awk -v FS="\t" '{ print $2 }' )
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
      # Set CNV-specific parameters
      case "$CNV" in
        DEL)
          highlight_bed=/opt/rCNV2/refs/UKBB_GD.Owen_2018.DEL.bed.gz
          highlight_title="Known DEL GDs (Owen 2018)"
          ;;
        DUP)
          highlight_bed=/opt/rCNV2/refs/UKBB_GD.Owen_2018.DUP.bed.gz
          highlight_title="Known DUP GDs (Owen 2018)"
          ;;
      esac

      # Count CNVs
      /opt/rCNV2/analysis/sliding_windows/count_cnvs_per_window.py \
        --fraction ${bin_overlap} \
        --pad-controls ${pad_controls} \
        -t $CNV \
        --hpo ${hpo} \
        -z \
        -o "$meta.${prefix}.${freq_code}.$CNV.sliding_window.counts.bed.gz" \
        $cnv_bed \
        ${binned_genome}
      tabix -f "$meta.${prefix}.${freq_code}.$CNV.sliding_window.counts.bed.gz"

      # Perform burden test
      /opt/rCNV2/analysis/sliding_windows/window_burden_test.R \
        --pheno-table ${metacohort_sample_table} \
        --cohort-name $meta \
        --case-hpo ${hpo} \
        --bgzip \
        "$meta.${prefix}.${freq_code}.$CNV.sliding_window.counts.bed.gz" \
        "$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz"
      tabix -f "$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz"

      # Generate Manhattan & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --p-col-name "fisher_phred_p" \
        --p-is-phred \
        --max-phred-p 100 \
        --cutoff ${p_cutoff} \
        --highlight-bed "$highlight_bed" \
        --highlight-name "$highlight_title" \
        --label-prefix "$CNV" \
        --title "$title" \
        "$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz" \
        "$meta.${prefix}.${freq_code}.$CNV.sliding_window"
    done

    # Generate Miami & QQ plots
    /opt/rCNV2/utils/plot_manhattan_qq.R \
      --miami \
      --p-col-name "fisher_phred_p" \
      --p-is-phred \
      --max-phred-p 100 \
      --cutoff ${p_cutoff} \
      --highlight-bed /opt/rCNV2/refs/UKBB_GD.Owen_2018.DUP.bed.gz \
      --highlight-name "Known DUP GDs (Owen 2018)" \
      --label-prefix "DUP" \
      --highlight-bed-2 /opt/rCNV2/refs/UKBB_GD.Owen_2018.DEL.bed.gz \
      --highlight-name-2 "Known DEL GDs (Owen 2018)" \
      --label-prefix-2 "DEL" \
      --title "$title" \
      "$meta.${prefix}.${freq_code}.DUP.sliding_window.stats.bed.gz" \
      "$meta.${prefix}.${freq_code}.DEL.sliding_window.stats.bed.gz" \
      "$meta.${prefix}.${freq_code}.sliding_window"
  done < ${metacohort_list}
done < refs/test_phenotypes.list




# Run phenotype permutation to determine empirical FDR cutoff
# Test/dev parameters
hpo="HP:0001250"
prefix="HP0001250"
meta="meta1"
freq_code="rCNV"
phenotype_list="refs/test_phenotypes.list"
metacohort_list="refs/rCNV_metacohort_list.txt"
metacohort_sample_table="refs/HPOs_by_metacohort.table.tsv"
binned_genome="windows/GRCh37.200kb_bins_10kb_steps.raw.bed.gz"
rCNV_bucket="gs://rcnv_project"
p_cutoff=0.00000385862
n_pheno_perms=20
meta_model_prefix="re"
i=1
bin_overlap=0.5
pad_controls=50000

# Copy all filtered CNV data, sliding windows, and other references 
# from the project Google Bucket (note: requires permissions)
gcloud auth login
mkdir cleaned_cnv/
gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/* cleaned_cnv/
mkdir windows/
gsutil -m cp -r gs://rcnv_project/cleaned_data/binned_genome/* windows/
mkdir refs/
gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/

# Count CNVs in cases and controls per phenotype, split by metacohort and CNV type
# Iterate over phenotypes
while read prefix hpo; do
  for i in $( seq 1 ${n_pheno_perms} ); do

    # Shuffle phenotypes for each metacohort CNV dataset, and restrict CNVs from
    # phenotype of relevance
    if ! [ -e shuffled_cnv/ ]; then
      mkdir shuffled_cnv/
    fi
    yes $i | head -n1000000 > seed_$i.txt
    while read meta cohorts; do

      cnvbed="cleaned_cnv/$meta.${freq_code}.bed.gz"

      # Determine CNV size quintiles
      for CNV in DEL DUP; do
        zcat $cnvbed \
        | awk -v CNV="$CNV" '{ if ($1 !~ "#" && $5==CNV) print $3-$2 }' \
        | /opt/rCNV2/utils/quantiles.py \
          --quantiles "0,0.2,0.4,0.6,0.8,1.0" \
          --no-header \
        > $meta.$CNV.quantiles.tsv
      done

      # Shuffle phenotypes, matching by CNV type and size quintile
      for CNV in DEL DUP; do
        for qr in $( seq 1 5 ); do
          smin=$( awk -v qr="$qr" '{ if (NR==qr) print $2 }' $meta.$CNV.quantiles.tsv )
          smax=$( awk -v qr="$qr" '{ if (NR==(qr+1)) print $2 + 1 }' $meta.$CNV.quantiles.tsv )
          zcat $cnvbed | sed '1d' \
          | awk -v smin="$smin" -v smax="$smax" -v CNV="$CNV" \
            '{ if ($3-$2>=smin && $3-$2<smax && $5==CNV) print $0 }' \
          > cnv_subset.bed
          paste <( cut -f1-5 cnv_subset.bed ) \
                <( cut -f6 cnv_subset.bed | shuf --random-source seed_$i.txt ) \
          | awk -v hpo=${hpo} '{ if ($NF ~ "HEALTHY_CONTROL" || $NF ~ hpo) print $0 }'
          rm cnv_subset.bed
        done
      done \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | cat <( tabix -H $cnvbed ) - \
      | bgzip -c \
      > shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed.gz
      tabix -f shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed.gz

    done < ${metacohort_list}

    # Iterate over CNV types
    for CNV in DEL DUP; do
      # Perform association test for each metacohort
      while read meta cohorts; do

        echo -e "[$( date )] Starting permutation $i for $CNV in ${hpo} from $meta...\n"

        cnv_bed="shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed.gz"

        # Count CNVs
        /opt/rCNV2/analysis/sliding_windows/count_cnvs_per_window.py \
          --fraction ${bin_overlap} \
          --pad-controls ${pad_controls} \
          -t $CNV \
          --hpo ${hpo} \
          -z \
          -o "$meta.${prefix}.${freq_code}.$CNV.sliding_window.counts.bed.gz" \
          $cnv_bed \
          ${binned_genome}

        # Perform burden test
        /opt/rCNV2/analysis/sliding_windows/window_burden_test.R \
          --pheno-table ${metacohort_sample_table} \
          --cohort-name $meta \
          --case-hpo ${hpo} \
          --bgzip \
          "$meta.${prefix}.${freq_code}.$CNV.sliding_window.counts.bed.gz" \
          "$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz"
          tabix -f "$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz"

      done < <( fgrep -v "mega" ${metacohort_list} )

      # Perform meta-analysis
      while read meta cohorts; do
        echo -e "$meta\t$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz"
      done < <( fgrep -v mega ${metacohort_list} ) \
      > ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.input.txt
      /opt/rCNV2/analysis/sliding_windows/window_meta_analysis.R \
        --model ${meta_model_prefix} \
        --p-is-phred \
        ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.input.txt \
        ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.perm_$i.bed
      bgzip -f ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.perm_$i.bed
      tabix -f ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.perm_$i.bed.gz

      # Copy results to output bucket
      gsutil cp \
        ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.perm_$i.bed.gz \
        "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/permutations/"
    done
  done
done < refs/test_phenotypes.list

# Gather all permutation results and compute FDR CDFs
mkdir perm_res/
while read prefix hpo; do
  echo -e "$prefix\n\n"
  gsutil -m cp \
    "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/permutations/**.stats.perm_*.bed.gz" \
    perm_res/
  for CNV in DEL DUP; do
    for i in $( seq 1 ${n_pheno_perms} ); do
      if [ -e perm_res/$prefix.${freq_code}.$CNV.sliding_window.meta_analysis.stats.perm_$i.bed.gz ]; then
        zcat perm_res/$prefix.${freq_code}.$CNV.sliding_window.meta_analysis.stats.perm_$i.bed.gz \
        | grep -ve '^#' \
        | awk '{ print $NF }' \
        | cat <( echo "$prefix.$CNV.$i" ) - \
        > perm_res/$prefix.${freq_code}.$CNV.sliding_window.meta_analysis.permuted_p_values.$i.txt
      fi
    done
  done
done < ${phenotype_list}
paste perm_res/*.sliding_window.meta_analysis.permuted_p_values.*.txt \
| gzip -c \
> ${freq_code}.permuted_pval_matrix.txt.gz

# Calculate empirical P-value cutoffs
# Genome-wide
fdr_table_suffix="empirical_genome_wide_pval"
fdr_target=${p_cutoff}
for CNV in DEL DUP; do
  /opt/rCNV2/analysis/sliding_windows/calc_empirical_genome_wide.R \
    --cnv ${CNV} \
    --fdr-target ${fdr_target} \
    --plot ${freq_code}.${CNV}.${fdr_table_suffix}_permutation_results.png \
    ${freq_code}.${CNV}.permuted_pval_matrix.txt.gz \
    ${metacohort_sample_table} \
    sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}
done
  
# 1% FDR
fdr_table_suffix="empirical_fdr_1pct_pval"
fdr_target=0.01
for CNV in DEL DUP; do
  /opt/rCNV2/analysis/sliding_windows/calc_empirical_genome_wide.R \
    --cnv ${CNV} \
    --fdr-target ${fdr_target} \
    --plot ${freq_code}.${CNV}.${fdr_table_suffix}_permutation_results.png \
    ${freq_code}.${CNV}.permuted_pval_matrix.txt.gz \
    ${metacohort_sample_table} \
    sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}
done
  
# 5% FDR
fdr_table_suffix="empirical_genome_wide_pval"
fdr_target=0.05
for CNV in DEL DUP; do
  /opt/rCNV2/analysis/sliding_windows/calc_empirical_genome_wide.R \
    --cnv ${CNV} \
    --fdr-target ${fdr_target} \
    --plot ${freq_code}.${CNV}.${fdr_table_suffix}_permutation_results.png \
    ${freq_code}.${CNV}.permuted_pval_matrix.txt.gz \
    ${metacohort_sample_table} \
    sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}
done




# Run meta-analysis for each phenotype
while read prefix hpo; do

  # Get metadata for meta-analysis
  mega_idx=$( head -n1 "${metacohort_sample_table}" \
              | sed 's/\t/\n/g' \
              | awk '{ if ($1=="mega") print NR }' )
  ncase=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
           | awk -v FS="\t" -v mega_idx="$mega_idx" '{ print $mega_idx }' \
           | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
  nctrl=$( fgrep -w "HEALTHY_CONTROL" "${metacohort_sample_table}" \
           | awk -v FS="\t" -v mega_idx="$mega_idx" '{ print $mega_idx }' \
           | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
  descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" '{ print $2 }' )
  title="$descrip (${hpo})\nMeta-analysis of $ncase cases and $nctrl controls"
  DEL_p_cutoff=$( awk -v hpo=${prefix} '{ if ($1==hpo) print $2 }' \
                  sliding_window.${freq_code}.DEL.empirical_genome_wide.hpo_cutoffs.tsv )
  DUP_p_cutoff=$( awk -v hpo=${prefix} '{ if ($1==hpo) print $2 }' \
                  sliding_window.${freq_code}.DEL.empirical_genome_wide.hpo_cutoffs.tsv )

  # Run meta-analysis for each CNV type
  for CNV in DEL DUP; do
    # Set CNV-specific parameters
    case "$CNV" in
      DEL)
        highlight_bed=/opt/rCNV2/refs/UKBB_GD.Owen_2018.DEL.bed.gz
        highlight_title="Known DEL GDs (Owen 2018)"
        meta_p_cutoff=$DEL_p_cutoff
        ;;
      DUP)
        highlight_bed=/opt/rCNV2/refs/UKBB_GD.Owen_2018.DUP.bed.gz
        highlight_title="Known DUP GDs (Owen 2018)"
        meta_p_cutoff=$DUP_p_cutoff
        ;;
    esac

    # Perform meta-analysis
    while read meta cohorts; do
      echo -e "$meta\t$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz"
    done < <( fgrep -v mega ${metacohort_list} ) \
    > ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.input.txt
    /opt/rCNV2/analysis/sliding_windows/window_meta_analysis.R \
      --or-corplot ${prefix}.${freq_code}.$CNV.sliding_window.or_corplot_grid.jpg \
      --model ${meta_model_prefix} \
      --p-is-phred \
      ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.input.txt \
      ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed
    bgzip -f ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed
    tabix -f ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed.gz

    # Generate Manhattan & QQ plots
    /opt/rCNV2/utils/plot_manhattan_qq.R \
      --p-col-name "meta_phred_p" \
      --p-is-phred \
      --cutoff $meta_p_cutoff \
      --highlight-bed "$highlight_bed" \
      --highlight-name "$highlight_title" \
      --label-prefix "$CNV" \
      --title "$title" \
      "${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed.gz" \
      "${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis"
  done

  # Generate Miami & QQ plots
  /opt/rCNV2/utils/plot_manhattan_qq.R \
    --miami \
    --p-col-name "meta_phred_p" \
    --p-is-phred \
    --cutoff $DUP_p_cutoff \
    --highlight-bed /opt/rCNV2/refs/UKBB_GD.Owen_2018.DUP.bed.gz \
    --highlight-name "Known DUP GDs (Owen 2018)" \
    --label-prefix "DUP" \
    --cutoff-2 $DEL_p_cutoff \
    --highlight-bed-2 /opt/rCNV2/refs/UKBB_GD.Owen_2018.DEL.bed.gz \
    --highlight-name-2 "Known DEL GDs (Owen 2018)" \
    --label-prefix-2 "DEL" \
    --title "$title" \
    "${prefix}.${freq_code}.DUP.sliding_window.meta_analysis.stats.bed.gz" \
    "${prefix}.${freq_code}.DEL.sliding_window.meta_analysis.stats.bed.gz" \
    "${prefix}.${freq_code}.sliding_window.meta_analysis"

  # Copy results to output bucket
  gsutil -m cp *.sliding_window.meta_analysis.stats.bed.gz* \
    "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/stats/"
  gsutil -m cp *.sliding_window.or_corplot_grid.jpg \
    "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/plots/"
  gsutil -m cp *.sliding_window.meta_analysis.*.png \
    "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/plots/"
done < refs/test_phenotypes.list


# Collapse all meta-analysis p-values into single matrix for visualizing calibration
mkdir meta_res/
while read prefix hpo; do
  echo -e "$prefix\n\n"
  gsutil -m cp \
    "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/stats/${prefix}.${freq_code}.*.sliding_window.meta_analysis.stats.bed.gz" \
    meta_res/
  for CNV in DEL DUP; do
      if [ -e meta_res/$prefix.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed.gz ]; then
        zcat meta_res/$prefix.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed.gz \
        | grep -ve '^#' \
        | awk '{ print $NF }' \
        | cat <( echo "$prefix.$CNV" ) - \
        > meta_res/$prefix.${freq_code}.$CNV.sliding_window.meta_analysis.p_values.txt
      fi
  done
done < ${phenotype_list}
paste meta_res/*.${freq_code}.$CNV.sliding_window.meta_analysis.p_values.txt \
| gzip -c \
> ${freq_code}.observed_pval_matrix.txt.gz




# Refine final set of significant segments
# Test/dev parameters
freq_code="rCNV"
CNV="DEL"
phenotype_list="refs/test_phenotypes.list"
metacohort_list="refs/rCNV_metacohort_list.txt"
binned_genome="windows/GRCh37.200kb_bins_10kb_steps.raw.bed.gz"
rCNV_bucket="gs://rcnv_project"
meta_p_cutoff=0.00000385862
meta_or_cutoff=2
meta_nominal_cohorts_cutoff=2
meta_model_prefix="re"
sig_window_pad=1000000
refine_max_cnv_size=3000000

# Download all meta-analysis stats files and necessary data
mkdir cleaned_cnv/
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/cnv/* cleaned_cnv/
mkdir stats/
gsutil -m cp \
  ${rCNV_bucket}/analysis/sliding_windows/**.${freq_code}.**.sliding_window.meta_analysis.stats.bed.gz \
  stats/
mkdir refs/
gsutil -m cp ${rCNV_bucket}/analysis/analysis_refs/* refs/
mkdir phenos/
gsutil -m cp ${rCNV_bucket}/cleaned_data/phenotypes/filtered/* phenos/

# Iterate over phenotypes and make matrix of p-values, odds ratios (lower 95% CI), and nominal sig cohorts
mkdir pvals/
mkdir ors/
mkdir nomsig/
while read pheno hpo; do
  zcat stats/$pheno.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.bed.gz \
  | awk -v FS="\t" '{ if ($1 !~ "#") print $NF }' \
  | cat <( echo "$pheno.${CNV}" ) - \
  > pvals/$pheno.${CNV}.pvals.txt
  zcat stats/$pheno.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.bed.gz \
  | awk -v FS="\t" '{ if ($1 !~ "#") print $6 }' \
  | cat <( echo "$pheno.${CNV}" ) - \
  > ors/$pheno.${CNV}.lnOR_lower.txt
  zcat stats/$pheno.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.bed.gz \
  | awk -v FS="\t" '{ if ($1 !~ "#") print $4 }' \
  | cat <( echo "$pheno.${CNV}" ) - \
  > nomsig/$pheno.${CNV}.nomsig_counts.txt
done < ${phenotype_list}
paste <( zcat ${binned_genome} | cut -f1-3 ) \
      pvals/*.${CNV}.pvals.txt \
| bgzip -c \
> ${CNV}.pval_matrix.bed.gz
paste <( zcat ${binned_genome} | cut -f1-3 ) \
      ors/*.${CNV}.lnOR_lower.txt \
| bgzip -c \
> ${CNV}.lnOR_lower_matrix.bed.gz
paste <( zcat ${binned_genome} | cut -f1-3 ) \
      nomsig/*.${CNV}.nomsig_counts.txt \
| bgzip -c \
> ${CNV}.nominal_cohort_counts.bed.gz

# Get matrix of window significance labels
/opt/rCNV2/analysis/sliding_windows/get_significant_windows.R \
  --pvalues ${CNV}.pval_matrix.bed.gz \
  --p-is-phred \
  --p-cutoffs sliding_window.${freq_code}.${CNV}.empirical_genome_wide.hpo_cutoffs.tsv \
  --odds-ratios ${CNV}.lnOR_lower_matrix.bed.gz \
  --or-is-ln \
  --min-or ${meta_or_cutoff} \
  --nominal-counts ${CNV}.nominal_cohort_counts.bed.gz \
  --min-nominal ${meta_nominal_cohorts_cutoff} \
  --out-prefix ${freq_code}.${CNV}. \
  ${binned_genome}
bgzip -f ${freq_code}.${CNV}.all_windows_labeled.bed
bgzip -f ${freq_code}.${CNV}.significant_windows.bed

# Define regions to be refined (sig windows padded by $sig_window_pad and merged)
zcat ${freq_code}.${CNV}.significant_windows.bed.gz \
| fgrep -v "#" \
| awk -v buf=${sig_window_pad} -v OFS="\t" '{ print $1, $2-buf, $3+buf }' \
| awk -v OFS="\t" '{ if ($2<0) $2=0; print $1, $2, $3 }' \
| sort -Vk1,1 -k2,2V -k3,3V \
| bedtools merge -i - \
| bgzip -c \
> ${freq_code}.${CNV}.sig_regions_to_refine.bed.gz

# Prep input file
while read meta; do
  echo -e "$meta\tcleaned_cnv/$meta.${freq_code}.bed.gz\tphenos/$meta.cleaned_phenos.txt"
done < <( cut -f1 ${metacohort_list} | fgrep -v "mega" )\
> window_refinement.${freq_code}_metacohort_info.tsv


# Refine associations within regions from above
for CNV in DEL DUP; do
  for contig in $( seq 1 22 ); do
    /opt/rCNV2/analysis/sliding_windows/refine_significant_regions.py \
      --cnv-type ${CNV} \
      --model ${meta_model_prefix} \
      --hpo-p-cutoffs sliding_window.${freq_code}.${CNV}.empirical_genome_wide.hpo_cutoffs.tsv \
      --p-cutoff-ladder sliding_window.${freq_code}.${CNV}.empirical_genome_wide.ncase_cutoff_ladder.tsv \
      --p-is-phred \
      --min-or-lower ${meta_or_cutoff} \
      --retest-min-or-lower ${meta_or_cutoff} \
      --max-cnv-size ${refine_max_cnv_size} \
      --min-nominal ${meta_nominal_cohorts_cutoff} \
      --credible-interval ${credible_interval} \
      --prefix "${freq_code}_${CNV}" \
      --log ${freq_code}.${CNV}.region_refinement.${contig}.log \
      regions_to_refine.bed.gz \
      ${metacohort_info_tsv} \
      pval_matrix.bed.gz \
      labeled_windows.bed.gz \
      ${freq_code}.${CNV}.final_regions.associations.${contig}.bed \
      ${freq_code}.${CNV}.final_regions.loci.${contig}.bed
    bgzip -f ${freq_code}.$CNV.final_regions.associations.${contig}.bed
    bgzip -f ${freq_code}.$CNV.final_regions.loci.${contig}.bed
  done
done

# Annotate final regions with genes
gsutil -m cp -r gs://rcnv_project/cleaned_data/genes ./
for CNV in DEL DUP; do
  /opt/rCNV2/analysis/sliding_windows/get_genes_per_region.py \
    -o ${freq_code}.$CNV.final_regions.loci.bed \
    ${freq_code}.$CNV.final_regions.loci.bed.gz \
    genes/gencode.v19.canonical.gtf.gz
    bgzip -f ${freq_code}.$CNV.final_regions.loci.bed
done

# Plot summary figures for final regions
/opt/rCNV2/analysis/sliding_windows/regions_summary.plot.R \
  -o "${freq_code}.final_regions." \
  ${freq_code}.DEL.final_regions.loci.bed.gz \
  ${freq_code}.DUP.final_regions.loci.bed.gz






