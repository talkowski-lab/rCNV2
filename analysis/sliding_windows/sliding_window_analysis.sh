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
p_cutoff=0.000003715428
meta_p_cutoff=0.000003715428
meta_model_prefix="fe"
bin_overlap=0.5
pad_controls=50000
max_manhattan_phred_p=30
# Test/dev parameters (anxiety, with bad case:control imbalance)
hpo="HP:0100852"
prefix="HP0100852"



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
        --max-phred-p ${max_manhattan_phred_p} \
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
      --max-phred-p ${max_manhattan_phred_p} \
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
p_cutoff=0.000003715428
n_pheno_perms=50
meta_model_prefix="fe"
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
        --spa \
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
p_val_column_name="meta_phred_p"
while read prefix hpo; do
  echo -e "\nSTARTING $prefix\n"
  gsutil -m cp \
    "${rCNV_bucket}/analysis/sliding_windows/$prefix/${freq_code}/permutations/$prefix.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.perm_*.bed.gz" \
    perm_res/
  for i in $( seq 1 ${n_pheno_perms} ); do
    p_idx=$( zcat perm_res/$prefix.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.perm_$i.bed.gz \
             | sed -n '1p' | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' \
             | fgrep -w ${p_val_column_name} | cut -f2 )
    zcat perm_res/$prefix.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.perm_$i.bed.gz \
    | grep -ve '^#' \
    | awk -v p_idx=$p_idx '{ print $(p_idx) }' \
    | cat <( echo "$prefix.${CNV}.$i" ) - \
    > perm_res/$prefix.${freq_code}.${CNV}.sliding_window.meta_analysis.permuted_p_values.$i.txt
  done
  rm perm_res/$prefix.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.perm_*.bed.gz
  echo -e "\nFINISHED $prefix\n"
done < ${phenotype_list}
echo -e "\nMAKING P-VALUE MATRIX\n"
paste perm_res/*.sliding_window.meta_analysis.permuted_p_values.*.txt \
| gzip -c \
> ${freq_code}.${CNV}.permuted_pval_matrix.txt.gz

# Calculate empirical P-value cutoffs
# Genome-wide
fdr_table_suffix="empirical_genome_wide_pval"
fdr_target=${p_cutoff}
echo -e "\nANALYZING P-VALUE MATRIX\n"
/opt/rCNV2/analysis/sliding_windows/calc_empirical_fdr.R \
  --cnv ${CNV} \
  --fdr-target ${fdr_target} \
  --linear-fit \
  --flat-ladder \
  --plot sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}_permutation_results.png \
  ${freq_code}.${CNV}.permuted_pval_matrix.txt.gz \
  ${metacohort_sample_table} \
  sliding_window.${freq_code}.${CNV}.${fdr_table_suffix}




# Test/dev parameters (all cases)
hpo="HP:0000118"
prefix="HP0000118"
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
p_cutoff=0.000003715428
meta_p_cutoff=0.000003715428
meta_model_prefix="fe"
bin_overlap=0.5
pad_controls=50000
# Test/dev parameters (anxiety, with bad case:control imbalance)
hpo="HP:0100852"
prefix="HP0100852"

# Copy necessary data for local testing (without running the above -- this is not in the WDL)
gsutil -m cp \
  ${rCNV_bucket}/analysis/analysis_refs/sliding_window.${freq_code}.*.empirical_genome_wide_pval.hpo_cutoffs.tsv \
  ./
gsutil -m cp \
  ${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/stats/meta**.stats.bed.gz* \
  ./



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
                  sliding_window.${freq_code}.DEL.empirical_genome_wide_pval.hpo_cutoffs.tsv )
  DUP_p_cutoff=$( awk -v hpo=${prefix} '{ if ($1==hpo) print $2 }' \
                  sliding_window.${freq_code}.DEL.empirical_genome_wide_pval.hpo_cutoffs.tsv )

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
      --spa \
      ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.input.txt \
      ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed
    bgzip -f ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed
    tabix -f ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed.gz

    # Generate Manhattan & QQ plots
    /opt/rCNV2/utils/plot_manhattan_qq.R \
      --p-col-name "meta_phred_p" \
      --p-is-phred \
      --max-phred-p ${max_manhattan_phred_p} \
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
    --max-phred-p ${max_manhattan_phred_p} \
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
# Download data
while read prefix hpo; do
  for CNV in DEL DUP; do
    echo -e "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/stats/${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed.gz"
  done
done < ${phenotype_list} \
| gsutil -m cp -I meta_res/
# Primary p-values
p_val_column_name="meta_phred_p"
while read prefix hpo; do
  echo -e "$prefix\n\n"
  for CNV in DEL DUP; do
      stats=meta_res/$prefix.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed.gz
      if [ -e $stats ]; then
        p_idx=$( zcat $stats | sed -n '1p' | sed 's/\t/\n/g' \
                 | awk -v OFS="\t" '{ print $1, NR }' \
                 | fgrep -w ${p_val_column_name} | cut -f2 )
        zcat $stats | grep -ve '^#' \
        | awk -v p_idx=$p_idx '{ print $(p_idx) }' \
        | cat <( echo "$prefix.$CNV" ) - \
        > meta_res/$prefix.${freq_code}.$CNV.sliding_window.meta_analysis.p_values.txt
      fi
  done
done < ${phenotype_list}
paste meta_res/*.${freq_code}.$CNV.sliding_window.meta_analysis.p_values.txt \
| gzip -c \
> ${freq_code}.observed_pval_matrix.txt.gz
# Secondary p-values
p_val_column_name="meta_phred_p_secondary"
while read prefix hpo; do
  echo -e "$prefix\n\n"
  for CNV in DEL DUP; do
      stats=meta_res/$prefix.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed.gz
      if [ -e $stats ]; then
        p_idx=$( zcat $stats | sed -n '1p' | sed 's/\t/\n/g' \
                 | awk -v OFS="\t" '{ print $1, NR }' \
                 | fgrep -w ${p_val_column_name} | cut -f2 )
        zcat $stats | grep -ve '^#' \
        | awk -v p_idx=$p_idx '{ print $(p_idx) }' \
        | cat <( echo "$prefix.$CNV" ) - \
        > meta_res/$prefix.${freq_code}.$CNV.sliding_window.meta_analysis.secondary_p_values.txt
      fi
  done
done < ${phenotype_list}
paste meta_res/*.${freq_code}.$CNV.sliding_window.meta_analysis.secondary_p_values.txt \
| gzip -c \
> ${freq_code}.observed_pval_matrix_secondary.txt.gz
# Make a representative set of bins, for GD comparisons
while read prefix hpo; do
  zcat meta_res/$prefix.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed.gz \
  | cut -f1-3 | bgzip -c \
  > ${freq_code}.sliding_window.meta_analysis.bins.bed.gz
done < <( head -n1 ${phenotype_list} )
# DEV NOTE: these p-values can be visualized with plot_sliding_window_meta_analysis_p_values.R,
#           which is currently just a code snippet referencing local filepaths




# Refine final set of significant segments
# Test/dev parameters
freq_code="rCNV"
CNV="DEL"
phenotype_list="refs/test_phenotypes.list"
metacohort_list="refs/rCNV_metacohort_list.txt"
metacohort_sample_table="refs/HPOs_by_metacohort.table.tsv"
rCNV_bucket="gs://rcnv_project"
meta_p_cutoffs_tsv="refs/sliding_window.rCNV.DEL.empirical_genome_wide_pval.hpo_cutoffs.tsv"
meta_secondary_p_cutoff=0.05
meta_nominal_cohorts_cutoff=2
sig_window_pad=200000
credset=0.99
gtf="gencode.v19.canonical.pext_filtered.gtf.gz"


# Download all meta-analysis stats files and necessary data
mkdir stats/
gsutil -m cp \
  ${rCNV_bucket}/analysis/sliding_windows/**.${freq_code}.**.sliding_window.meta_analysis.stats.bed.gz \
  stats/
mkdir refs/
gsutil -m cp \
  ${rCNV_bucket}/analysis/analysis_refs/* \
  ${rCNV_bucket}/refs/GRCh37.cytobands.bed.gz \
  refs/

# Apply an initial loose mask per HPO to P<0.01 regions ±sig_window_pad 
# to reduce I/O time reading sumstats files in refinement
# Also add all known GD regions for null variance estimation
while read prefix hpo; do
  echo -e "Subsetting $prefix summary stats prior to refinement"
  allstats="stats/$prefix.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.bed.gz"
  primary_p_idx=$( zcat $allstats | sed -n '1p' | sed 's/\t/\n/g' \
                   | awk -v FS="\t" '{ if ($1=="meta_phred_p") print NR }' )
  zcat $allstats \
  | grep -ve '^#' \
  | awk -v FS="\t" -v OFS="\t" -v idx=$primary_p_idx \
    '{ if ($(idx) >= 2 && $(idx) != "NA") print $1, $2, $3 }' \
  | cat - <( zcat /opt/rCNV2/refs/UKBB_GD.Owen_2018.bed.gz \
             | grep -ve '^#' | cut -f1-3 ) \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  | awk -v OFS="\t" -v dist=${sig_window_pad} \
    '{ if ($2-dist < 1) print $1, "1", $3+dist; else print $1, $2-dist, $3+dist }' \
  | bedtools intersect -wa -u -header \
    -a $allstats \
    -b - \
  | bgzip -c \
  > stats/$prefix.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.subset.bed.gz
  tabix -f stats/$prefix.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.subset.bed.gz
done < ${phenotype_list}

# Write tsv inputs
while read prefix hpo; do
  for wrapper in 1; do
    echo "$hpo"
    echo "stats/$prefix.${freq_code}.${CNV}.sliding_window.meta_analysis.stats.subset.bed.gz"
    awk -v x=$prefix -v FS="\t" '{ if ($1==x) print $2 }' \
      ${meta_p_cutoffs_tsv}
  done | paste -s
done < ${phenotype_list} \
> ${freq_code}.${CNV}.segment_refinement.stats_input.tsv
echo "/opt/rCNV2/refs/UKBB_GD.Owen_2018.${CNV}.bed.gz" \
> known_causal_loci_lists.${CNV}.tsv

# Refine significant segments
/opt/rCNV2/analysis/sliding_windows/refine_significant_regions.py \
  --cnv ${CNV} \
  --secondary-p-cutoff ${meta_secondary_p_cutoff} \
  --min-nominal ${meta_nominal_cohorts_cutoff} \
  --secondary-or-nominal \
  --credible-sets ${credset} \
  --distance ${sig_window_pad} \
  --known-causal-loci-list known_causal_loci_lists.${CNV}.tsv \
  --cytobands refs/GRCh37.cytobands.bed.gz \
  --sig-loci-bed ${freq_code}.${CNV}.final_segments.loci.pregenes.bed \
  --sig-assoc-bed ${freq_code}.${CNV}.final_segments.associations.pregenes.bed \
  ${freq_code}.${CNV}.segment_refinement.stats_input.tsv \
  ${metacohort_sample_table}

# Annotate final regions with genes & sort by coordinates
for entity in loci associations; do
  /opt/rCNV2/analysis/sliding_windows/get_genes_per_region.py \
    -o ${freq_code}.${CNV}.final_segments.$entity.bed \
    ${freq_code}.${CNV}.final_segments.$entity.pregenes.bed \
    ${gtf}
done



# Plot summary figures for final regions
/opt/rCNV2/analysis/sliding_windows/regions_summary.plot.R \
  -o "${freq_code}.final_segments." \
  ${freq_code}.DEL.final_segments.loci.bed \
  ${freq_code}.DUP.final_segments.loci.bed






