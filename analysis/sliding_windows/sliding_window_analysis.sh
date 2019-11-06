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
p_cutoff=0.0000771724
meta_p_cutoff=0.0000192931
bin_overlap=0.5
pad_controls=50000


# Count CNVs in cases and controls per phenotype, split by metacohort and CNV type
# Iterate over phenotypes
while read prefix hpo; do
  descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" '{ print $2 }' )
  # Iterate over metacohorts
  while read meta cohorts; do
    cnv_bed="cleaned_cnv/$meta.$freq_code.bed.gz"
    meta_idx=$( head -n1 "${metacohort_sample_table}" \
                | sed 's/\t/\n/g' \
                | awk -v meta="$meta" '{ if ($1==meta) print NR }' )
    ncase=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" -v meta_idx="$meta_idx" '{ print $meta_idx }' \
             | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
    nctrl=$( fgrep -w "HEALTHY_CONTROL" "${metacohort_sample_table}" \
             | awk -v FS="\t" -v meta_idx="$meta_idx" '{ print $meta_idx }' \
             | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
    title="$descrip (${hpo})\n$ncase cases vs $nctrl controls in '${meta}' cohort"

    # Iterate over CNV types
    for CNV in CNV DEL DUP; do
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
        *)
          highlight_bed=/opt/rCNV2/refs/UKBB_GD.Owen_2018.bed.gz
          highlight_title="Known GDs (Owen 2018)"
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
        ${cnv_bed} \
        ${binned_genome}

      # Perform burden test
      /opt/rCNV2/analysis/sliding_windows/window_burden_test.R \
        --pheno-table refs/HPOs_by_metacohort.table.tsv \
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
        --cutoff ${p_cutoff} \
        --highlight-bed "$highlight_bed" \
        --highlight-name "$highlight_title" \
        --label-prefix "$CNV" \
        --title "$title" \
        "$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz" \
        "$meta.${prefix}.${freq_code}.$CNV.sliding_window"
    done

    # Generate Miami & QQ plots
    title="$descrip (${hpo})\nCohort '${meta}': $ncase cases vs $nctrl controls"
    /opt/rCNV2/utils/plot_manhattan_qq.R \
      --miami \
      --p-col-name "fisher_phred_p" \
      --p-is-phred \
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




# Run meta-analysis for each CNV type
while read prefix hpo; do
  descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" '{ print $2 }' )
  mega_idx=$( head -n1 "${metacohort_sample_table}" \
              | sed 's/\t/\n/g' \
              | awk '{ if ($1=="mega") print NR }' )
  ncase=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
           | awk -v FS="\t" -v mega_idx="$mega_idx" '{ print $mega_idx }' \
           | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
  nctrl=$( fgrep -w "HEALTHY_CONTROL" "${metacohort_sample_table}" \
           | awk -v FS="\t" -v mega_idx="$mega_idx" '{ print $mega_idx }' \
           | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
  title="$descrip (${hpo})\nMeta-analysis of $ncase cases and $nctrl controls"
  for CNV in DEL DUP CNV; do
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
      *)
        highlight_bed=/opt/rCNV2/refs/UKBB_GD.Owen_2018.bed.gz
        highlight_title="Known GDs (Owen 2018)"
        ;;
    esac

    # Perform meta-analysis
    while read meta cohorts; do
      echo -e "$meta\t$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz"
    done < <( fgrep -v mega ${metacohort_list} ) \
    > ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.input.txt
    /opt/rCNV2/analysis/sliding_windows/window_meta_analysis.R \
      --or-corplot ${prefix}.${freq_code}.$CNV.sliding_window.or_corplot_grid.jpg \
      --model mh \
      ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.input.txt \
      ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed
    bgzip -f ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed
    tabix -f ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed.gz

    # Generate Manhattan & QQ plots
    /opt/rCNV2/utils/plot_manhattan_qq.R \
      --p-col-name "meta_phred_p" \
      --p-is-phred \
      --cutoff ${meta_p_cutoff} \
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
    --cutoff ${meta_p_cutoff} \
    --highlight-bed /opt/rCNV2/refs/UKBB_GD.Owen_2018.DUP.bed.gz \
    --highlight-name "Known DUP GDs (Owen 2018)" \
    --label-prefix "DUP" \
    --highlight-bed-2 /opt/rCNV2/refs/UKBB_GD.Owen_2018.DEL.bed.gz \
    --highlight-name-2 "Known DEL GDs (Owen 2018)" \
    --label-prefix-2 "DEL" \
    --title "$title" \
    "${prefix}.${freq_code}.DUP.sliding_window.meta_analysis.stats.bed.gz" \
    "${prefix}.${freq_code}.DEL.sliding_window.meta_analysis.stats.bed.gz" \
    "${prefix}.${freq_code}.sliding_window.meta_analysis"
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
p_cutoff=0.0000771724
n_pheno_perms=10
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

  mkdir shuffled_cnv/

  for i in $( seq 1 ${n_pheno_perms} ); do

    # Shuffle phenotypes for each metacohort CNV dataset, and restrict CNVs from
    # phenotype of relevance
    yes $i | head -n1000000 > seed_$i.txt
    while read meta cohorts; do
      cnvbed="cleaned_cnv/$meta.${freq_code}.bed.gz"
      tabix -H $cnvbed > shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed
      paste <( zcat $cnvbed | sed '1d' | cut -f1-5 ) \
            <( zcat $cnvbed | sed '1d' | cut -f6 \
               | shuf --random-source seed_$i.txt ) \
      | awk -v hpo=${hpo} '{ if ($NF ~ "HEALTHY_CONTROL" || $NF ~ hpo) print $0 }' \
      >> shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed
      bgzip -f shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed
      tabix -f shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed.gz
    done < ${metacohort_list}

    # Iterate over metacohorts & CNV types
    while read meta cohorts; do
      for CNV in DEL DUP CNV; do

        echo -e "[$( date )] Starting permutation $i for $CNV in $hpo from $meta...\n"

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
          --pheno-table refs/HPOs_by_metacohort.table.tsv \
          --cohort-name $meta \
          --case-hpo ${hpo} \
          --bgzip \
          "$meta.${prefix}.${freq_code}.$CNV.sliding_window.counts.bed.gz" \
          "$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz"
          tabix -f "$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz"

        # Perform meta-analysis
        while read meta cohorts; do
          echo -e "$meta\t$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz"
        done < <( fgrep -v mega ${metacohort_list} ) \
        > ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.input.txt
        /opt/rCNV2/analysis/sliding_windows/window_meta_analysis.R \
          --or-corplot ${prefix}.${freq_code}.$CNV.sliding_window.or_corplot_grid.jpg \
          --model mh \
          ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.input.txt \
          ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.perm_$i.bed
        bgzip -f ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.perm_$i.bed
        tabix -f ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.perm_$i.bed.gz

        # Copy results to output bucket
        gsutil cp \
          ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.perm_$i.bed.gz \
          "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/permutations/"
      done
    done < ${metacohort_list}
  done

done < refs/test_phenotypes.list




# Refine final set of significant segments
# Test/dev parameters (seizures)
freq_code="rCNV"
CNV="DEL"
phenotype_list="refs/test_phenotypes.list"
metacohort_list="refs/rCNV_metacohort_list.txt"
metacohort_sample_table="refs/HPOs_by_metacohort.table.tsv"
binned_genome="windows/GRCh37.200kb_bins_10kb_steps.raw.bed.gz"
rCNV_bucket="gs://rcnv_project"
meta_p_cutoff=0.0000192931
sig_window_pad=1000000

# Download all meta-analysis stats files and necessary data
mkdir cleaned_cnv/
gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/* cleaned_cnv/
mkdir stats/
gsutil -m cp \
  ${rCNV_bucket}/analysis/sliding_windows/**.${freq_code}.**.sliding_window.meta_analysis.stats.bed.gz \
  stats/
mkdir windows/
gsutil -m cp -r gs://rcnv_project/cleaned_data/binned_genome/* windows/
mkdir refs/
gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/

# Iterate over phenotypes and make matrix of p-values and odds ratios (lower 95% CI)
mkdir pvals/
mkdir ors/
while read pheno hpo; do
  for CNV in DEL DUP; do
    zcat stats/$pheno.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed.gz \
    | awk -v FS="\t" '{ if ($1 !~ "#") print $NF }' \
    | cat <( echo "$pheno.$CNV" ) - \
    > pvals/$pheno.$CNV.pvals.txt
    zcat stats/$pheno.${freq_code}.$CNV.sliding_window.meta_analysis.stats.bed.gz \
    | awk -v FS="\t" '{ if ($1 !~ "#") print $5 }' \
    | cat <( echo "$pheno.$CNV" ) - \
    > pvals/$pheno.$CNV.or_lower.txt
  done
done < ${phenotype_list}
for CNV in DEL DUP; do
  paste <( zcat windows/GRCh37.200kb_bins_10kb_steps.raw.bed.gz \
           | cut -f1-3 ) \
        pvals/*.$CNV.pvals.txt \
  | bgzip -c \
  > $CNV.pval_matrix.bed.gz
  paste <( zcat windows/GRCh37.200kb_bins_10kb_steps.raw.bed.gz \
           | cut -f1-3 ) \
        pvals/*.$CNV.or_lower.txt \
  | bgzip -c \
  > $CNV.or_lower_matrix.bed.gz
done

# Iterate over phenotypes and count number of individual cohorts with nominal 
# significance per window
mkdir nomsig/
while read pheno hpo; do
  
done < ${phenotype_list}

# Identify regions to be refined (get significant windows, and pad by $sig_window_pad)
for CNV in DEL DUP; do
  /opt/rCNV2/analysis/sliding_windows/get_significant_windows.R \
    --p-is-phred \
    --cutoff ${meta_p_cutoff} \
    --pad-windows ${sig_window_pad} \
    $CNV.pval_matrix.bed.gz \
  | bedtools merge -i - \
  | bgzip -c \
  > $CNV.sig_regions_to_refine.bed.gz
done

# Refine associations within regions from above
while read meta; do
  echo -e "$meta\tcleaned_cnv/$meta.${freq_code}.bed.gz"
done < <( cut -f1 ${metacohort_list} | fgrep -v "mega" )\
> window_refinement.${freq_code}_input.tsv
for CNV in DEL DUP; do
  /opt/rCNV2/analysis/sliding_windows/refine_significant_regions.py \
    --cnv-type $CNV \
    --p-is-phred \
    --cutoff ${meta_p_cutoff} \
    --prefix "${freq_code}_$CNV" \
    $CNV.sig_regions_to_refine.bed.gz \
    window_refinement.${freq_code}_input.tsv \
    $CNV.pval_matrix.bed.gz \
    refs/HPOs_by_metacohort.table.tsv
done






