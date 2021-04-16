######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control CNV burdens per cis-regulatory block (CRB)


import "https://api.firecloud.org/ga4gh/v1/tools/rCNV:scattered_crb_burden_perm_test/versions/3/plain-WDL/descriptor" as scattered_perm


workflow crb_burden_analysis {
  File phenotype_list
  File metacohort_list
  File metacohort_sample_table
  File crbs
  File crb_elements
  String noncoding_filter
  Int n_pheno_perms
  Int pad_controls
  Int min_probes_per_crb
  Float min_frac_controls_probe_exclusion
  String meta_model_prefix
  Float min_element_ovr
  Float min_frac_all_elements
  Float p_cutoff
  Int max_manhattan_phred_p
  Float meta_secondary_p_cutoff
  Int meta_nominal_cohorts_cutoff
  String rCNV_bucket
  String fisher_cache_string
  String perm_cache_string
  String meta_cache_string

  Array[Array[String]] phenotypes = read_tsv(phenotype_list)

  Array[String] cnv_types = ["DEL", "DUP"]

  # Determine conditional cohort exclusion list based on probe density
  call build_exclusion_list {
    input:
      crbs=crbs,
      crbs_prefix=basename(crbs, '.bed.gz'),
      min_probes_per_crb=min_probes_per_crb,
      min_frac_controls_probe_exclusion=min_frac_controls_probe_exclusion,
      metacohort_list=metacohort_list,
      rCNV_bucket=rCNV_bucket,
      freq_code=freq_code
  }

  # Scatter over phenotypes
  scatter ( pheno in phenotypes ) {
    # Run assocation tests per phenotype for rCNV
    call burden_test as rCNV_burden_test {
      input:
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        noncoding_filter=noncoding_filter,
        crbs=crbs,
        crb_elements=crb_elements,
        pad_controls=pad_controls,
        min_element_ovr=min_element_ovr,
        min_frac_all_elements=min_frac_all_elements,
        p_cutoff=p_cutoff,
        max_manhattan_phred_p=max_manhattan_phred_p,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0],
        cache_string=fisher_cache_string
    }
    call coding_burden_test as rCNV_coding_burden_test {
      input:
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        crbs=crbs,
        crb_elements=crb_elements,
        pad_controls=pad_controls,
        min_element_ovr=min_element_ovr,
        min_frac_all_elements=min_frac_all_elements,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0],
        cache_string=fisher_cache_string
    }

    # Permute phenotypes to estimate empirical FDR
    call scattered_perm.scattered_crb_burden_perm_test as rCNV_perm_test {
      input:
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        exclusion_bed=build_exclusion_list.exclusion_bed,
        freq_code="rCNV",
        noncoding_filter=noncoding_filter,
        crbs=crbs,
        crb_elements=crb_elements,
        pad_controls=pad_controls,
        min_element_ovr=min_element_ovr,
        min_frac_all_elements=min_frac_all_elements,
        p_cutoff=p_cutoff,
        n_pheno_perms=n_pheno_perms,
        meta_model_prefix=meta_model_prefix,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0],
        cache_string=perm_cache_string
    }
  }


  # Determine appropriate primary P-value thresholds for meta-analysis
  scatter ( cnv in cnv_types ) {
    # Genome-wide, primary
    call calc_meta_p_cutoff as calc_genome_wide_cutoffs {
      input:
        phenotype_list=phenotype_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        noncoding_filter=noncoding_filter,
        CNV=cnv,
        n_pheno_perms=n_pheno_perms,
        fdr_target=p_cutoff,
        rCNV_bucket=rCNV_bucket,
        dummy_completion_markers=rCNV_perm_test.completion_marker,
        fdr_table_suffix="empirical_genome_wide_pval",
        p_val_column_name="meta_phred_p"
    }
  }

  # Perform meta-analysis of rCNV association statistics
  scatter ( pheno in phenotypes ) {
    call meta_analysis as rCNV_meta_analysis {
      input:
        stats_beds=rCNV_burden_test.stats_beds,
        stats_bed_idxs=rCNV_burden_test.stats_bed_idxs,
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        exclusion_bed=build_exclusion_list.exclusion_bed,
        freq_code="rCNV",
        noncoding_filter=noncoding_filter,
        meta_p_cutoff_tables=calc_genome_wide_cutoffs.bonferroni_cutoff_table,
        max_manhattan_phred_p=max_manhattan_phred_p,
        meta_model_prefix=meta_model_prefix,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0],
        cache_string=meta_cache_string
    }
    call coding_meta_analysis as rCNV_coding_meta_analysis {
      input:
        stats_beds=rCNV_coding_burden_test.stats_beds,
        stats_bed_idxs=rCNV_coding_burden_test.stats_bed_idxs,
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        exclusion_bed=build_exclusion_list.exclusion_bed,
        freq_code="rCNV",
        noncoding_filter=noncoding_filter,
        meta_model_prefix=meta_model_prefix,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0],
        cache_string=meta_cache_string
    }
  }

  output {}
}


# Build list of cohorts to exclude per CRB based on inadequate probe density
task build_exclusion_list {
  File crbs
  String crbs_prefix
  Int min_probes_per_crb
  Float min_frac_controls_probe_exclusion
  File metacohort_list
  String rCNV_bucket
  String freq_code

  command <<<
    set -euo pipefail

    # Download control probesets
    gsutil -m cp -r \
      ${rCNV_bucket}/cleaned_data/control_probesets \
      ./

    # Make inputs for conditional exclusion script
    zcat ${crbs} | cut -f1-4 | bgzip -c > crb_coords.bed.gz
    for file in control_probesets/*bed.gz; do
      echo -e "$file\t$( basename $file | sed 's/\.bed\.gz//g' )"
    done > probeset_tracks.tsv

    # Build conditional exclusion list
    /opt/rCNV2/data_curation/other/probe_based_exclusion.py \
      --outfile ${crbs_prefix}.cohort_exclusion.bed.gz \
      --probecounts-outfile ${crbs_prefix}.probe_counts.bed.gz \
      --frac-pass-outfile ${crbs_prefix}.frac_passing.bed.gz \
      --min-probes ${min_probes_per_crb} \
      --min-frac-samples ${min_frac_controls_probe_exclusion} \
      --keep-n-columns 4 \
      --bgzip \
      crb_coords.bed.gz \
      probeset_tracks.tsv \
      control_probesets/rCNV.control_counts_by_array.tsv \
      <( fgrep -v mega ${metacohort_list} )
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:db7a75beada57d8e2649ce132581f675eb47207de489c3f6ac7f3452c51ddb6e"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    File exclusion_bed = "${crbs_prefix}.cohort_exclusion.bed.gz"
    File probe_counts = "${crbs_prefix}.probe_counts.bed.gz"
    File frac_passing = "${crbs_prefix}.frac_passing.bed.gz"
  }  
}


# Run burden test for a single phenotype for all metacohorts
task burden_test {
  String hpo
  File metacohort_list
  File metacohort_sample_table
  String freq_code
  File crbs
  File crb_elements
  String noncoding_filter
  Int pad_controls
  Float min_element_ovr
  Float min_frac_all_elements
  Float p_cutoff
  Int max_manhattan_phred_p
  String rCNV_bucket
  String prefix
  String cache_string

  command <<<
    set -e

    # Copy CNV data and other references
    mkdir cleaned_cnv/
    gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/noncoding/* cleaned_cnv/
    mkdir refs/
    gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/

    # Set HPO-specific parameters
    descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
               | awk -v FS="\t" '{ print $2 }' )

    # Iterate over metacohorts
    while read meta cohorts; do
      echo $meta

      # Set metacohort-specific parameters
      cnv_bed="cleaned_cnv/$meta.${freq_code}.${noncoding_filter}_noncoding.bed.gz"
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

        # Count CNVs
        /opt/rCNV2/analysis/noncoding/count_cnvs_per_crb.py \
          --cnvs $cnv_bed \
          --crbs ${crbs} \
          --elements ${crb_elements} \
          --pad-controls ${pad_controls} \
          --min-element-ovr ${min_element_ovr} \
          --min-frac-all-elements ${min_frac_all_elements} \
          -t $CNV \
          --hpo ${hpo} \
          -z \
          -o "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.counts.bed.gz" 
        tabix -f "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.counts.bed.gz"

        # Perform burden test
        /opt/rCNV2/analysis/generic_scripts/fisher_test_single_cohort.R \
          --pheno-table ${metacohort_sample_table} \
          --cohort-name $meta \
          --case-hpo ${hpo} \
          --keep-n-columns 4 \
          --bgzip \
          "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.counts.bed.gz" \
          "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.stats.bed.gz"
        tabix -f "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.stats.bed.gz"

        # Generate Manhattan & QQ plots
        /opt/rCNV2/utils/plot_manhattan_qq.R \
          --p-col-name "fisher_phred_p" \
          --p-is-phred \
          --max-phred-p ${max_manhattan_phred_p} \
          --cutoff ${p_cutoff} \
          --label-prefix "$CNV" \
          --title "$title" \
          "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.stats.bed.gz" \
          "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden"
      done

      # Generate Miami & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --miami \
        --p-col-name "fisher_phred_p" \
        --p-is-phred \
        --max-phred-p ${max_manhattan_phred_p} \
        --cutoff ${p_cutoff} \
        --label-prefix "DUP" \
        --label-prefix-2 "DEL" \
        --title "$title" \
        "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.DUP.crb_burden.stats.bed.gz" \
        "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.DEL.crb_burden.stats.bed.gz" \
        "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.crb_burden"
    done < ${metacohort_list}


    # Copy results to output bucket
    gsutil -m cp *.crb_burden.counts.bed.gz* \
      "${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/counts/"
    gsutil -m cp *.crb_burden.stats.bed.gz* \
      "${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/stats/"
    gsutil -m cp *.crb_burden.*.png \
      "${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/plots/"

    echo "${cache_string}" > completion.txt
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:a486e234f2ddc9c58e7402a5e5c1600f4b02a9c57ab054402e1035dca1774050"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    Array[File] stats_beds = glob("*.crb_burden.stats.bed.gz")
    Array[File] stats_bed_idxs = glob("*.crb_burden.stats.bed.gz.tbi")
    File completion_token = "completion.txt"
  }
}


# Run unfiltered burden test for a single phenotype for all metacohorts
task coding_burden_test {
  String hpo
  File metacohort_list
  File metacohort_sample_table
  String freq_code
  File crbs
  File crb_elements
  Int pad_controls
  Float min_element_ovr
  Float min_frac_all_elements
  String rCNV_bucket
  String prefix
  String cache_string

  command <<<
    set -e

    # Copy CNV data and other references
    mkdir cleaned_cnv/
    gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/* cleaned_cnv/
    mkdir refs/
    gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/

    # Iterate over metacohorts
    while read meta cohorts; do
      echo $meta
      cnv_bed="cleaned_cnv/$meta.${freq_code}.bed.gz"

      # Iterate over CNV types
      for CNV in DEL DUP; do
        echo $CNV

        # Count CNVs
        /opt/rCNV2/analysis/noncoding/count_cnvs_per_crb.py \
          --cnvs $cnv_bed \
          --crbs ${crbs} \
          --elements ${crb_elements} \
          --pad-controls ${pad_controls} \
          --min-element-ovr ${min_element_ovr} \
          --min-frac-all-elements ${min_frac_all_elements} \
          -t $CNV \
          --hpo ${hpo} \
          -z \
          -o "$meta.${prefix}.${freq_code}.$CNV.crb_burden.counts.bed.gz" 
        tabix -f "$meta.${prefix}.${freq_code}.$CNV.crb_burden.counts.bed.gz"

        # Perform burden test
        /opt/rCNV2/analysis/generic_scripts/fisher_test_single_cohort.R \
          --pheno-table ${metacohort_sample_table} \
          --cohort-name $meta \
          --case-hpo ${hpo} \
          --keep-n-columns 4 \
          --bgzip \
          "$meta.${prefix}.${freq_code}.$CNV.crb_burden.counts.bed.gz" \
          "$meta.${prefix}.${freq_code}.$CNV.crb_burden.stats.bed.gz"
        tabix -f "$meta.${prefix}.${freq_code}.$CNV.crb_burden.stats.bed.gz"
      done
    done < ${metacohort_list}

    # Copy results to output bucket
    gsutil -m cp *.crb_burden.counts.bed.gz* \
      "${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/counts/"
    gsutil -m cp *.crb_burden.stats.bed.gz* \
      "${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/stats/"

    echo "${cache_string}" > completion.txt
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:a486e234f2ddc9c58e7402a5e5c1600f4b02a9c57ab054402e1035dca1774050"
    preemptible: 1
    memory: "16 GB"
    disks: "local-disk 200 HDD"
    bootDiskSizeGb: "20"
  }

  output {
    Array[File] stats_beds = glob("*.crb_burden.stats.bed.gz")
    Array[File] stats_bed_idxs = glob("*.crb_burden.stats.bed.gz.tbi")
    File completion_token = "completion.txt"
  }
}

# Aggregate all permutation results to determine empirical P-value cutoff
task calc_meta_p_cutoff {
  File phenotype_list
  File metacohort_sample_table
  String freq_code
  String CNV
  String noncoding_filter
  Int n_pheno_perms
  Float fdr_target
  String rCNV_bucket
  Array[File] dummy_completion_markers #Must delocalize something or Cromwell will bypass permutation test
  String fdr_table_suffix
  String p_val_column_name

  command <<<
    set -e 

    # Gather all permutation results and compute FDR CDFs
    mkdir perm_res/
    while read prefix hpo; do
      gsutil -m cp \
        "${rCNV_bucket}/analysis/crb_burden/$prefix/${freq_code}/permutations/$prefix.${freq_code}.${noncoding_filter}_noncoding.${CNV}.crb_burden.meta_analysis.stats.perm_*.bed.gz" \
        perm_res/
      for i in $( seq 1 ${n_pheno_perms} ); do
        p_idx=$( zcat perm_res/$prefix.${freq_code}.${noncoding_filter}_noncoding.${CNV}.crb_burden.meta_analysis.stats.perm_$i.bed.gz \
                 | sed -n '1p' | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' \
                 | fgrep -w ${p_val_column_name} | cut -f2 )
        zcat perm_res/$prefix.${freq_code}.${noncoding_filter}_noncoding.${CNV}.crb_burden.meta_analysis.stats.perm_$i.bed.gz \
        | grep -ve '^#' \
        | awk -v p_idx=$p_idx '{ print $(p_idx) }' \
        | cat <( echo "$prefix.${CNV}.$i" ) - \
        > perm_res/$prefix.${freq_code}.${noncoding_filter}_noncoding.${CNV}.crb_burden.meta_analysis.permuted_p_values.$i.txt
      done
      rm perm_res/$prefix.${freq_code}.${noncoding_filter}_noncoding.${CNV}.crb_burden.meta_analysis.stats.perm_*.bed.gz
    done < ${phenotype_list}
    paste perm_res/*.crb_burden.meta_analysis.permuted_p_values.*.txt \
    | gzip -c \
    > ${freq_code}.${noncoding_filter}_noncoding.${CNV}.permuted_pval_matrix.txt.gz

    # Analyze p-values and compute FDR
    /opt/rCNV2/analysis/sliding_windows/calc_empirical_fdr.R \
      --cnv ${CNV} \
      --fdr-target ${fdr_target} \
      --linear-fit \
      --flat-ladder \
      --plot crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.${fdr_table_suffix}_permutation_results.png \
      ${freq_code}.${noncoding_filter}_noncoding.${CNV}.permuted_pval_matrix.txt.gz \
      ${metacohort_sample_table} \
      crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.${fdr_table_suffix}

    # Also produce an optional table of flat Bonferroni P-value cutoffs
    awk -v pval=${fdr_target} -v FS="\t" -v OFS="\t" \
      '{ print $1, pval }' ${phenotype_list} \
    > crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.bonferroni_pval.hpo_cutoffs.tsv

    # Copy cutoff tables to output bucket
    gsutil -m cp crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.${fdr_table_suffix}*tsv \
      "${rCNV_bucket}/analysis/analysis_refs/"
    gsutil -m cp crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.bonferroni_pval.hpo_cutoffs.tsv \
      "${rCNV_bucket}/analysis/analysis_refs/"
    gsutil -m cp crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.${fdr_table_suffix}_permutation_results.png \
      "${rCNV_bucket}/analysis/crb_burden/empirical_fdr_plots/"
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:a486e234f2ddc9c58e7402a5e5c1600f4b02a9c57ab054402e1035dca1774050"
    preemptible: 1
    memory: "32 GB"
    disks: "local-disk 275 HDD"
    bootDiskSizeGb: "40"
  }

  output {
    File perm_results_plot = "crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.${fdr_table_suffix}_permutation_results.png"
    File p_cutoff_table = "crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.${fdr_table_suffix}.hpo_cutoffs.tsv"
    File p_cutoff_ladder = "crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.${fdr_table_suffix}.ncase_cutoff_ladder.tsv"
    File bonferroni_cutoff_table = "crb_burden.${freq_code}.${noncoding_filter}_noncoding.${CNV}.bonferroni_pval.hpo_cutoffs.tsv "
  }  
}

# Run meta-analysis across metacohorts for a single phenotype
task meta_analysis {
  Array[Array[File]] stats_beds
  Array[Array[File]] stats_bed_idxs
  String hpo
  File metacohort_list
  File metacohort_sample_table
  File exclusion_bed
  String freq_code
  String noncoding_filter
  Array[File] meta_p_cutoff_tables
  Int max_manhattan_phred_p
  String meta_model_prefix
  String rCNV_bucket
  String prefix
  String cache_string

  command <<<
    set -e

    # Copy burden counts & gene coordinates
    find / -name "*${prefix}.${freq_code}.${noncoding_filter}_noncoding.*.crb_burden.stats.bed.gz*" \
    | xargs -I {} mv {} ./
    find / -name "*crb_burden.${freq_code}.${noncoding_filter}_noncoding.*.bonferroni_pval.hpo_cutoffs.tsv*" \
    | xargs -I {} mv {} ./
    mkdir refs/
    gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/

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
                    crb_burden.${freq_code}.${noncoding_filter}_noncoding.DEL.bonferroni_pval.hpo_cutoffs.tsv )
    DUP_p_cutoff=$( awk -v hpo=${prefix} '{ if ($1==hpo) print $2 }' \
                    crb_burden.${freq_code}.${noncoding_filter}_noncoding.DUP.bonferroni_pval.hpo_cutoffs.tsv )

    # Set HPO-specific parameters
    descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
               | awk -v FS="\t" '{ print $2 }' )

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

      # Perform meta-analysis of CNV counts
      while read meta cohorts; do
        echo -e "$meta\t$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.stats.bed.gz"
      done < <( fgrep -v mega ${metacohort_list} ) \
      > ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.input.txt
      /opt/rCNV2/analysis/generic_scripts/meta_analysis.R \
        --or-corplot ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.or_corplot_grid.jpg \
        --model ${meta_model_prefix} \
        --conditional-exclusion ${exclusion_bed} \
        --p-is-phred \
        --spa \
        --keep-n-columns 4 \
        ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.input.txt \
        ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.bed
      bgzip -f ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.bed
      tabix -f ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.bed.gz

      # Generate Manhattan & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --p-col-name "meta_phred_p" \
        --p-is-phred \
        --max-phred-p ${max_manhattan_phred_p} \
        --cutoff $meta_p_cutoff \
        --label-prefix "$CNV" \
        --title "$title" \
        "${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.bed.gz" \
        "${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis"
    done

    # Generate Miami & QQ plots
    /opt/rCNV2/utils/plot_manhattan_qq.R \
      --miami \
      --p-col-name "meta_phred_p" \
      --p-is-phred \
      --max-phred-p ${max_manhattan_phred_p} \
      --cutoff $DUP_p_cutoff \
      --label-prefix "DUP" \
      --cutoff-2 $DEL_p_cutoff \
      --label-prefix-2 "DEL" \
      --title "$title" \
      "${prefix}.${freq_code}.${noncoding_filter}_noncoding.DUP.crb_burden.meta_analysis.stats.bed.gz" \
      "${prefix}.${freq_code}.${noncoding_filter}_noncoding.DEL.crb_burden.meta_analysis.stats.bed.gz" \
      "${prefix}.${freq_code}.${noncoding_filter}_noncoding.crb_burden.meta_analysis"

    # Copy results to output bucket
    gsutil -m cp *.crb_burden.*meta_analysis.stats.bed.gz* \
      "${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/stats/"
    gsutil -m cp *.crb_burden.*or_corplot_grid.jpg \
      "${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/plots/"
    gsutil -m cp *.crb_burden.*meta_analysis.*.png \
      "${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/plots/"

    # Must delocalize completion marker to prevent caching of final step
    echo "${cache_string}" > completion.txt
  >>>

  output {
    File completion_token = "completion.txt"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:a486e234f2ddc9c58e7402a5e5c1600f4b02a9c57ab054402e1035dca1774050"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }
}



# Run unfiltered meta-analysis (including coding CNVs) across metacohorts for a single phenotype
task coding_meta_analysis {
  Array[Array[File]] stats_beds
  Array[Array[File]] stats_bed_idxs
  String hpo
  File metacohort_list
  File metacohort_sample_table
  File exclusion_bed
  String freq_code
  String noncoding_filter
  String meta_model_prefix
  String rCNV_bucket
  String prefix
  String cache_string

  command <<<
    set -e

    # Copy burden counts & gene coordinates
    find / -name "*${prefix}.${freq_code}.*.crb_burden.stats.bed.gz*" \
    | xargs -I {} mv {} ./
    mkdir refs/
    gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/

    # Run meta-analysis for each CNV type
    for CNV in DEL DUP; do

      # Perform meta-analysis of CNV counts
      while read meta cohorts; do
        echo -e "$meta\t$meta.${prefix}.${freq_code}.$CNV.crb_burden.stats.bed.gz"
      done < <( fgrep -v mega ${metacohort_list} ) \
      > ${prefix}.${freq_code}.$CNV.crb_burden.meta_analysis.input.txt

      /opt/rCNV2/analysis/generic_scripts/meta_analysis.R \
        --or-corplot ${prefix}.${freq_code}.$CNV.crb_burden.or_corplot_grid.jpg \
        --model ${meta_model_prefix} \
        --conditional-exclusion ${exclusion_bed} \
        --p-is-phred \
        --keep-n-columns 4 \
        --spa \
        ${prefix}.${freq_code}.$CNV.crb_burden.meta_analysis.input.txt \
        ${prefix}.${freq_code}.$CNV.crb_burden.meta_analysis.stats.bed
      bgzip -f ${prefix}.${freq_code}.$CNV.crb_burden.meta_analysis.stats.bed
      tabix -f ${prefix}.${freq_code}.$CNV.crb_burden.meta_analysis.stats.bed.gz

    done

    # Copy results to output bucket
    gsutil -m cp *.crb_burden.*meta_analysis.stats.bed.gz* \
      "${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/stats/"

    # Must delocalize completion marker to prevent caching of final step
    echo "${cache_string}" > completion.txt
  >>>

  output {
    File completion_token = "completion.txt"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:a486e234f2ddc9c58e7402a5e5c1600f4b02a9c57ab054402e1035dca1774050"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }
}

