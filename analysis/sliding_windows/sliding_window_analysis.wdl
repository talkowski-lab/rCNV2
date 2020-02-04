######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control CNV burdens in sliding windows, genome-wide


import "https://api.firecloud.org/ga4gh/v1/tools/rCNV:scattered_sliding_window_perm_test/versions/10/plain-WDL/descriptor" as scattered_perm


workflow sliding_window_analysis {
  File phenotype_list
  File metacohort_list
  File metacohort_sample_table
  File binned_genome
  Float bin_overlap
  Int pad_controls
  Float p_cutoff
  Int n_pheno_perms
  String meta_model_prefix
  Float meta_or_cutoff
  Int meta_nominal_cohorts_cutoff
  Float credible_interval
  Int sig_window_pad
  Int refine_max_cnv_size
  String rCNV_bucket

  Array[Array[String]] phenotypes = read_tsv(phenotype_list)

  # Scatter over phenotypes
  scatter ( pheno in phenotypes ) {
    # Run rCNV assocation tests per phenotype
    call burden_test as rCNV_burden_test {
      input:
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        binned_genome=binned_genome,
        bin_overlap=bin_overlap,
        pad_controls=pad_controls,
        p_cutoff=p_cutoff,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0]
    }

    # # Permute phenotypes to estimate empirical FDR
    # call scattered_perm.scattered_sliding_window_perm_test as rCNV_perm_test {
    #   input:
    #     hpo=pheno[1],
    #     metacohort_list=metacohort_list,
    #     metacohort_sample_table=metacohort_sample_table,
    #     freq_code="rCNV",
    #     binned_genome=binned_genome,
    #     bin_overlap=bin_overlap,
    #     pad_controls=pad_controls,
    #     p_cutoff=p_cutoff,
    #     n_pheno_perms=n_pheno_perms,
    #     meta_model_prefix=meta_model_prefix,
    #     rCNV_bucket=rCNV_bucket,
    #     prefix=pheno[0]
    # }
  }

  # DEV NOTE: REBUILDING PERMUTATION ANALYSIS
  
  # # Determine appropriate genome-wide P-value threshold for meta-analysis
  # call calc_meta_p_cutoff as rCNV_calc_meta_p_cutoff {
  #   input:
  #     phenotype_list=phenotype_list,
  #     freq_code="rCNV",
  #     n_pheno_perms=n_pheno_perms,
  #     rCNV_bucket=rCNV_bucket,
  #     dummy_completion_markers=rCNV_perm_test.completion_marker
  # }

  # Perform meta-analysis of rCNV association statistics
  scatter ( pheno in phenotypes ) {
    call meta_analysis as rCNV_meta_analysis {
      input:
        stats_beds=rCNV_burden_test.stats_beds,
        stats_bed_idxs=rCNV_burden_test.stats_bed_idxs,
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        meta_p_cutoff=p_cutoff,
        meta_model_prefix=meta_model_prefix,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0]
        # meta_p_cutoff=rCNV_calc_meta_p_cutoff.meta_p_cutoff,
    }
  }

  # Refine minimal credible regions
  scatter ( cnv in ["DEL", "DUP"] ) {
    call refine_regions as refine_rCNV_regions {
      input:
        completion_tokens=rCNV_meta_analysis.completion_token,
        phenotype_list=phenotype_list,
        metacohort_list=metacohort_list,
        binned_genome=binned_genome,
        freq_code="rCNV",
        CNV=cnv,
        meta_p_cutoff=p_cutoff,
        meta_or_cutoff=meta_or_cutoff,
        meta_nominal_cohorts_cutoff=meta_nominal_cohorts_cutoff,
        meta_model_prefix=meta_model_prefix,
        credible_interval=credible_interval,
        sig_window_pad=sig_window_pad,
        refine_max_cnv_size=refine_max_cnv_size,
        rCNV_bucket=rCNV_bucket
        # meta_p_cutoff=rCNV_calc_meta_p_cutoff.meta_p_cutoff,
    }
  }

  # Plot summary metrics for final credible regions
  call plot_region_summary as plot_rCNV_regions {
    input:
      freq_code="rCNV",
      DEL_regions=refine_rCNV_regions.final_loci[0],
      DUP_regions=refine_rCNV_regions.final_loci[1],
      rCNV_bucket=rCNV_bucket
  }

  output {
    Array[File] final_sig_regions = refine_rCNV_regions.final_loci
    Array[File] final_sig_associations = refine_rCNV_regions.final_associations
  }
}


# Run burden test for a single phenotype for all metacohorts
task burden_test {
  String hpo
  File metacohort_list
  File metacohort_sample_table
  String freq_code
  File binned_genome
  Float bin_overlap
  Int pad_controls
  Float p_cutoff
  String rCNV_bucket
  String prefix

  command <<<
    set -e

    # Copy CNV data
    mkdir cleaned_cnv/
    gsutil -m cp ${rCNV_bucket}/cleaned_data/cnv/* cleaned_cnv/

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

    # Copy results to output bucket
    gsutil -m cp *.sliding_window.stats.bed.gz* \
      "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/stats/"
    gsutil -m cp *.sliding_window.*.png \
      "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/plots/"
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:44d9ed4679ba2ed87097211fcd70f127e9f8ab34b9192b13036a7c00152b18dd"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    Array[File] stats_beds = glob("*.sliding_window.stats.bed.gz")
    Array[File] stats_bed_idxs = glob("*.sliding_window.stats.bed.gz.tbi")
    # Array[File] plots = glob("*.sliding_window.*.png")
  }
}


# Aggregate all genome-wide permutation results to determine empirical P-value cutoff
task calc_meta_p_cutoff {
  File phenotype_list
  String freq_code
  Int n_pheno_perms
  String rCNV_bucket
  Array[File] dummy_completion_markers #Must delocalize something or Cromwell will bypass permutation test

  command <<<
    # Gather all permutation results and compute FDR CDFs
    mkdir perm_res/
    while read prefix hpo; do
      echo -e "$prefix\n\n"
      gsutil -m cp \
        "${rCNV_bucket}/analysis/sliding_windows/$prefix/${freq_code}/permutations/**.stats.perm_*.bed.gz" \
        perm_res/
      for CNV in DEL DUP; do
        for i in $( seq 1 ${n_pheno_perms} ); do
          zcat perm_res/$prefix.${freq_code}.$CNV.sliding_window.meta_analysis.stats.perm_$i.bed.gz \
          | grep -ve '^#' \
          | awk '{ print $NF }' \
          | cat <( echo "$prefix.$CNV.$i" ) - \
          > perm_res/$prefix.${freq_code}.$CNV.sliding_window.meta_analysis.permuted_p_values.$i.txt
        done
      done
    done < ${phenotype_list}
    paste perm_res/*.sliding_window.meta_analysis.permuted_p_values.*.txt \
    | gzip -c \
    > ${freq_code}.permuted_pval_matrix.txt.gz

    # Calculate empirical FDR
    /opt/rCNV2/analysis/sliding_windows/calc_empirical_fdr.R \
      --plot ${freq_code}.FDR_permutation_results.png \
      ${freq_code}.permuted_pval_matrix.txt.gz \
      ${freq_code}.empirical_fdr_cutoffs.tsv

    tail -n1 ${freq_code}.empirical_fdr_cutoffs.tsv | cut -f3 \
    > meta_p_cutoff.txt
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:44d9ed4679ba2ed87097211fcd70f127e9f8ab34b9192b13036a7c00152b18dd"
    preemptible: 1
    memory: "16 GB"
    disks: "local-disk 100 HDD"
    bootDiskSizeGb: "20"
  }

  output {
    File perm_results_plot = "${freq_code}.FDR_permutation_results.png"
    Float meta_p_cutoff = read_float("meta_p_cutoff.txt")
  }  
}


# Run meta-analysis across metacohorts for a single phenotype
task meta_analysis {
  Array[Array[File]] stats_beds
  Array[Array[File]] stats_bed_idxs
  String hpo
  File metacohort_list
  File metacohort_sample_table
  String freq_code
  Float meta_p_cutoff
  String meta_model_prefix
  String rCNV_bucket
  String prefix

  command <<<
    set -e

    # Copy burden stats
    find / -name "*${prefix}.${freq_code}.*.sliding_window.stats.bed.gz*" \
    | xargs -I {} mv {} ./
    # gsutil -m cp \
    #   ${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/stats/** \
    #   ./

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

    # Run meta-analysis for each CNV type
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

    # Copy results to output bucket
    gsutil -m cp *.sliding_window.meta_analysis.stats.bed.gz* \
      "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/stats/"
    gsutil -m cp *.sliding_window.or_corplot_grid.jpg \
      "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/plots/"
    gsutil -m cp *.sliding_window.meta_analysis.*.png \
      "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/plots/"

    # Must delocalize completion marker to prevent caching of final step
    echo "DONE" > completion.txt
  >>>

  output {
    File completion_token = "completion.txt"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:44d9ed4679ba2ed87097211fcd70f127e9f8ab34b9192b13036a7c00152b18dd"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }
}


# Refine associated regions to minimal credible regions
task refine_regions {
  Array[File] completion_tokens # Must delocalize something from meta-analysis step to prevent caching
  File phenotype_list
  File metacohort_list
  File binned_genome
  String freq_code
  String CNV
  Float meta_p_cutoff
  Float meta_or_cutoff
  Int meta_nominal_cohorts_cutoff
  String meta_model_prefix
  Float credible_interval
  Int sig_window_pad
  Int refine_max_cnv_size
  String rCNV_bucket

  command <<<
    set -e

    # Download all meta-analysis stats files and necessary data
    mkdir cleaned_cnv/
    gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/* cleaned_cnv/
    mkdir stats/
    gsutil -m cp \
      ${rCNV_bucket}/analysis/sliding_windows/**.${freq_code}.**.sliding_window.meta_analysis.stats.bed.gz \
      stats/
    mkdir refs/
    gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/
    mkdir phenos/
    gsutil -m cp gs://rcnv_project/cleaned_data/phenotypes/filtered/* phenos/

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
      --max-p ${meta_p_cutoff} \
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

    # Refine associations within regions from above
    while read meta; do
      echo -e "$meta\tcleaned_cnv/$meta.${freq_code}.bed.gz\tphenos/$meta.cleaned_phenos.txt"
    done < <( cut -f1 ${metacohort_list} | fgrep -v "mega" )\
    > window_refinement.${freq_code}_metacohort_info.tsv
    /opt/rCNV2/analysis/sliding_windows/refine_significant_regions.py \
      --cnv-type ${CNV} \
      --model ${meta_model_prefix} \
      --min-p ${meta_p_cutoff} \
      --p-is-phred \
      --min-or-lower 0 \
      --retest-min-or-lower ${meta_or_cutoff} \
      --max-cnv-size ${refine_max_cnv_size} \
      --min-nominal ${meta_nominal_cohorts_cutoff} \
      --credible-interval ${credible_interval} \
      --prefix "${freq_code}_${CNV}" \
      --log ${freq_code}.${CNV}.region_refinement.log \
      ${freq_code}.${CNV}.sig_regions_to_refine.bed.gz \
      window_refinement.${freq_code}_metacohort_info.tsv \
      ${CNV}.pval_matrix.bed.gz \
      ${freq_code}.${CNV}.all_windows_labeled.bed.gz \
      ${freq_code}.${CNV}.final_regions.associations.bed \
      ${freq_code}.${CNV}.final_regions.loci.bed
    bgzip -f ${freq_code}.${CNV}.final_regions.associations.bed
    bgzip -f ${freq_code}.${CNV}.final_regions.loci.bed
    
    # Annotate final regions with genes
    gsutil -m cp -r gs://rcnv_project/cleaned_data/genes ./
    /opt/rCNV2/analysis/sliding_windows/get_genes_per_region.py \
      -o ${freq_code}.${CNV}.final_regions.loci.bed \
      ${freq_code}.${CNV}.final_regions.loci.bed.gz \
      genes/gencode.v19.canonical.gtf.gz
    bgzip -f ${freq_code}.${CNV}.final_regions.loci.bed

    # Copy results to output bucket
    gsutil -m cp ${freq_code}.${CNV}.final_regions.*.bed.gz \
      "${rCNV_bucket}/results/sliding_windows/"
    gsutil -m cp ${freq_code}.${CNV}.region_refinement.log \
      "${rCNV_bucket}/results/sliding_windows/"
  >>>

  output {
    File final_loci = "${freq_code}.${CNV}.final_regions.loci.bed.gz"
    File final_associations = "${freq_code}.${CNV}.final_regions.associations.bed.gz"
    File logfile = "${freq_code}.${CNV}.region_refinement.log"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:44d9ed4679ba2ed87097211fcd70f127e9f8ab34b9192b13036a7c00152b18dd"
    preemptible: 1
    memory: "16 GB"
    bootDiskSizeGb: "20"
    disks: "local-disk 100 HDD"
  }
}


task plot_region_summary {
  String freq_code
  File DEL_regions
  File DUP_regions
  String rCNV_bucket

  command <<<
    /opt/rCNV2/analysis/sliding_windows/regions_summary.plot.R \
      -o "${freq_code}.final_regions." \
      ${DEL_regions} \
      ${DUP_regions}

    gsutil -m cp \
      ./*.jpg \
      "${rCNV_bucket}/results/sliding_windows/plots/"
    gsutil -m cp \
      "${freq_code}.final_regions.multipanel_summary.jpg" \
      gs://rcnv_project/public/
    gsutil acl ch -u AllUsers:R gs://rcnv_project/public/*.jpg
  >>>

  output {}

  runtime {
    docker: "talkowski/rcnv@sha256:44d9ed4679ba2ed87097211fcd70f127e9f8ab34b9192b13036a7c00152b18dd"
    preemptible: 1
  }
}

