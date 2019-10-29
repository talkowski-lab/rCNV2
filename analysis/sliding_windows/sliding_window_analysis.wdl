######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control CNV burdens in sliding windows, genome-wide


workflow sliding_window_analysis {
  File phenotype_list
  File metacohort_list
  File metacohort_sample_table
  File binned_genome
  Float bin_overlap
  Int pad_controls
  Float p_cutoff
  Float meta_p_cutoff
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
    
    # # Run uCNV assocation tests per phenotype
    # call burden_test as uCNV_burden_test {
    #   input:
    #     hpo=pheno[1],
    #     metacohort_list=metacohort_list,
    #     metacohort_sample_table=metacohort_sample_table,
    #     freq_code="uCNV",
    #     binned_genome=binned_genome,
    #     bin_overlap=bin_overlap,
    #     rCNV_bucket=rCNV_bucket,
    #     prefix=pheno[0]
    # }

    # Perform meta-analysis of rCNV association statistics
    call meta_analysis as rCNV_meta_analysis {
      input:
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        meta_p_cutoff=meta_p_cutoff,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0]
    }
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
    docker: "talkowski/rcnv@sha256:d1bb9b68186d6bf1c1995885e5d677a276c8b2c9dbed72e10a2e33f24cb1edd5"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    # Array[File] stats_beds = glob("*.sliding_window.stats.bed.gz")
    # Array[File] stats_bed_idxs = glob("*.sliding_window.stats.bed.gz.tbi")
    # Array[File] plots = glob("*.sliding_window.*.png")
  }
}


# Run meta-analysis across metacohorts for a single phenotype
task meta_analysis {
  String hpo
  File metacohort_list
  File metacohort_sample_table
  String freq_code
  Float meta_p_cutoff
  String rCNV_bucket
  String prefix

  command <<<
    set -e

    # Copy burden stats
    gsutil -m cp \
      ${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/stats/** \
      ./

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
      "${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis"

    # Copy results to output bucket
    gsutil -m cp *.sliding_window.meta_analysis.stats.bed.gz* \
      "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/stats/"
    gsutil -m cp *.sliding_window.meta_analysis.*.jpg \
      "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/plots/"
    gsutil -m cp *.sliding_window.meta_analysis.*.png \
      "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/plots/"
  >>>

  output {}

  runtime {
    docker: "talkowski/rcnv@sha256:d1bb9b68186d6bf1c1995885e5d677a276c8b2c9dbed72e10a2e33f24cb1edd5"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }  
}
