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
  Float p_cutoff
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
  Float p_cutoff
  String rCNV_bucket
  String prefix

  command <<<
    # Copy CNV data
    mkdir cleaned_cnv/
    gsutil cp ${rCNV_bucket}/cleaned_data/cnv/* cleaned_cnv/

    # Iterate over metacohorts
    while read meta cohorts; do

      # Set metacohort parameters
      cnv_bed="cleaned_cnv/$meta.$freq_code.bed.gz"
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
      title="$descrip (${hpo})\n$ncase cases vs $nctrl controls in '${meta}' cohort"

      # Iterate over CNV types
      for CNV in CNV DEL DUP; do
        # Set CNV-specific parameters
        case "$CNV" in
          DEL)
            highlight_bed=/opt/rCNV2/refs/UKBB_GD.Owen_2018.DEL.bed.gz
            ;;
          DUP)
            highlight_bed=/opt/rCNV2/refs/UKBB_GD.Owen_2018.DUP.bed.gz
            ;;
          *)
            highlight_bed=/opt/rCNV2/refs/UKBB_GD.Owen_2018.bed.gz
            ;;
        esac

        # Count CNVs
        /opt/rCNV2/analysis/sliding_windows/count_cnvs_per_window.py \
          --fraction ${bin_overlap} \
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
          --cutoff ${p_cutoff} \
          --highlight-bed "$highlight_bed" \
          --highlight-name "Known GDs (Owen 2018)" \
          --title "$title" \
          "$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz" \
          "$meta.${prefix}.${freq_code}.$CNV.sliding_window"
      done

      # Generate Miami & QQ plots
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

    # Copy results to output bucket
    gsutil cp *.sliding_window.stats.bed.gz* \
      "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/stats/"
    gsutil cp *.sliding_window.*.png \
      "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/plots/"
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:c2b45b470f6b897269be2aa5bb82bf90d3f3aa626bfa702c038b19e2568afe26"
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

