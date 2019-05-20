######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Filter raw CNV data to rare and ultra-rare subsets


import "https://api.firecloud.org/ga4gh/v1/tools/rCNV:filter_cnvs_singleCohort/versions/36/plain-WDL/descriptor" as filter_single


workflow filter_CNV_data {
  Array[String] cohorts
  Array[Int] sample_sizes
  Array[File] raw_CNVs
  File contiglist
  String rCNV_bucket

  Array[Int] cohort_idxs = range(length(cohorts))

  # Scatter over cohorts
  scatter ( i in cohort_idxs ) {
    # Call subworkflow to perform CNV filtering for a single cohort
    call filter_single.filter_cnvs_singleCohort as filter_cohort {
      input:
        cohort=cohorts[i],
        sample_size=sample_sizes[i],
        raw_CNVs=raw_CNVs[i],
        contiglist=contiglist,
        rCNV_bucket=rCNV_bucket
    }
  }

  # Make master CASE and CTRL subsets
  call combine_subsets as combine_case_rCNVs {
    input:
      beds=filter_cohort.case_rCNVs,
      bed_idxs=filter_cohort.case_rCNVs_idx,
      output_bucket="${rCNV_bucket}/cleaned_data/cnv",
      prefix="ALL.rCNV.CASE"
  }
  call combine_subsets as combine_control_rCNVs {
    input:
      beds=filter_cohort.control_rCNVs,
      bed_idxs=filter_cohort.control_rCNVs_idx,
      output_bucket="${rCNV_bucket}/cleaned_data/cnv",
      prefix="ALL.rCNV.CTRL"
  }
  call combine_subsets as combine_case_urCNVs {
    input:
      beds=filter_cohort.case_urCNVs,
      bed_idxs=filter_cohort.case_urCNVs_idx,
      output_bucket="${rCNV_bucket}/cleaned_data/cnv",
      prefix="ALL.urCNV.CASE"
  }
  call combine_subsets as combine_control_urCNVs {
    input:
      beds=filter_cohort.control_urCNVs,
      bed_idxs=filter_cohort.control_urCNVs_idx,
      output_bucket="${rCNV_bucket}/cleaned_data/cnv",
      prefix="ALL.urCNV.CTRL"
  }

  output {
    Array[File] rCNVs = filter_cohort.rCNVs
    Array[File] rCNVs_idx = filter_cohort.rCNVs_idx
    File case_rCNVs = combine_case_rCNVs.combined_bed
    File case_rCNVs_idx = combine_case_rCNVs.combined_bed_idx
    File control_rCNVs = combine_control_rCNVs.combined_bed
    File control_rCNVs_idx = combine_control_rCNVs.combined_bed_idx
    Array[File] uCNVs = filter_cohort.uCNVs
    Array[File] uCNVs_idx = filter_cohort.uCNVs_idx
    File case_urCNVs = combine_case_urCNVs.combined_bed
    File case_urCNVs_idx = combine_case_urCNVs.combined_bed_idx
    File control_urCNVs = combine_control_urCNVs.combined_bed
    File control_urCNVs_idx = combine_control_urCNVs.combined_bed_idx
  }
}


# Merge all case or control subsets across cohorts
task combine_subsets {
  Array[File] beds
  Array[File] bed_idxs
  String output_bucket
  String prefix

  command <<<
    tabix -H ${beds[0]} > header.txt
    zcat ${sep=" " beds} \
    | fgrep -v "#" \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V -k6,6V \
    | cat header.txt - \
    | bgzip -c \
    > "${prefix}.bed.gz"
    tabix -f "${prefix}.bed.gz"
    # Copy to google bucket
    gsutil cp "${prefix}.bed.gz" ${output_bucket}/
    gsutil cp "${prefix}.bed.gz.tbi" ${output_bucket}/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:745e7ad66ae7218035cec2c23488a521b11f313989f8312e85ca1429b38c5903"
    preemptible: 1
  }

  output {
    File combined_bed = "${prefix}.bed.gz"
    File combined_bed_idx = "${prefix}.bed.gz.tbi"
  }
}

