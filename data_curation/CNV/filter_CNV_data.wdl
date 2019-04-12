######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Filter raw CNV data to rare and ultra-rare subsets


import "https://api.firecloud.org/ga4gh/v1/tools/rCNV:filter_cnvs_singleCohort/versions/23/plain-WDL/descriptor" as filter_single


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

  output {
    Array[File] rCNVs = filter_cohort.rCNVs
    Array[File] rCNVs_idx = filter_cohort.rCNVs_idx
    Array[File] uCNVs = filter_cohort.uCNVs
    Array[File] uCNVs_idx = filter_cohort.uCNVs_idx
  }
}

