######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Filter raw CNV data to rare and ultra-rare subsets

import "https://api.firecloud.org/ga4gh/v1/tools/rCNV:filter_cnvs_singleCohort/versions/47/plain-WDL/descriptor" as filter_single


workflow filter_CNV_data {
  Array[String] cohorts
  Array[Int] sample_sizes
  Array[File] raw_CNVs
  File metacohort_list
  File contiglist
  String rCNV_bucket

  Array[Array[String]] metacohorts = read_tsv(metacohort_list)
  Array[Int] cohort_idxs = range(length(cohorts))

  # Scatter over cohorts and filter CNVs
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

  # Scatter over metacohorts and combine filtered CNVs
  scatter ( metacohort in metacohorts ) {
    call combine_subsets as combine_meta_rCNVs {
      input:
        beds=filter_cohort.rCNVs,
        bed_idxs=filter_cohort.rCNVs_idx,
        cohorts=metacohort[1],
        output_bucket="${rCNV_bucket}/cleaned_data/cnv",
        prefix="${metacohort[0]}.rCNV"
    }
    call combine_subsets as combine_meta_uCNVs {
      input:
        beds=filter_cohort.uCNVs,
        bed_idxs=filter_cohort.uCNVs_idx,
        cohorts=metacohort[1],
        output_bucket="${rCNV_bucket}/cleaned_data/cnv",
        prefix="${metacohort[0]}.uCNV"
    }
  }

  output {
    Array[File] rCNVs = filter_cohort.rCNVs
    Array[File] rCNVs_idx = filter_cohort.rCNVs_idx
    Array[File] uCNVs = filter_cohort.uCNVs
    Array[File] uCNVs_idx = filter_cohort.uCNVs_idx
  }
}


# Merge all case or control subsets across cohorts
task combine_subsets {
  Array[File] beds
  Array[File] bed_idxs
  String cohorts
  String output_bucket
  String prefix

  File bedlist = write_tsv(beds)
  File idxlist = write_tsv(bed_idxs)

  command <<<
    tabix -H ${beds[0]} > header.txt
    echo "${cohorts}" | sed 's/;/\n/g' > cohorts.list
    while read cohort; do
      find / -name "$cohort.*.bed.gz" \
      | xargs -I {} zcat {} \
      | fgrep -v "#"
    done < cohorts.list \
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
    docker: "talkowski/rcnv@sha256:f2fde8ddf20b69e25b48f19194b6d4a716132c5214f3b2a65854ebe0efc64da5"
    preemptible: 1
  }

  output {
    File combined_bed = "${prefix}.bed.gz"
    File combined_bed_idx = "${prefix}.bed.gz.tbi"
  }
}

