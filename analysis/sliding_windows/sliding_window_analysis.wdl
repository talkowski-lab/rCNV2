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
  File binned_genome
  String rCNV_bucket

  File phenotypes = read_tsv(phenotype_list)

  # Scatter over phenotypes
  scatter ( pheno in phenotypes ) {
    # Run rCNV assocation tests per phenotype
    call burden_test as rCNV_burden_test {
      input:
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        freq_code="rCNV",
        binned_genome=binned_genome,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0],
    }
    # Run uCNV assocation tests per phenotype
    call burden_test as uCNV_burden_test {}
  }

  output {}
}


# Run burden test for a single phenotype for all metacohorts
task burden_test {
  String hpo
  File metacohort_list
  String freq_code
  File binned_genome
  String rCNV_bucket
  String prefix

  command <<<
    # Copy CNV data
    mkdir cleaned_cnv/
    gsutil cp ${rCNV_bucket}/cleaned_data/cnv/* cleaned_cnv/
    # Iterate over metacohorts
    while read meta cohorts; do
      cnv_bed="cleaned_cnv/$meta.$freq_code.bed.gz"
      # Iterate over CNV types
      for CNV in CNV DEL DUP; do
        # Count CNVs
        /opt/rCNV2/analysis/sliding_windows/count_cnvs_per_window.py \
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
      done
    done < ${metacohort_list}
    # Copy results to output bucket
    gsutil cp *.sliding_window.stats.bed.gz* \
      "${rCNV_bucket}/analysis/sliding_window/${prefix}/"
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:26ca8f4b57f8c074bc08c42fab76dd7d157bebf041ff233efa7b6a5de438da9c"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    Array[File] stats_beds = glob("*.sliding_window.stats.bed.gz")
    Array[File] stats_bed_idxs = glob("*.sliding_window.stats.bed.gz.tbi")
  }
}

