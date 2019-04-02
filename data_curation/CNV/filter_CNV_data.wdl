######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Filter raw CNV data to rare and ultra-rare subsets


workflow filter_CNV_data {
	String rCNV_bucket
  File contiglist

  Array[String] contigs = read_lines(contiglist)

  # Scatter over contigs
  scatter ( contig in contigs ) {
    # Generate rare and ultra-rare CNVs for each cohort per contig
    call filter_cnvs_singleChrom as filter_rare {
      input:
        contig=contig[0],
        rCNV_bucket=rCNV_bucket,
        max_freq=0.01,
        CNV_suffix=rCNV
    }
    call filter_cnvs_singleChrom as filter_ultrarare {
      input:
        contig=contig[0],
        rCNV_bucket=rCNV_bucket,
        max_freq=0.0001,
        CNV_suffix=urCNV
    }
  }

  # Merge results across chromosomes


  output {}

}


# Task to filter CNVs for all cohorts on a single chromosome
task filter_cnvs_singleChrom {
  String contig
  String rCNV_bucket
  Float max_freq
  String CNV_suffix

  command <<<
    # Copy all raw CNV data
    gsutil cp -r "$rCNV_bucket"/raw_data/cnv ./
    gsutil cp -r "$rCNV_bucket"/refs ./

    # Make master BED file of all raw CNV data
    for bed in cnv/*bed.gz; do
      tabix "$bed" ${contig} 
    done \
    | sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
    | bgzip -c \
    > all_raw_cnvs.bed.gz
    allcohorts_nsamp=$( fgrep -v "#" /opt/rCNV2/refs/rCNV_sample_counts.txt \
                        | awk '{ sum+=$2 }END{ print sum }' )

    # Filter each cohort
    while read cohort N; do
      /opt/rCNV2/data_curation/CNV/filter_cnv_bed.py \
        --chr ${contig} \
        --minsize 50000 \
        --maxsize 10000000 \
        --nsamp $N \
        --maxfreq ${max_freq} \
        --recipoverlap 0.5 \
        --dist 50000 \
        --blacklist refs/GRCh37.segDups_plus_simpleRepeats.bed.gz \
        --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
        --blacklist refs/GRCh37.Nmask.bed.gz \
        --xcov 0.3 \
        --allcohorts all_raw_cnvs.bed.gz \
        --allcohorts_nsamp $allcohorts_nsamp \
        --gnomad refs/gnomAD_v2_SV_MASTER.sites.vcf.gz \
        --gnomad-af-field POPMAX_AF \
        --bgzip \
        cnv/$cohort.raw.bed.gz \
        $cohort.${contig}.${CNV_suffix}.bed.gz
    done < <( fgrep -v "#" /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f1-2 )
  >>>

  runtime {
    docker: "talkowsi/rcnv@sha256:38d6e3bf72295d5dfc94c7f1f61f5891b18a293c0684a972c87ff433045635d0"
    preemptible: 1
    memory: "8 GB"
    disks: "local-disk 30 SSD"
  }

  output {
    Array[File] filtered_cnvs = glob("*.${CNV_suffix}.bed.gz")
  }
}