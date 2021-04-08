######################
#    rCNV Project    #
######################

# Copyright (c) 2019-2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Create and annotate sequential genome-wide bins used for association testing

workflow create_genome_bins {
  Int binsize
  Int stepsize
  Float blacklist_cov
  String rCNV_bucket

  Array[Array[String]] contigs = read_tsv(contiglist)

  # Create bins
  call create_raw_bins {
    input:
      binsize=binsize,
      stepsize=stepsize,
      bl_cov=blacklist_cov,
      rCNV_bucket=rCNV_bucket
  }

  # Annotate bins
  call annotate_bins {
    input:
      bins=create_raw_bins.bins,
      binsize=binsize,
      stepsize=stepsize
  }

  output {
    File raw_bins = create_raw_bins.bins
    File raw_bins_idx = create_raw_bins.bins_idx
    File annotated_bins = annotate_bins.annotated_bins
    File annotated_bins_idx = annotate_bins.annotated_bins_idx
  }
}


# Create bins
task create_raw_bins {
  Int binsize
  Int stepsize
  Float bl_cov
  String rCNV_bucket

  command <<<
    # Gather reference files
    mkdir refs/
    gsutil -m cp -r gs://rcnv_project/refs/GRCh37.Nmask.autosomes.bed.gz refs/
    gsutil -m cp -r gs://rcnv_project/refs/GRCh37.somatic_hypermutable_sites.bed.gz refs/
    gsutil -m cp -r gs://rcnv_project/refs/GRCh37.autosomes.genome refs/

    # Create bins
    athena make-bins -z \
      -x refs/GRCh37.Nmask.autosomes.bed.gz \
      -x refs/GRCh37.somatic_hypermutable_sites.bed.gz \
      -s ${stepsize}000 \
      --blacklist-cov ${bl_cov} \
      refs/GRCh37.autosomes.genome \
      ${binsize}000 \
      GRCh37.${binsize}kb_bins_${stepsize}kb_steps.raw.bed.gz
    tabix -f GRCh37.${binsize}kb_bins_${stepsize}kb_steps.raw.bed.gz

    # Move copy to master rCNV bucket
    gsutil cp GRCh37.${binsize}kb_bins_${stepsize}kb_steps.raw.bed.gz* \
      ${rCNV_bucket}/cleaned_data/binned_genome/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:db7a75beada57d8e2649ce132581f675eb47207de489c3f6ac7f3452c51ddb6e"
    preemptible: 1
    disks: "local-disk 50 SSD"
  }

  output {
    File bins = "GRCh37.${binsize}kb_bins_${stepsize}kb_steps.raw.bed.gz"
    File bins_idx = "GRCh37.${binsize}kb_bins_${stepsize}kb_steps.raw.bed.gz.tbi"
  }
}


# Annotate bins with probe counts
task annotate_bins {
  File bins
  Int binsize
  Int stepsize

  command <<<
    set -euo pipefail

    # Build athena input
    gsutil cp -r gs://rcnv_project/cleaned_data/control_probesets ./
    for file in control_probesets/*bed.gz; do
      echo -e "$file\tcount\t$( basename $file | sed 's/\.bed\.gz//g' )"
    done > probeset_tracks.athena.tsv

    # Annotate
    athena annotate-bins \
      --track-list probeset_tracks.athena.tsv \
      --bgzip \
      ${bins} \
      GRCh37.${binsize}kb_bins_${stepsize}kb_steps.annotated.bed.gz
    tabix -f GRCh37.${binsize}kb_bins_${stepsize}kb_steps.annotated.bed.gz

    # Move copy to master rCNV bucket
    gsutil cp GRCh37.${binsize}kb_bins_${stepsize}kb_steps.annotated.bed.gz* \
      gs://rcnv_project/cleaned_data/binned_genome/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:db7a75beada57d8e2649ce132581f675eb47207de489c3f6ac7f3452c51ddb6e"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 30 SSD"
    bootDiskSizeGb: "20"
  }

  output {
    File annotated_bins = "GRCh37.${binsize}kb_bins_${stepsize}kb_steps.annotated.bed.gz"
    File annotated_bins = "GRCh37.${binsize}kb_bins_${stepsize}kb_steps.annotated.bed.gz.tbi"
  }
}
