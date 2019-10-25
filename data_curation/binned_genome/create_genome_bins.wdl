######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Create and annotate sequential genome-wide bins used for association testing

workflow create_genome_bins {
  Int binsize
  Int stepsize
  File contiglist
  String rCNV_bucket

  Array[Array[String]] contigs = read_tsv(contiglist)

  # Create bins
  call create_raw_bins {
    input:
      binsize=binsize,
      stepsize=stepsize,
      rCNV_bucket=rCNV_bucket
  }

  # Annotate bins per-chromosome
  scatter ( contig in contigs ) {
    call annotate_bins {
      input:
        bins=create_raw_bins.bins,
        binsize=binsize,
        stepsize=stepsize,
        contig=contig[0]
    }
  }
  call cat_annotated_bins {
    input:
      bins=annotate_bins.annotated_bins,
      binsize=binsize,
      stepsize=stepsize,
      rCNV_bucket=rCNV_bucket
  }

  # Decompose annotations
  call decompose_annotations {
    input:
      bins=cat_annotated_bins.merged_bins,
      binsize=binsize,
      stepsize=stepsize,
      rCNV_bucket=rCNV_bucket
  }

  output {
    File raw_bins = create_raw_bins.bins
    File annotated_bins = cat_annotated_bins.merged_bins
    File eigen_bins = decompose_annotations.eigen_bins
    File bin_eigenfeature_stats = decompose_annotations.stats
  }
}


# Create bins
task create_raw_bins {
  Int binsize
  Int stepsize
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
      --buffer ${binsize}000 \
      refs/GRCh37.autosomes.genome \
      ${binsize}000 \
      GRCh37.${binsize}kb_bins_${stepsize}kb_steps.raw.bed.gz

    # Move copy to master rCNV bucket
    gsutil cp GRCh37.${binsize}kb_bins_${stepsize}kb_steps.raw.bed.gz \
      ${rCNV_bucket}/cleaned_data/binned_genome/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:d566dc85f0cdca47d93a0cf1cd412d10e6d48b3b5e0ce27f96fea9a6472ce88e"
    preemptible: 1
    disks: "local-disk 50 SSD"
  }

  output {
    File bins = "GRCh37.${binsize}kb_bins_${stepsize}kb_steps.raw.bed.gz"
  }
}


# Annotate bins
task annotate_bins {
  File bins
  Int binsize
  Int stepsize
  String contig

  command <<<
    # Gather reference files
    mkdir refs/
    gsutil -m cp -r gs://rcnv_project/refs/Affy_UKBB_axiom_probes.bed.gz* refs/
    gsutil -m cp -r gs://rcnv_project/GRCh37_ref_build/* refs/

    # Annotate bins
    athena annotate-bins -z \
      --include-chroms ${contig} \
      -t refs/Affy_UKBB_axiom_probes.bed.gz -a count -n ukbbAxiom_probes \
      --ucsc-list /opt/rCNV2/data_curation/binned_genome/genome_bin_ucsc_annotations.tsv \
      --fasta refs/GRCh37.primary_assembly.fa \
      --ucsc-ref hg19 \
      ${bins} \
      GRCh37.${binsize}kb_bins_${stepsize}kb_steps.annotated.${contig}.bed.gz
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:d566dc85f0cdca47d93a0cf1cd412d10e6d48b3b5e0ce27f96fea9a6472ce88e"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 30 SSD"
    bootDiskSizeGb: "20"
  }

  output {
    File annotated_bins = "GRCh37.${binsize}kb_bins_${stepsize}kb_steps.annotated.${contig}.bed.gz"
  }
}


# Merge an array of annotated bins
task cat_annotated_bins {
  Array[File] bins
  Int binsize
  Int stepsize
  String rCNV_bucket

  command <<<
    # Get header row from first file
    zcat ${bins[0]} \
    | head -n1 \
    > header.txt

    # Merge array of bins and add header
    zcat ${sep=" " bins} \
    | fgrep -v "#" \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | cat header.txt - \
    | bgzip -c \
    > GRCh37.${binsize}kb_bins_${stepsize}kb_steps.annotated.bed.gz

    # Move copy to master rCNV bucket
    gsutil cp GRCh37.${binsize}kb_bins_${stepsize}kb_steps.annotated.bed.gz \
      ${rCNV_bucket}/cleaned_data/binned_genome/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:d566dc85f0cdca47d93a0cf1cd412d10e6d48b3b5e0ce27f96fea9a6472ce88e"
    preemptible: 1
  }

  output {
    File merged_bins = "GRCh37.${binsize}kb_bins_${stepsize}kb_steps.annotated.bed.gz"
  }
}


# Perform Eigendecomposition on a set of annotated bins
task decompose_annotations {
  File bins
  Int binsize
  Int stepsize
  String rCNV_bucket

  command <<<
    athena eigen-bins -z \
      -e 10 \
      --sqrt-transform max_recomb_rate \
      --sqrt-transform mean_recomb_rate \
      --log-transform affy6_probes \
      --log-transform affy5_probes \
      --log-transform illOmni_probes \
      --log-transform ill1M_probes \
      --log-transform ill650_probes \
      --log-transform ill550_probes \
      --log-transform ukbbAxiom_probes \
      --log-transform rmsk_LINE \
      --log-transform rmsk_SINE \
      --log-transform rmsk_LTR \
      --log-transform rmsk_DNA_repeats \
      --log-transform rmsk_simple_repeats \
      --log-transform rmsk_low_complexity_repeats \
      --log-transform segdups \
      --log-transform max_segdup_identity \
      --log-transform simple_repeats \
      --log-transform self_chain \
      --stats GRCh37.${binsize}kb_bins_${stepsize}kb_steps.eigenfeature_stats.txt \
      ${bins} \
      GRCh37.${binsize}kb_bins_${stepsize}kb_steps.annotated.eigen.bed.gz

    # Move copy to master rCNV bucket
    gsutil cp GRCh37.${binsize}kb_bins_${stepsize}kb_steps.annotated.eigen.bed.gz \
      ${rCNV_bucket}/cleaned_data/binned_genome/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:d566dc85f0cdca47d93a0cf1cd412d10e6d48b3b5e0ce27f96fea9a6472ce88e"
    preemptible: 1
    memory: "8 GB"
  }

  output {
    File eigen_bins = "GRCh37.${binsize}kb_bins_${stepsize}kb_steps.annotated.eigen.bed.gz"
    File stats = "GRCh37.${binsize}kb_bins_${stepsize}kb_steps.eigenfeature_stats.txt"
  }
}
