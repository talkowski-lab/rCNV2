######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Collect gene-level metadata

workflow get_gene_metadata {
  File gtf
  String gtf_prefix
  File contiglist
  String rCNV_bucket

  Array[String] contigs = read_lines(contiglist)

  # Scatter over contigs and gather metadata for all genes per contig
  scatter ( contig in contigs ) {
    call get_genomic_data {
      input:
        gtf=gtf,
        prefix=gtf_prefix,
        contig=contig,
        rCNV_bucket=rCNV_bucket
    }
  }

  # Merge across contigs & upload to rCNV bucket
  call merge_metadata as merge_genomic_data {
    input:
      data_shards=get_genomic_data.metadata_table,
      prefix="${gtf_prefix}.genomic_features",
      rCNV_bucket=rCNV_bucket
  }

  output {
    genes_genomic_metadata = merge_genomic_data.merged_data
  }
}


# Collect genomic metadata for all genes from a single contig
task get_genomic_data {
  File gtf
  String prefix
  String contig
  String rCNV_bucket

  command <<<
    set -e

    # Download necessary data & references
    mkdir refs/
    gsutil -m cp ${rCNV_bucket}/analysis/analysis_refs/** refs/
    gsutil -m cp ${rCNV_bucket}/refs/** refs/
    wget ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
    gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
    bgzip Homo_sapiens.GRCh37.dna.primary_assembly.fa
    samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz

    # Subset input files to chromosome of interest
    tabix -f ${gtf}
    tabix -h ${gtf} ${contig} | bgzip -c > subset.gtf.gz
    samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz ${contig} \
    | bgzip -c \
    > ref.fa.gz
    samtools faidx ref.fa.gz

    # Prep athena input file
    for wrapper in 1; do 
      echo -e "refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz\tcoverage\trepeat_cov"
      echo -e "refs/GRCh37.somatic_hypermutable_sites.bed.gz\tcoverage\thypermutable_cov"
      echo -e "https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/wgEncodeCrgMapabilityAlign100mer.bw\tmap-mean\talign_100mer"
    done > gene_features.athena_tracklist.tsv

    # Collect genomic metadata
    /opt/rCNV2/data_curation/gene/get_gene_features.py \
      --get-genomic \
      --centro-telo-bed refs/GRCh37.centromeres_telomeres.bed.gz \
      --ref-fasta ref.fa.gz \
      --athena-tracks gene_features.athena_tracklist.tsv \
      --outbed ${prefix}.genomic_features.bed.gz \
      --bgzip \
      subset.gtf.gz
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:f15f25e75adf9a055be738659b436f911e55e362e31bd0c5fdc0b6e1e10b4c6f"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: "20"
  }

  output {
   File metadata_table = "${prefix}.genomic_features.bed.gz"
  }
}


# Generic function to merge metadata tables across chromosomes
task merge_metadata {
  Array[File] data_shards
  String prefix
  String rCNV_bucket

  command <<<
    set -e

    zcat ${data_shards[0]} \
    | head -n1 \
    > header.txt

    zcat ${sep=" " data_shards} \
    | grep -ve '^#' \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | cat header.txt - \
    | bgzip -c \
    > ${prefix}.bed.gz
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:f15f25e75adf9a055be738659b436f911e55e362e31bd0c5fdc0b6e1e10b4c6f"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: "20"
  }

  output {
   File merged_data = "${prefix}.bed.gz"
  }
}

