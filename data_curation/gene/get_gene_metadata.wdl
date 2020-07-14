######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Collect gene-level metadata

workflow get_gene_metadata {
  File gtf
  File ref_fasta
  File ref_fasta_idx
  Float eigenfeatures_min_var_exp
  String gtf_prefix
  File contiglist
  String rCNV_bucket

  Array[Array[String]] contigs = read_tsv(contiglist)

  # Scatter over contigs and gather metadata for all genes per contig
  scatter ( contig in contigs ) {
    call get_genomic_data {
      input:
        gtf=gtf,
        ref_fasta=ref_fasta,
        ref_fasta_idx=ref_fasta_idx,
        prefix=gtf_prefix,
        contig=contig[0],
        rCNV_bucket=rCNV_bucket
    }
    call get_expression_data {
      input:
        gtf=gtf,
        prefix=gtf_prefix,
        contig=contig[0],
        rCNV_bucket=rCNV_bucket
    }
    call get_chromatin_data {
      input:
        gtf=gtf,
        prefix=gtf_prefix,
        contig=contig[0],
        rCNV_bucket=rCNV_bucket
    }
    call get_constraint_data {
      input:
        gtf=gtf,
        ref_fasta=ref_fasta,
        ref_fasta_idx=ref_fasta_idx,
        prefix=gtf_prefix,
        contig=contig[0],
        rCNV_bucket=rCNV_bucket
    }
    call get_variation_data {
      input:
        gtf=gtf,
        ref_fasta=ref_fasta,
        ref_fasta_idx=ref_fasta_idx,
        prefix=gtf_prefix,
        contig=contig[0],
        rCNV_bucket=rCNV_bucket
    }
    call join_data as join_data_no_variation {
      input:
        metadata_tables=[get_genomic_data.metadata_table, get_expression_data.metadata_table, get_chromatin_data.metadata_table, get_constraint_data.metadata_table],
        prefix="${gtf_prefix}.merged_features.no_variation.${contig[0]}"
    }
    call join_data as join_data_with_variation {
      input:
        metadata_tables=[get_genomic_data.metadata_table, get_expression_data.metadata_table, get_chromatin_data.metadata_table, get_constraint_data.metadata_table, get_variation_data.metadata_table],
        prefix="${gtf_prefix}.merged_features.${contig[0]}"
    }
  }

  # Concatenate across contigs & upload to rCNV bucket
  call cat_metadata as merge_genomic_data {
    input:
      data_shards=get_genomic_data.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.genomic_features",
      eigen_prefix="genomic_eigenfeature",
      rCNV_bucket=rCNV_bucket
  }
  call cat_metadata as merge_expression_data {
    input:
      data_shards=get_expression_data.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.expression_features",
      eigen_prefix="expression_eigenfeature",
      rCNV_bucket=rCNV_bucket
  }
  call cat_metadata as merge_chromatin_data {
    input:
      data_shards=get_chromatin_data.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.chromatin_features",
      eigen_prefix="chromatin_eigenfeature",
      rCNV_bucket=rCNV_bucket
  }
  call cat_metadata as merge_constraint_data {
    input:
      data_shards=get_constraint_data.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.constraint_features",
      eigen_prefix="constraint_eigenfeature",
      rCNV_bucket=rCNV_bucket
  }
  call cat_metadata as merge_variation_data {
    input:
      data_shards=get_variation_data.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.variation_features",
      eigen_prefix="variation_eigenfeature",
      rCNV_bucket=rCNV_bucket
  }
  call cat_metadata as merge_joined_data_no_variation {
    input:
      data_shards=join_data_no_variation.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.all_features.no_variation",
      eigen_prefix="joined_eigenfeature_no_variation",
      rCNV_bucket=rCNV_bucket
  }
  call cat_metadata as merge_joined_data {
    input:
      data_shards=join_data_with_variation.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.all_features",
      eigen_prefix="joined_eigenfeature",
      rCNV_bucket=rCNV_bucket
  }

  # Outputs
  output {
    File genes_genomic_metadata = merge_genomic_data.merged_data
    File genes_genomic_eigen_metadata = merge_genomic_data.merged_data_eigen
    File genes_expression_metadata = merge_expression_data.merged_data
    File genes_expression_eigen_metadata = merge_expression_data.merged_data_eigen
    File genes_chromatin_metadata = merge_chromatin_data.merged_data
    File genes_chromatin_eigen_metadata = merge_chromatin_data.merged_data_eigen
    File genes_constraint_metadata = merge_constraint_data.merged_data
    File genes_constraint_eigen_metadata = merge_constraint_data.merged_data_eigen
    File genes_variation_metadata = merge_variation_data.merged_data
    File genes_variation_eigen_metadata = merge_variation_data.merged_data_eigen
    File genes_joined_no_variation_metadata = merge_joined_data_no_variation.merged_data
    File genes_joined_no_variation_eigen_metadata = merge_joined_data_no_variation.merged_data_eigen
    File genes_joined_metadata = merge_joined_data.merged_data
    File genes_joined_eigen_metadata = merge_joined_data.merged_data_eigen
  }
}


# Collect genomic metadata for all genes from a single contig
task get_genomic_data {
  File gtf
  File ref_fasta
  File ref_fasta_idx
  String prefix
  String contig
  String rCNV_bucket

  command <<<
    set -e

    # Download necessary data & references
    wget https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
    mkdir refs/
    gsutil -m cp ${rCNV_bucket}/analysis/analysis_refs/** refs/
    gsutil -m cp ${rCNV_bucket}/refs/** refs/

    # Subset input files to chromosome of interest
    tabix -f ${gtf}
    tabix -h ${gtf} ${contig} | bgzip -c > subset.gtf.gz
    samtools faidx ${ref_fasta} ${contig} \
    | bgzip -c \
    > ref.fa.gz
    samtools faidx ref.fa.gz

    # Prep athena input file
    for wrapper in 1; do 
      echo -e "refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz\tcoverage\trepeat_cov"
      echo -e "https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/wgEncodeCrgMapabilityAlign100mer.bw\tmap-mean\talign_100mer"
    done > gene_features.athena_tracklist.tsv

    # Collect genomic metadata
    /opt/rCNV2/data_curation/gene/get_gene_features.py \
      --get-genomic \
      --centro-telo-bed refs/GRCh37.centromeres_telomeres.bed.gz \
      --ref-fasta ref.fa.gz \
      --athena-tracks gene_features.athena_tracklist.tsv \
      --gnomad-constraint gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
      --no-scaling \
      --outbed ${prefix}.genomic_features.${contig}.bed.gz \
      --bgzip \
      subset.gtf.gz
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:5e785c158825c02088d4625c76aa34fae9a503f171487825e7dd25b1496097ca"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: "20"
  }

  output {
   File metadata_table = "${prefix}.genomic_features.${contig}.bed.gz"
  }
}


# Collect expression metadata for all genes from a single contig
task get_expression_data {
  File gtf
  String prefix
  String contig
  String rCNV_bucket

  command <<<
    set -e

    # Download precomputed GTEx expression matrices (generated by preprocess_GTEx.py)
    gsutil -m cp -r \
      ${rCNV_bucket}/cleaned_data/genes/annotations/gtex_stats \
      ./

    # Subset input files to chromosome of interest
    tabix -f ${gtf}
    tabix -h ${gtf} ${contig} | bgzip -c > subset.gtf.gz

    # Collect genomic metadata
    /opt/rCNV2/data_curation/gene/get_gene_features.py \
      --get-expression \
      --gtex-medians gtex_stats/*.GTEx_v7_expression_stats.median.tsv.gz \
      --gtex-mads gtex_stats/*.GTEx_v7_expression_stats.mad.tsv.gz \
      --gtex-pca gtex_stats/*.GTEx_v7_expression_stats.pca.tsv.gz \
      --outbed ${prefix}.expression_features.${contig}.bed.gz \
      --bgzip \
      subset.gtf.gz
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:91695cf7e7a3074040a489eecfae6ae006cfccbe19d741da9c0324ba7916758b"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: "20"
  }

  output {
   File metadata_table = "${prefix}.expression_features.${contig}.bed.gz"
  }
}


# Collect chromatin metadata for all genes from a single contig
task get_chromatin_data {
  File gtf
  String prefix
  String contig
  String rCNV_bucket

  command <<<
    set -e

    # Download precomputed Roadmap chromatin data (generated by preprocess_REP.py)
    gsutil -m cp -r \
      ${rCNV_bucket}/cleaned_data/genes/annotations/roadmap_stats \
      ./

    # Subset input files to chromosome of interest
    tabix -f ${gtf}
    tabix -h ${gtf} ${contig} | bgzip -c > subset.gtf.gz

    # Collect genomic metadata
    /opt/rCNV2/data_curation/gene/get_gene_features.py \
      --get-chromatin \
      --roadmap-means roadmap_stats/gencode.v19.canonical.pext_filtered.REP_chromatin_stats.mean.tsv.gz \
      --roadmap-sds roadmap_stats/gencode.v19.canonical.pext_filtered.REP_chromatin_stats.sd.tsv.gz \
      --roadmap-pca roadmap_stats/gencode.v19.canonical.pext_filtered.REP_chromatin_stats.pca.tsv.gz \
      --outbed ${prefix}.chromatin_features.${contig}.bed.gz \
      --bgzip \
      subset.gtf.gz
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:91695cf7e7a3074040a489eecfae6ae006cfccbe19d741da9c0324ba7916758b"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: "20"
  }

  output {
   File metadata_table = "${prefix}.chromatin_features.${contig}.bed.gz"
  }
}


# Collect constraint metadata for all genes from a single contig
task get_constraint_data {
  File gtf
  File ref_fasta
  File ref_fasta_idx
  String prefix
  String contig
  String rCNV_bucket

  command <<<
    set -e

    # Download necessary data & references
    wget https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
    wget http://genic-intolerance.org/data/RVIS_Unpublished_ExACv2_March2017.txt
    gsutil -m cp ${rCNV_bucket}/cleaned_data/genes/annotations/EDS.Wang_2018.tsv.gz ./
    wget https://storage.googleapis.com/gnomad-public/legacy/exac_browser/forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz
    wget https://doi.org/10.1371/journal.pgen.1001154.s002

    # Subset input files to chromosome of interest
    tabix -f ${gtf}
    tabix -h ${gtf} ${contig} | bgzip -c > subset.gtf.gz
    samtools faidx ${ref_fasta} ${contig} \
    | bgzip -c \
    > ref.fa.gz
    samtools faidx ref.fa.gz

    # Collect genomic metadata
    /opt/rCNV2/data_curation/gene/get_gene_features.py \
      --get-constraint \
      --ref-fasta ref.fa.gz \
      --gnomad-constraint gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
      --exac-cnv forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz \
      --rvis-tsv RVIS_Unpublished_ExACv2_March2017.txt \
      --eds-tsv EDS.Wang_2018.tsv.gz \
      --hi-tsv journal.pgen.1001154.s002 \
      --outbed ${prefix}.constraint_features.${contig}.bed.gz \
      --bgzip \
      subset.gtf.gz
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:91695cf7e7a3074040a489eecfae6ae006cfccbe19d741da9c0324ba7916758b"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: "20"
  }

  output {
   File metadata_table = "${prefix}.constraint_features.${contig}.bed.gz"
  }
}


# Collect variation metadata for all genes from a single contig
task get_variation_data {
  File gtf
  File ref_fasta
  File ref_fasta_idx
  String prefix
  String contig
  String rCNV_bucket

  command <<<
    set -e

    # Download necessary data & references
    wget https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
    mkdir refs/
    gsutil -m cp ${rCNV_bucket}/analysis/paper/data/misc/** refs/

    # Subset input files to chromosome of interest
    tabix -f ${gtf}
    tabix -h ${gtf} ${contig} | bgzip -c > subset.gtf.gz
    samtools faidx ${ref_fasta} ${contig} \
    | bgzip -c \
    > ref.fa.gz
    samtools faidx ref.fa.gz

    # Collect genomic metadata
    /opt/rCNV2/data_curation/gene/get_gene_features.py \
      --get-variation \
      --gnomad-constraint gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
      --ddd-dnms refs/ddd_dnm_counts.tsv.gz \
      --asc-dnms refs/asc_dnm_counts.tsv.gz \
      --asc-unaffected-dnms refs/asc_dnm_counts.unaffecteds.tsv.gz \
      --gnomad-svs refs/gnomad_sv_nonneuro_counts.tsv.gz \
      --redin-bcas refs/redin_bca_counts.tsv.gz \
      --outbed ${prefix}.variation_features.${contig}.bed.gz \
      --bgzip \
      subset.gtf.gz
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:91695cf7e7a3074040a489eecfae6ae006cfccbe19d741da9c0324ba7916758b"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: "20"
  }

  output {
   File metadata_table = "${prefix}.variation_features.${contig}.bed.gz"
  }
}


# Join metadata for multiple conditions
task join_data {
  Array[File] metadata_tables
  String prefix

  command <<<
    set -e

    /opt/rCNV2/data_curation/gene/join_gene_metadata.R \
      ${sep=" " metadata_tables} \
    | bgzip -c \
    > ${prefix}.bed.gz
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:91695cf7e7a3074040a489eecfae6ae006cfccbe19d741da9c0324ba7916758b"
    preemptible: 1
    disks: "local-disk 50 SSD"
    bootDiskSizeGb: "20"
  }

  output {
    File metadata_table = "${prefix}.bed.gz"
  }
}


# Generic function to merge metadata tables across chromosomes & upload to rCNV_bucket
task cat_metadata {
  Array[File] data_shards
  Float eigenfeatures_min_var_exp
  String prefix
  String eigen_prefix
  String rCNV_bucket

  command <<<
    set -e

    # Localize feature transformation metadata
    gsutil -m cp \
      ${rCNV_bucket}/analysis/analysis_refs/gene_feature_transformations.tsv \
      ./

    zcat ${data_shards[0]} \
    | head -n1 \
    > header.txt

    zcat ${sep=" " data_shards} \
    | grep -ve '^#' \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | cat header.txt - \
    | bgzip -c \
    > ${prefix}.bed.gz

    athena eigen-bins \
      --min-variance ${eigenfeatures_min_var_exp} \
      --fill-missing mean \
      --skip-columns 4 \
      --max-components 1000 \
      --transformations-tsv gene_feature_transformations.tsv \
      --prefix ${eigen_prefix} \
      --bgzip \
      ${prefix}.bed.gz \
      ${prefix}.eigenfeatures.bed.gz

    gsutil -m cp \
      ${prefix}*bed.gz \
      ${rCNV_bucket}/cleaned_data/genes/metadata/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:91ee6094cfae27626def40ff4557f30701cb4df99e2221991cc87a5361aaf796"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: "20"
  }

  output {
   File merged_data = "${prefix}.bed.gz"
   File merged_data_eigen = "${prefix}.eigenfeatures.bed.gz"
  }
}

