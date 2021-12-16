######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
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
  String athena_cloud_docker
  String rCNV_git_hash
  String rCNV_docker

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
        rCNV_bucket=rCNV_bucket,
        athena_cloud_docker=athena_cloud_docker,
        rCNV_git_hash=rCNV_git_hash
    }
    call get_expression_data {
      input:
        gtf=gtf,
        prefix=gtf_prefix,
        contig=contig[0],
        rCNV_bucket=rCNV_bucket,
        athena_cloud_docker=athena_cloud_docker,
        rCNV_git_hash=rCNV_git_hash
    }
    call get_chromatin_data {
      input:
        gtf=gtf,
        prefix=gtf_prefix,
        contig=contig[0],
        rCNV_bucket=rCNV_bucket,
        athena_cloud_docker=athena_cloud_docker,
        rCNV_git_hash=rCNV_git_hash
    }
    call get_protein_data {
      input:
        gtf=gtf,
        prefix=gtf_prefix,
        contig=contig[0],
        rCNV_bucket=rCNV_bucket,
        athena_cloud_docker=athena_cloud_docker,
        rCNV_git_hash=rCNV_git_hash
    }
    call get_constraint_data {
      input:
        gtf=gtf,
        ref_fasta=ref_fasta,
        ref_fasta_idx=ref_fasta_idx,
        prefix=gtf_prefix,
        contig=contig[0],
        rCNV_bucket=rCNV_bucket,
        athena_cloud_docker=athena_cloud_docker,
        rCNV_git_hash=rCNV_git_hash
    }
    call get_variation_data {
      input:
        gtf=gtf,
        ref_fasta=ref_fasta,
        ref_fasta_idx=ref_fasta_idx,
        prefix=gtf_prefix,
        contig=contig[0],
        rCNV_bucket=rCNV_bucket,
        athena_cloud_docker=athena_cloud_docker,
        rCNV_git_hash=rCNV_git_hash
    }
    call join_data as join_data_no_variation {
      input:
        metadata_tables=[get_genomic_data.metadata_table, get_expression_data.metadata_table, 
                         get_chromatin_data.metadata_table, get_protein_data.metadata_table, 
                         get_constraint_data.metadata_table],
        prefix="${gtf_prefix}.merged_features.no_variation.${contig[0]}",
        rCNV_docker=rCNV_docker
    }
    call join_data as join_data_with_variation {
      input:
        metadata_tables=[get_genomic_data.metadata_table, get_expression_data.metadata_table, 
                         get_chromatin_data.metadata_table, get_protein_data.metadata_table, 
                         get_constraint_data.metadata_table, get_variation_data.metadata_table],
        prefix="${gtf_prefix}.merged_features.${contig[0]}",
        rCNV_docker=rCNV_docker
    }
  }

  # Concatenate across contigs & upload to rCNV bucket
  call cat_metadata as merge_genomic_data {
    input:
      data_shards=get_genomic_data.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.genomic_features",
      eigen_prefix="genomic_eigenfeature",
      rCNV_bucket=rCNV_bucket,
      athena_cloud_docker=athena_cloud_docker
  }
  call cat_metadata as merge_expression_data {
    input:
      data_shards=get_expression_data.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.expression_features",
      eigen_prefix="expression_eigenfeature",
      rCNV_bucket=rCNV_bucket,
      athena_cloud_docker=athena_cloud_docker
  }
  call cat_metadata as merge_chromatin_data {
    input:
      data_shards=get_chromatin_data.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.chromatin_features",
      eigen_prefix="chromatin_eigenfeature",
      rCNV_bucket=rCNV_bucket,
      athena_cloud_docker=athena_cloud_docker
  }
  call cat_metadata as merge_protein_data {
    input:
      data_shards=get_protein_data.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.protein_features",
      eigen_prefix="protein_eigenfeature",
      rCNV_bucket=rCNV_bucket,
      athena_cloud_docker=athena_cloud_docker
  }
  call cat_metadata as merge_constraint_data {
    input:
      data_shards=get_constraint_data.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.constraint_features",
      eigen_prefix="constraint_eigenfeature",
      rCNV_bucket=rCNV_bucket,
      athena_cloud_docker=athena_cloud_docker
  }
  call cat_metadata as merge_variation_data {
    input:
      data_shards=get_variation_data.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.variation_features",
      eigen_prefix="variation_eigenfeature",
      rCNV_bucket=rCNV_bucket,
      athena_cloud_docker=athena_cloud_docker
  }
  call cat_metadata as merge_joined_data_no_variation {
    input:
      data_shards=join_data_no_variation.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.all_features.no_variation",
      eigen_prefix="joined_eigenfeature_no_variation",
      rCNV_bucket=rCNV_bucket,
      athena_cloud_docker=athena_cloud_docker
  }
  call cat_metadata as merge_joined_data {
    input:
      data_shards=join_data_with_variation.metadata_table,
      eigenfeatures_min_var_exp=eigenfeatures_min_var_exp,
      prefix="${gtf_prefix}.all_features",
      eigen_prefix="joined_eigenfeature",
      rCNV_bucket=rCNV_bucket,
      athena_cloud_docker=athena_cloud_docker
  }

  # Outputs
  output {
    File genes_genomic_metadata = merge_genomic_data.merged_data
    File genes_genomic_eigen_metadata = merge_genomic_data.merged_data_eigen
    File genes_expression_metadata = merge_expression_data.merged_data
    File genes_expression_eigen_metadata = merge_expression_data.merged_data_eigen
    File genes_chromatin_metadata = merge_chromatin_data.merged_data
    File genes_chromatin_eigen_metadata = merge_chromatin_data.merged_data_eigen
    File genes_protein_metadata = merge_protein_data.merged_data
    File genes_protein_eigen_metadata = merge_protein_data.merged_data_eigen
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
  String athena_cloud_docker
  String rCNV_git_hash

  command <<<
    set -e

    # Clone rCNV repo
    cwd=`pwd`
    cd /opt/ && \
    git clone https://github.com/talkowski-lab/rCNV2.git && \
    cd rCNV2 && \
    git checkout "${rCNV_git_hash}" && \
    cd $cwd

    # Download necessary data & references
    mkdir refs/
    gsutil -m cp \
      ${rCNV_bucket}/analysis/analysis_refs/** \
      ${rCNV_bucket}/refs/** \
      ${rCNV_bucket}/cleaned_data/genes/annotations/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
      refs/

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
      --gnomad-constraint refs/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
      --no-scaling \
      --outbed ${prefix}.genomic_features.${contig}.bed.gz \
      --bgzip \
      subset.gtf.gz
  >>>

  runtime {
    docker: "${athena_cloud_docker}"
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
  String athena_cloud_docker
  String rCNV_git_hash

  command <<<
    set -e

    # Clone rCNV repo
    cwd=`pwd`
    cd /opt/ && \
    git clone https://github.com/talkowski-lab/rCNV2.git && \
    cd rCNV2 && \
    git checkout "${rCNV_git_hash}" && \
    cd $cwd

    # Download precomputed GTEx expression matrices (generated by preprocess_GTEx.py)
    gsutil -m cp -r \
      ${rCNV_bucket}/cleaned_data/genes/annotations/gtex_stats \
      ./

    # Subset input files to chromosome of interest
    tabix -f ${gtf}
    tabix -h ${gtf} ${contig} | bgzip -c > subset.gtf.gz

    # Collect expression metadata
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
    docker: "${athena_cloud_docker}"
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
  String athena_cloud_docker
  String rCNV_git_hash

  command <<<
    set -e

    # Clone rCNV repo
    cwd=`pwd`
    cd /opt/ && \
    git clone https://github.com/talkowski-lab/rCNV2.git && \
    cd rCNV2 && \
    git checkout "${rCNV_git_hash}" && \
    cd $cwd

    # Download precomputed Roadmap chromatin data (generated by preprocess_REP.py)
    gsutil -m cp -r \
      ${rCNV_bucket}/cleaned_data/genes/annotations/roadmap_stats \
      ./

    # Download Episcores
    gsutil -m cp \
      ${rCNV_bucket}/cleaned_data/genes/annotations/episcore.Han_2018.tsv.gz \
      ./

    # Subset input files to chromosome of interest
    tabix -f ${gtf}
    tabix -h ${gtf} ${contig} | bgzip -c > subset.gtf.gz

    # Collect chromatin metadata
    /opt/rCNV2/data_curation/gene/get_gene_features.py \
      --get-chromatin \
      --roadmap-means roadmap_stats/gencode.v19.canonical.pext_filtered.REP_chromatin_stats.mean.tsv.gz \
      --roadmap-sds roadmap_stats/gencode.v19.canonical.pext_filtered.REP_chromatin_stats.sd.tsv.gz \
      --roadmap-pca roadmap_stats/gencode.v19.canonical.pext_filtered.REP_chromatin_stats.pca.tsv.gz \
      --episcore-tsv episcore.Han_2018.tsv.gz \
      --outbed ${prefix}.chromatin_features.${contig}.bed.gz \
      --bgzip \
      subset.gtf.gz
  >>>

  runtime {
    docker: "${athena_cloud_docker}"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: "20"
  }

  output {
   File metadata_table = "${prefix}.chromatin_features.${contig}.bed.gz"
  }
}


# Collect protein metadata for all genes from a single contig
task get_protein_data {
  File gtf
  String prefix
  String contig
  String rCNV_bucket
  String athena_cloud_docker
  String rCNV_git_hash

  command <<<
    set -e

    # Clone rCNV repo
    cwd=`pwd`
    cd /opt/ && \
    git clone https://github.com/talkowski-lab/rCNV2.git && \
    cd rCNV2 && \
    git checkout "${rCNV_git_hash}" && \
    cd $cwd

    # Download precomputed UniProt data (generated by preprocess_uniprot_data.py)
    gsutil -m cp \
      ${rCNV_bucket}/cleaned_data/genes/annotations/UniProt_features.cleaned.tsv.gz \
      ./

    # Subset input files to chromosome of interest
    tabix -f ${gtf}
    tabix -h ${gtf} ${contig} | bgzip -c > subset.gtf.gz

    # Collect protein metadata
    /opt/rCNV2/data_curation/gene/get_gene_features.py \
      --get-protein \
      --uniprot-tsv UniProt_features.cleaned.tsv.gz \
      --outbed ${prefix}.protein_features.${contig}.bed.gz \
      --bgzip \
      subset.gtf.gz
  >>>

  runtime {
    docker: "${athena_cloud_docker}"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 100 SSD"
    bootDiskSizeGb: "20"
  }

  output {
   File metadata_table = "${prefix}.protein_features.${contig}.bed.gz"
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
  String athena_cloud_docker
  String rCNV_git_hash

  command <<<
    set -e

    # Clone rCNV repo
    cwd=`pwd`
    cd /opt/ && \
    git clone https://github.com/talkowski-lab/rCNV2.git && \
    cd rCNV2 && \
    git checkout "${rCNV_git_hash}" && \
    cd $cwd

    # Download necessary data & references
    wget http://genic-intolerance.org/data/RVIS_Unpublished_ExACv2_March2017.txt
    gsutil -m cp \
      ${rCNV_bucket}/cleaned_data/genes/annotations/EDS.Wang_2018.tsv.gz \
      ${rCNV_bucket}/cleaned_data/genes/annotations/sHet.Cassa_2017.tsv.gz \
      ${rCNV_bucket}/cleaned_data/genes/annotations/CCDG_DS_scores.Abel_2020.tsv.gz \
      ${rCNV_bucket}/cleaned_data/genes/annotations/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
      ${rCNV_bucket}/cleaned_data/genes/annotations/exac-final-cnv.gene.scores071316 \
      ./
    wget https://doi.org/10.1371/journal.pgen.1001154.s002

    # Subset input files to chromosome of interest
    tabix -f ${gtf}
    tabix -h ${gtf} ${contig} | bgzip -c > subset.gtf.gz
    samtools faidx ${ref_fasta} ${contig} \
    | bgzip -c \
    > ref.fa.gz
    samtools faidx ref.fa.gz

    # Collect constraint metadata
    /opt/rCNV2/data_curation/gene/get_gene_features.py \
      --get-constraint \
      --ref-fasta ref.fa.gz \
      --gnomad-constraint gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
      --exac-cnv exac-final-cnv.gene.scores071316 \
      --rvis-tsv RVIS_Unpublished_ExACv2_March2017.txt \
      --eds-tsv EDS.Wang_2018.tsv.gz \
      --hi-tsv journal.pgen.1001154.s002 \
      --shet-tsv sHet.Cassa_2017.tsv.gz \
      --ccdg-tsv CCDG_DS_scores.Abel_2020.tsv.gz \
      --outbed ${prefix}.constraint_features.${contig}.bed.gz \
      --bgzip \
      subset.gtf.gz
  >>>

  runtime {
    docker: "${athena_cloud_docker}"
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
  String athena_cloud_docker
  String rCNV_git_hash

  command <<<
    set -e

    # Clone rCNV repo
    cwd=`pwd`
    cd /opt/ && \
    git clone https://github.com/talkowski-lab/rCNV2.git && \
    cd rCNV2 && \
    git checkout "${rCNV_git_hash}" && \
    cd $cwd

    # Download necessary data & references
    mkdir refs/
    gsutil -m cp \
      ${rCNV_bucket}/analysis/paper/data/misc/** \
      ${rCNV_bucket}/cleaned_data/genes/annotations/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
      refs/

    # Subset input files to chromosome of interest
    tabix -f ${gtf}
    tabix -h ${gtf} ${contig} | bgzip -c > subset.gtf.gz
    samtools faidx ${ref_fasta} ${contig} \
    | bgzip -c \
    > ref.fa.gz
    samtools faidx ref.fa.gz

    # Collect variation metadata
    /opt/rCNV2/data_curation/gene/get_gene_features.py \
      --get-variation \
      --gnomad-constraint refs/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
      --ddd-dnms refs/ddd_dnm_counts.tsv.gz \
      --asc-dnms refs/fu_asc_spark_dnm_counts.tsv.gz \
      --asc-unaffected-dnms refs/fu_asc_spark_dnm_counts.unaffecteds.tsv.gz \
      --gnomad-svs refs/gnomad_sv_nonneuro_counts.tsv.gz \
      --redin-bcas refs/redin_bca_counts.tsv.gz \
      --outbed ${prefix}.variation_features.${contig}.bed.gz \
      --bgzip \
      subset.gtf.gz
  >>>

  runtime {
    docker: "${athena_cloud_docker}"
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
  String rCNV_docker

  command <<<
    set -e

    /opt/rCNV2/data_curation/gene/join_gene_metadata.R \
      ${sep=" " metadata_tables} \
    | bgzip -c \
    > ${prefix}.bed.gz
  >>>

  runtime {
    docker: "${rCNV_docker}"
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
  String athena_cloud_docker

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
      -o ${prefix}.eigenfeatures.bed.gz \
      ${prefix}.bed.gz
      

    gsutil -m cp \
      ${prefix}*bed.gz \
      ${rCNV_bucket}/cleaned_data/genes/metadata/
  >>>

  runtime {
    docker: "${athena_cloud_docker}"
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

