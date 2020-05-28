######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Parallelize segment permutations (matched by gene) for rCNV2 paper analyses


workflow segment_permutation {
  Int total_n_perms
  Int perms_per_shard
  String perm_prefix
  String rCNV_bucket

  call enumerate_seeds {
    input:
      total_n_perms=total_n_perms,
      perms_per_shard=perms_per_shard
  }

  scatter ( seed in enumerate_seeds.seeds ) {
    # Permute segments and literature GDs
    call perm_shard_bygene {
      input:
        seed=seed,
        n_perms=perms_per_shard,
        perm_prefix="${perm_prefix}.bygene.starting_seed_${seed}",
        rCNV_bucket=rCNV_bucket
    }
    call perm_shard_bygene_litGDs {
      input:
        seed=seed,
        n_perms=perms_per_shard,
        perm_prefix="${perm_prefix}.lit_GDs.bygene.starting_seed_${seed}",
        rCNV_bucket=rCNV_bucket
    }
    # Annotate separately to allow for caching of permutation step
    call annotate_shard as annotate_segs {
      input:
        perm_table=perm_shard_bygene.perm_table,
        hpomatch=true,
        perm_prefix="${perm_prefix}.bygene.starting_seed_${seed}",
        rCNV_bucket=rCNV_bucket
    }
    call annotate_shard as annotate_litGDs {
      input:
        perm_table=perm_shard_bygene_litGDs.perm_table,
        hpomatch=false,
        perm_prefix="${perm_prefix}.lit_GDs.bygene.starting_seed_${seed}",
        rCNV_bucket=rCNV_bucket
    }
  }

  call merge_perms as merge_segs {
    input:
      annotated_tables=annotate_segs.annotated_table,
      perm_prefix="${perm_prefix}.${total_n_perms}_permuted_segments_bygene",
      rCNV_bucket=rCNV_bucket
  }
  call merge_perms as merge_litGDs {
    input:
      annotated_tables=annotate_litGDs.annotated_table,
      perm_prefix="${perm_prefix}.lit_GDs.${total_n_perms}_permuted_segments_bygene",
      rCNV_bucket=rCNV_bucket
  }
}


# Enumerate all permutation seeds
task enumerate_seeds {
  Int total_n_perms
  Int perms_per_shard

  command <<<
    set -e

    # Make seeds
    seq 1 ${perms_per_shard} ${total_n_perms} > seeds.txt
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:82b781b4374b85258457084abe4ca7b9d12c4f9b64471b7f336fc16279f742fb"
    preemptible: 1
  }

  output {
    Array[String] seeds = read_lines("seeds.txt")
  }
}


# Single shard of segment permutation (matching on number of genes)
task perm_shard_bygene {
  String seed
  Int n_perms
  String perm_prefix
  String rCNV_bucket

  command <<<
    set -e

    # Localize necessary reference files
    mkdir refs
    gsutil -m cp \
      ${rCNV_bucket}/refs/GRCh37.autosomes.genome \
      ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
      refs/
    gsutil -m cp \
      ${rCNV_bucket}/results/segment_association/* \
      ./

    # Reformat gene coordinates
    zcat refs/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
    | cut -f1-4 \
    | bgzip -c \
    > refs/gencode.v19.canonical.pext_filtered.bed.gz

    # Permute segments, matching on number of genes
    /opt/rCNV2/analysis/paper/scripts/gene_association/shuffle_gene_blocks.py \
      --genome refs/GRCh37.autosomes.genome \
      --n-perms ${n_perms} \
      --first-seed ${seed} \
      --outfile ${perm_prefix}.tsv.gz \
      --gzip \
      <( zcat rCNV.final_segments.loci.bed.gz \
         | awk -v FS="\t" -v OFS="\t" '{ print $4, $5, $NF }' ) \
      refs/gencode.v19.canonical.pext_filtered.bed.gz
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:82b781b4374b85258457084abe4ca7b9d12c4f9b64471b7f336fc16279f742fb"
    preemptible: 1
  }

  output {
    File perm_table = "${perm_prefix}.tsv.gz"
  }  
}


# Single shard of segment permutation for literature GDs (matching on number of genes)
task perm_shard_bygene_litGDs {
  String seed
  Int n_perms
  String perm_prefix
  String rCNV_bucket

  command <<<
    set -e

    # Localize necessary reference files
    mkdir refs
    gsutil -m cp \
      ${rCNV_bucket}/refs/GRCh37.autosomes.genome \
      ${rCNV_bucket}/cleaned_data/genes/metadata/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
      ${rCNV_bucket}/analysis/paper/data/large_segments/lit_GDs.*.bed.gz \
      refs/
    gsutil -m cp \
      ${rCNV_bucket}/results/segment_association/* \
      ./

    # Reformat gene coordinates
    zcat refs/gencode.v19.canonical.pext_filtered.all_features.bed.gz \
    | cut -f1-4 \
    | bgzip -c \
    > refs/gencode.v19.canonical.pext_filtered.bed.gz

    # Pool & reformat GDs for permutation
    zcat refs/lit_GDs.*.bed.gz \
    | fgrep -v "#" \
    | sort -Vk1,1 -k2,2n -k3,3n -k5,5V \
    | awk -v FS="\t" -v OFS="\t" \
      '{ print $4, $5, $NF }' \
    | cat <( echo -e "#block_id\tcnv\tgenes" ) - \
    > lit_GDs.pooled.genes.tsv

    # Permute segments, matching on number of genes
    /opt/rCNV2/analysis/paper/scripts/gene_association/shuffle_gene_blocks.py \
      --genome refs/GRCh37.autosomes.genome \
      --n-perms ${n_perms} \
      --first-seed ${seed} \
      --outfile ${perm_prefix}.tsv.gz \
      --gzip \
      lit_GDs.pooled.genes.tsv \
      refs/gencode.v19.canonical.pext_filtered.bed.gz
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:82b781b4374b85258457084abe4ca7b9d12c4f9b64471b7f336fc16279f742fb"
    preemptible: 1
  }

  output {
    File perm_table = "${perm_prefix}.tsv.gz"
  }  
}

# Single shard of annotation
task annotate_shard {
  File perm_table
  Boolean hpomatch
  String perm_prefix
  String rCNV_bucket

  command <<<
    set -e

    # Localize necessary references
    gsutil -m cp -r \
      ${rCNV_bucket}/cleaned_data/genes/gene_lists \
      ./
    mkdir refs
    gsutil -m cp \
      ${rCNV_bucket}/analysis/analysis_refs/test_phenotypes.list \
      refs/
    gsutil -m cp \
      ${rCNV_bucket}/results/segment_association/* \
      ./

    # Build necessary inputs
    echo -e "gnomAD_constrained\tgene_lists/gnomad.v2.1.1.lof_constrained.genes.list" > genelists_to_annotate.tsv
    echo -e "gnomAD_tolerant\tgene_lists/gnomad.v2.1.1.mutation_tolerant.genes.list" >> genelists_to_annotate.tsv
    echo -e "CLinGen_HI\tgene_lists/ClinGen.hmc_haploinsufficient.genes.list" >> genelists_to_annotate.tsv
    echo -e "CLinGen_TS\tgene_lists/ClinGen.hmc_triplosensitive.genes.list" >> genelists_to_annotate.tsv
    echo -e "DECIPHER_LoF\tgene_lists/DDG2P.hmc_lof.genes.list" >> genelists_to_annotate.tsv
    echo -e "DECIPHER_GoF\tgene_lists/DDG2P.hmc_gof.genes.list" >> genelists_to_annotate.tsv
    echo -e "OMIM\tgene_lists/HP0000118.HPOdb.genes.list" >> genelists_to_annotate.tsv
    while read nocolon hpo; do
      echo -e "$hpo\tgene_lists/$nocolon.HPOdb.genes.list"
    done < refs/test_phenotypes.list \
    > hpo_genelists.tsv
    zcat rCNV.final_segments.loci.bed.gz \
    | grep -ve '^#' \
    | awk -v FS="\t" -v OFS="\t" '{ print $4, $15 }' \
    > segment_hpos.tsv
    echo -e "ASC\trefs/asc_dnm_counts.tsv.gz" > dnm_counts_to_annotate.tsv
    echo -e "DDD\trefs/ddd_dnm_counts.tsv.gz" >> dnm_counts_to_annotate.tsv

    # Annotate permuted table
    extra_args=""
    if(hpomatch){
      extra_args="--hpo-genelists hpo_genelists.tsv --segment-hpos segment_hpos.tsv"
    }
    /opt/rCNV2/analysis/paper/scripts/large_segments/annotate_shuffled_seg_gene_blocks.py \
      --gene-sets genelists_to_annotate.tsv \
      "$extra_args" \
      --dnm-tsvs dnm_counts_to_annotate.tsv \
      --outfile ${perm_prefix}.annotated.tsv.gz \
      --gzip \
      ${perm_table}
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:82b781b4374b85258457084abe4ca7b9d12c4f9b64471b7f336fc16279f742fb"
    preemptible: 1
  }

  output {
    File annotated_table = "${perm_prefix}.annotated.tsv.gz"
  }  
}


# Merge annotated permutation results
task merge_perms {
  Array[File] annotated_tables
  String perm_prefix
  String rCNV_bucket

  command <<<
    set -e
    
    # Write header
    zcat ${annotated_tables[0]} \
    | sed -n '1p' \
    > header.tsv

    # Merge and sort permutation results
    zcat ${sep=" " annotated_tables} \
    | grep -ve '^#' \
    | sort -Vk1,1 \
    | cat header.tsv - \
    | gzip -c \
    > ${perm_prefix}.tsv.gz

    # Copy to rCNV bucket
    gsutil -m cp \
      ${perm_prefix}.tsv.gz \
      ${rCNV_bucket}/analysis/paper/data/large_segments/permutations/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:82b781b4374b85258457084abe4ca7b9d12c4f9b64471b7f336fc16279f742fb"
    preemptible: 1
    disks: "local-disk 200 SSD"
  }

  output {
    File merged_results = "${perm_prefix}.tsv.gz"
  }
}
