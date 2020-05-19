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
    call perm_shard_bygene {
      input:
        seed=seed,
        n_perms=perms_per_shard,
        perm_prefix="${perm_prefix}.bygene.starting_seed_${seed}",
        rCNV_bucket=rCNV_bucket
    }
    call annotate_shard {
      input:
        perm_table=perm_shard.perm_table,
        perm_prefix="${perm_prefix}.bygene.starting_seed_${seed}",
        rCNV_bucket=rCNV_bucket
    }
  }

  call merge_perms {
    input:
      perm_tables=annotate_shard.annotated_table,
      perm_prefix="${perm_prefix}.${total_n_perms}_permuted_segments_bygene",
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
    # docker: "talkowski/rcnv@sha256:af27d859af7d4675ae014cb66d1f48bc8ceb7033d4f25677cfc8ec6baac47751"
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
    # docker: "talkowski/rcnv@sha256:77adc010fe1501bf69f75988e7251f2a7b1c80fa266a9caea3a4b6f87ea3ff0f"
    preemptible: 1
  }

  output {
    File perm_table = "${perm_prefix}.tsv.gz"
  }  
}


# Single shard of annotation
task annotate_shard {
  File perm_table
  String perm_prefix
  String rCNV_bucket

  command <<<
    set -e

    # Localize necessary references
    gsutil -m cp \
      ${rCNV_bucket}/cleaned_data/gene_lists \
      ./
    mkdir refs
    gsutil -m cp \
      ${rCNV_bucket}/analysis/analysis_refs/test_phenotypes.list \
      refs/
    gsutil -m cp \
      ${rCNV_bucket}/results/segment_association/* \
      ./

    # Build necessary inputs
    cat << EOF > genelists_to_annotate.tsv
gnomAD_constrained${TAB}gene_lists/gnomad.v2.1.1.lof_constrained.genes.list
gnomAD_tolerant${TAB}gene_lists/gnomad.v2.1.1.mutation_tolerant.genes.list
CLinGen_HI${TAB}gene_lists/ClinGen.hmc_haploinsufficient.genes.list
CLinGen_TS${TAB}gene_lists/ClinGen.hmc_triplosensitive.genes.list
DECIPHER_LoF${TAB}gene_lists/DDG2P.hmc_lof.genes.list
DECIPHER_GoF${TAB}gene_lists/DDG2P.hmc_gof.genes.list
OMIM${TAB}gene_lists/HP0000118.HPOdb.genes.list
EOF
    while read nocolon hpo; do
      echo -e "${hpo}\tgene_lists/${nocolon}.HPOdb.genes.list"
    done < refs/test_phenotypes.list \
    > hpo_genelists.tsv
    zcat rCNV.final_segments.loci.bed.gz \
    | grep -ve '^#' \
    | awk -v FS="\t" -v OFS="\t" '{ print $4, $15 }' \
    > segment_hpos.tsv

    # Annotate permuted table
    /opt/rCNV2/analysis/paper/scripts/large_segments/annotate_shuffled_seg_gene_blocks.py \
      --gene-sets genelists_to_annotate.tsv \
      --hpo-genelists hpo_genelists.tsv \
      --segment-hpos segment_hpos.tsv \
      --outfile ${perm_prefix}.annotated.tsv.gz \
      --gzip \
      ${perm_prefix}.tsv.gz
  >>>

  runtime {
    # docker: "talkowski/rcnv@sha256:77adc010fe1501bf69f75988e7251f2a7b1c80fa266a9caea3a4b6f87ea3ff0f"
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
      ${perm_prefix}.bed.gz \
      ${rCNV_bucket}/analysis/paper/data/large_segments/permutations/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:77adc010fe1501bf69f75988e7251f2a7b1c80fa266a9caea3a4b6f87ea3ff0f"
    preemptible: 1
    disks: "local-disk 200 SSD"
  }

  output {
    File merged_results = "${perm_prefix}.tsv.gz"
  }
}
