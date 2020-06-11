######################
#    rCNV Project    #
######################

# Copyright (c) 2019-2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Compile permuted FDR matrices from sliding window meta-analyses for rCNV2


workflow get_matrices {
  File phenotype_list
  Float fdr_target
  String rCNV_bucket
  String prefix

  Array[Array[String]] phenotypes = read_tsv(phenotype_list)

  # Scatter over phenotypes & collect empirical fdrs
  scatter ( pheno in phenotypes ) {
    call get_fdrs_perPheno {
      input:
        pheno=pheno[0],
        fdr_target=fdr_target,
        rCNV_bucket=rCNV_bucket
    }
  }

  # Merge all phenotypes into one matrix per CNV type
  call merge_matrices as merge_DEL {
    input:
      fdr_vectors=get_fdrs_perPheno.del_fdrs,
      prefix=prefix,
      CNV="DEL",
      rCNV_bucket=rCNV_bucket
  }
  call merge_matrices as merge_DUP {
    input:
      fdr_vectors=get_fdrs_perPheno.dup_fdrs,
      prefix=prefix,
      CNV="DUP",
      rCNV_bucket=rCNV_bucket
  }

  output {}
}


# Get permuted FDRs for DEL and DUP for a single phenotype
task get_fdrs_perPheno {
  String pheno
  Float fdr_target
  String rCNV_bucket

  command <<<
    set -e

    mkdir perm_res/

    for CNV in DEL DUP; do
      echo -e "\nSTARTING $CNV\n"

      gsutil -m cp \
        "${rCNV_bucket}/analysis/sliding_windows/${pheno}/rCNV/permutations/${pheno}.rCNV.$CNV.sliding_window.meta_analysis.stats.perm_*.bed.gz" \
        perm_res/

      n_pheno_perms=$( find perm_res/ -name "${pheno}.rCNV.$CNV.sliding_window.meta_analysis.stats.perm_*.bed.gz" \
                       | awk -v FS="_" '{ print $NF }' | cut -f1 -d\. | sort -nrk1,1 | head -n1 )

      for i in $( seq 1 $n_pheno_perms ); do
        /opt/rCNV2/analysis/paper/scripts/large_segments/calc_permuted_fdr.single_pheno.R \
          --fdr-target ${fdr_target} \
          perm_res/${pheno}.rCNV.$CNV.sliding_window.meta_analysis.stats.perm_$i.bed.gz \
          perm_res/${pheno}.rCNV.$CNV.sliding_window.meta_analysis.stats.empirical_fdr.$i.txt
      done

      cat <( echo -e "${pheno}.$CNV" ) \
        perm_res/${pheno}.rCNV.$CNV.sliding_window.meta_analysis.stats.empirical_fdr.*.txt \
      > ${pheno}.rCNV.$CNV.sliding_window.meta_analysis.stats.empirical_fdrs.tsv
    done
  >>>

  output {
    File del_fdrs = "${pheno}.rCNV.DEL.sliding_window.meta_analysis.stats.empirical_fdrs.tsv"
    File dup_fdrs = "${pheno}.rCNV.DUP.sliding_window.meta_analysis.stats.empirical_fdrs.tsv"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:c2898e554c146620b4e8bd5f763fb4169370d943b969289867c79a0e061da766"
    preemptible: 1
    disks: "local-disk 200 SSD"
  }  
}


# Merge P-value matrices across all phenotypes
task merge_matrices {
  Array[File] fdr_vectors
  String prefix
  String CNV
  String rCNV_bucket

  command <<<
    set -e

    # Relocate fdr vectors
    echo -e "${sep="\n" fdr_vectors}" \
    | xargs -I {} mv {} ./

    # Paste matrices
    paste *${CNV}.sliding_window.meta_analysis.stats.empirical_fdrs.tsv \
    | gzip -c \
    > ${prefix}.rCNV.${CNV}.sliding_window.meta_analysis.stats.permuted_fdrs.tsv.gz

    # Copy to gs bucket
    gsutil -m cp \
      ${prefix}.rCNV.${CNV}.sliding_window.meta_analysis.stats.permuted_fdrs.tsv.gz \
      ${rCNV_bucket}/analysis/sliding_windows/permuted_pvalue_matrices/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:c2898e554c146620b4e8bd5f763fb4169370d943b969289867c79a0e061da766"
    preemptible: 1
    disks: "local-disk 50 SSD"
  }  
}

