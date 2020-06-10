######################
#    rCNV Project    #
######################

# Copyright (c) 2019-2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Compile P-value matrices from permuted sliding window meta-analyses for rCNV2


workflow get_matrices {
  File phenotype_list
  Array[String]? cnvs = ['DEL', 'DUP']
  String rCNV_bucket
  String prefix

  Array[Array[String]] phenotypes = read_tsv(phenotype_list)

  # Scatter over phenotypes & collect matrices
  scatter ( pheno in phenotypes ) {
    call get_pheno_matrices {
      pheno=pheno[0],
      rCNV_bucket=rCNV_bucket
    }
  }

  # # Merge all phenotypes into one matrix per CNV type
  # scatter ( cnv in cnvs ) {

  # }

  output {}
}


# Get P-value matrices for DEL and DUP for a single phenotype
task get_pheno_matrices {
  String pheno
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

        p_idx=$( zcat perm_res/${pheno}.rCNV.$CNV.sliding_window.meta_analysis.stats.perm_$i.bed.gz \
                 | sed -n '1p' | sed 's/\t/\n/g' | awk -v OFS="\t" '{ if ($1=="meta_phred_p") print NR }' )

        zcat perm_res/${pheno}.rCNV.$CNV.sliding_window.meta_analysis.stats.perm_$i.bed.gz \
        | grep -ve '^#' \
        | awk -v p_idx=$p_idx '{ print $(p_idx) }' \
        | cat <( echo "${pheno}.$CNV.$i" ) - \
        > perm_res/${pheno}.rCNV.$CNV.sliding_window.meta_analysis.stats.permuted_p_values.$i.txt
      done

      paste perm_res/${pheno}.rCNV.$CNV.sliding_window.meta_analysis.stats.permuted_p_values.*.txt \
      | gzip -c \
      > ${pheno}.rCNV.$CNV.sliding_window.meta_analysis.stats.permuted_p_values.tsv.gz
    done
  >>>

  output {
    del_matrix = "${pheno}.rCNV.DEL.sliding_window.meta_analysis.stats.permuted_p_values.tsv.gz"
    dup_matrix = "${pheno}.rCNV.DUP.sliding_window.meta_analysis.stats.permuted_p_values.tsv.gz"
  }

  runtime {
    docker: "talkowski/rcnv@sha256:2176cfce20a878f1dd3288e8e534c57689c58663d0354ef68b5c7c18bf34d4b7"
    preemptible: 1
    disks: "local-disk 200 SSD"
  }  
}


# # Merge P-value matrices across all phenotypes
# task merge_matrices {
#   Array[File] matrices
#   String prefix
#   String CNV
#   String rCNV_bucket

#   command <<<
#     set -e


#   >>>
# }

