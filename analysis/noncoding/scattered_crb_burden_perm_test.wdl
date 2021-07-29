######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Permutation subroutine for case-control CNV burden testing of cis-regulatory blocks (CRBs)


workflow scattered_crb_burden_perm_test {
  String hpo
  File metacohort_list
  File metacohort_sample_table
  File exclusion_bed
  String freq_code
  File crbs
  File crb_elements
  String noncoding_filter
  Int pad_controls
  Float min_element_ovr
  Float min_frac_all_elements
  Float p_cutoff
  Int n_pheno_perms
  String meta_model_prefix
  String rCNV_bucket
  String rCNV_docker
  String prefix
  String cache_string
  # Note: passing cache_string as a WDL variable is required to manually 
  # override caching since rCNV data is drawn directly from GCP bucket (and not)
  # read as WDL input

  scatter ( idx in range(n_pheno_perms) ) {
    call permuted_burden_test as rCNV_perm_test {
      input:
        hpo=hpo,
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        exclusion_bed=exclusion_bed,
        freq_code=freq_code,
        noncoding_filter=noncoding_filter,
        crbs=crbs,
        crb_elements=crb_elements,
        pad_controls=pad_controls,
        min_element_ovr=min_element_ovr,
        min_frac_all_elements=min_frac_all_elements,
        p_cutoff=p_cutoff,
        meta_model_prefix=meta_model_prefix,
        perm_idx=idx,
        rCNV_bucket=rCNV_bucket,
        rCNV_docker=rCNV_docker,
        prefix=prefix,
        cache_string=cache_string
    }
  }

  output {
    File completion_marker = rCNV_perm_test.completion_marker[0]
  }
}


# Permute phenotype labels to determine empirical FDR for a single phenotype
task permuted_burden_test {
  String hpo
  File metacohort_list
  File metacohort_sample_table
  File exclusion_bed
  String freq_code
  File crbs
  File crb_elements
  String noncoding_filter
  Int pad_controls
  Float min_element_ovr
  Float min_frac_all_elements
  Float p_cutoff
  String meta_model_prefix
  Int perm_idx
  String rCNV_bucket
  String rCNV_docker
  String prefix
  String cache_string

  command <<<
    set -e

    i=$(( ${perm_idx} + 1 ))

    # Copy CNV data and other references
    mkdir cleaned_cnv/
    gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/noncoding/* cleaned_cnv/
    mkdir refs/
    gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/

    # Shuffle phenotypes for each metacohort CNV dataset, and restrict CNVs from
    # phenotype of relevance
    if ! [ -e shuffled_cnv/ ]; then
      mkdir shuffled_cnv/
    fi
    yes $i | head -n1000000 > seed_$i.txt
    while read meta cohorts; do

      cnvbed="cleaned_cnv/$meta.${freq_code}.${noncoding_filter}_noncoding.bed.gz"

      # Determine CNV size quintiles
      for CNV in DEL DUP; do
        zcat $cnvbed \
        | awk -v CNV="$CNV" '{ if ($1 !~ "#" && $5==CNV) print $3-$2 }' \
        | /opt/rCNV2/utils/quantiles.py \
          --quantiles "0,0.2,0.4,0.6,0.8,1.0" \
          --no-header \
        > $meta.$CNV.quantiles.tsv
      done

      # Shuffle phenotypes, matching by CNV type and size quintile
      for CNV in DEL DUP; do
        for qr in $( seq 1 5 ); do
          smin=$( awk -v qr="$qr" '{ if (NR==qr) print $2 }' $meta.$CNV.quantiles.tsv )
          smax=$( awk -v qr="$qr" '{ if (NR==(qr+1)) print $2 + 1 }' $meta.$CNV.quantiles.tsv )
          zcat $cnvbed | sed '1d' \
          | awk -v smin="$smin" -v smax="$smax" -v CNV="$CNV" \
            '{ if ($3-$2>=smin && $3-$2<smax && $5==CNV) print $0 }' \
          > cnv_subset.bed
          paste <( cut -f1-5 cnv_subset.bed ) \
                <( cut -f6 cnv_subset.bed | shuf --random-source seed_$i.txt ) \
          | awk -v hpo=${hpo} '{ if ($NF ~ "HEALTHY_CONTROL" || $NF ~ hpo) print $0 }'
          rm cnv_subset.bed
        done
      done \
      | sort -Vk1,1 -k2,2n -k3,3n \
      | cat <( tabix -H $cnvbed ) - \
      | bgzip -c \
      > shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed.gz
      tabix -f shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed.gz

    done < ${metacohort_list}

    # Iterate over CNV types
    for CNV in DEL DUP; do

      # Perform association test for each metacohort
      while read meta cohorts; do

        echo -e "[$( date )] Starting permutation $i for $CNV in ${hpo} from $meta...\n"

        cnv_bed="shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed.gz"

        # Count CNVs
        /opt/rCNV2/analysis/noncoding/count_cnvs_per_crb.py \
          --cnvs $cnv_bed \
          --crbs ${crbs} \
          --elements ${crb_elements} \
          --pad-controls ${pad_controls} \
          --min-element-ovr ${min_element_ovr} \
          --min-frac-all-elements ${min_frac_all_elements} \
          -t $CNV \
          --hpo ${hpo} \
          -z \
          -o "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.counts.bed.gz"
        tabix -f "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.counts.bed.gz"

        # Perform burden test
        /opt/rCNV2/analysis/generic_scripts/fisher_test_single_cohort.R \
          --pheno-table ${metacohort_sample_table} \
          --cohort-name $meta \
          --case-hpo ${hpo} \
          --keep-n-columns 4 \
          --bgzip \
          "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.counts.bed.gz" \
          "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.stats.bed.gz"
        tabix -f "$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.stats.bed.gz"

      done < <( fgrep -v "mega" ${metacohort_list} )

      # Perform meta-analysis
      while read meta cohorts; do
        echo -e "$meta\t$meta.${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.stats.bed.gz"
      done < <( fgrep -v mega ${metacohort_list} ) \
      > ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.input.txt
      /opt/rCNV2/analysis/generic_scripts/meta_analysis.R \
        --model ${meta_model_prefix} \
        --conditional-exclusion ${exclusion_bed} \
        --p-is-phred \
        --keep-n-columns 4 \
        --spa \
        ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.input.txt \
        ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.perm_$i.bed
      bgzip -f ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.perm_$i.bed
      tabix -f ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.perm_$i.bed.gz

      # Copy results to output bucket
      gsutil cp \
        ${prefix}.${freq_code}.${noncoding_filter}_noncoding.$CNV.crb_burden.meta_analysis.stats.perm_$i.bed.gz \
        "${rCNV_bucket}/analysis/crb_burden/${prefix}/${freq_code}/permutations/"

    done

    echo "${cache_string}" > complete.txt
  >>>

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  # Note: must delocalize _something_ otherwise Cromwell will bypass this step
  output {
    File completion_marker = "complete.txt"
  }
}
