######################
#    rCNV Project    #
######################

# Copyright (c) 2019-2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Permutation subroutine for case-control CNV burden testing of genes


workflow scattered_gene_burden_perm_test {
  String hpo
  File metacohort_list
  File metacohort_sample_table
  String freq_code
  File gtf
  Int pad_controls
  String weight_mode
  Float min_cds_ovr_del
  Float min_cds_ovr_dup
  Int max_genes_per_cnv
  Float p_cutoff
  Int n_pheno_perms
  String meta_model_prefix
  String rCNV_bucket
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
        freq_code=freq_code,
        gtf=gtf,
        pad_controls=pad_controls,
        weight_mode=weight_mode,
        min_cds_ovr_del=min_cds_ovr_del,
        min_cds_ovr_dup=min_cds_ovr_dup,
        max_genes_per_cnv=max_genes_per_cnv,
        p_cutoff=p_cutoff,
        meta_model_prefix=meta_model_prefix,
        perm_idx=idx,
        rCNV_bucket=rCNV_bucket,
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
  String freq_code
  File gtf
  Int pad_controls
  String weight_mode
  Float min_cds_ovr_del
  Float min_cds_ovr_dup
  Int max_genes_per_cnv
  Float p_cutoff
  String meta_model_prefix
  Int perm_idx
  String rCNV_bucket
  String prefix
  String cache_string

  command <<<
    set -e

    i=$(( ${perm_idx} + 1 ))

    # Copy CNV data and analysis refs
    mkdir cleaned_cnv/
    gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/* cleaned_cnv/
    gsutil -m cp -r gs://rcnv_project/cleaned_data/genes ./
    mkdir refs/
    gsutil -m cp ${rCNV_bucket}/refs/GRCh37.*.bed.gz refs/
    gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/

    # Shuffle phenotypes for each metacohort CNV dataset, and restrict CNVs from
    # phenotype of relevance
    if ! [ -e shuffled_cnv/ ]; then
      mkdir shuffled_cnv/
    fi
    yes $i | head -n1000000 > seed_$i.txt
    while read meta cohorts; do

      cnvbed="cleaned_cnv/$meta.${freq_code}.bed.gz"

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

      # Set CNV-specific parameters
      case "$CNV" in
        DEL)
          min_cds_ovr=${min_cds_ovr_del}
          ;;
        DUP)
          min_cds_ovr=${min_cds_ovr_dup}
          ;;
      esac

      # Perform association test for each metacohort
      while read meta cohorts; do

        echo -e "[$( date )] Starting permutation $i for $CNV in ${hpo} from $meta...\n"

        cnv_bed="shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed.gz"

        # Count CNVs
        /opt/rCNV2/analysis/genes/count_cnvs_per_gene.py \
          --pad-controls ${pad_controls} \
          --weight-mode ${weight_mode} \
          --min-cds-ovr $min_cds_ovr \
          --max-genes ${max_genes_per_cnv} \
          -t $CNV \
          --hpo ${hpo} \
          --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
          --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
          --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
          -z \
          --verbose \
          -o "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.bed.gz" \
          $cnv_bed \
          ${gtf}
        tabix -f "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.bed.gz"

        # Perform burden test
        /opt/rCNV2/analysis/genes/gene_burden_test.R \
          --pheno-table ${metacohort_sample_table} \
          --cohort-name $meta \
          --cnv $CNV \
          --case-hpo ${hpo} \
          --bgzip \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.bed.gz" \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz"
        tabix -f "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz"

      done < <( fgrep -v "mega" ${metacohort_list} )

      # Perform meta-analysis
      while read meta cohorts; do
        echo -e "$meta\t$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz"
      done < <( fgrep -v mega ${metacohort_list} ) \
      > ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.input.txt
      /opt/rCNV2/analysis/genes/gene_meta_analysis.R \
        --model ${meta_model_prefix} \
        --p-is-phred \
        --spa \
        ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.input.txt \
        ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.perm_$i.bed
      bgzip -f ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.perm_$i.bed
      tabix -f ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.perm_$i.bed.gz

      # Copy results to output bucket
      gsutil cp \
        ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.perm_$i.bed.gz \
        "${rCNV_bucket}/analysis/gene_burden/${prefix}/${freq_code}/permutations/"

    done
    echo "${cache_string}" > complete.txt
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:fde16e782393b0e6319dc904e59137dcede0ec0d16e634ccb91a4b8de7a1565f"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  # Note: must delocalize _something_ otherwise Cromwell will bypass this step
  output {
    File completion_marker = "complete.txt"
  }
}