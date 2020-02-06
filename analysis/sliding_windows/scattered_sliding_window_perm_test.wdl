######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control CNV burdens in sliding windows, genome-wide


workflow scattered_sliding_window_perm_test {
  String hpo
  File metacohort_list
  File metacohort_sample_table
  String freq_code
  File binned_genome
  Float bin_overlap
  Int pad_controls
  Float p_cutoff
  Int n_pheno_perms
  String meta_model_prefix
  String rCNV_bucket
  String prefix

	scatter ( idx in range(n_pheno_perms) ) {
	  call permuted_burden_test as rCNV_perm_test {
	    input:
	      hpo=hpo,
	      metacohort_list=metacohort_list,
	      metacohort_sample_table=metacohort_sample_table,
	      freq_code=freq_code,
	      binned_genome=binned_genome,
	      bin_overlap=bin_overlap,
	      pad_controls=pad_controls,
	      p_cutoff=p_cutoff,
        meta_model_prefix=meta_model_prefix,
	      perm_idx=idx,
	      rCNV_bucket=rCNV_bucket,
	      prefix=prefix
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
  File binned_genome
  Float bin_overlap
  Int pad_controls
  Float p_cutoff
  String meta_model_prefix
  Int perm_idx
  String rCNV_bucket
  String prefix

  command <<<
    set -e

    i=$(( ${perm_idx} + 1 ))

    # Copy CNV data
    mkdir cleaned_cnv/
    gsutil -m cp ${rCNV_bucket}/cleaned_data/cnv/* cleaned_cnv/

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
      # Perform association test for each metacohort
      while read meta cohorts; do

        echo -e "[$( date )] Starting permutation $i for $CNV in ${hpo} from $meta...\n"

        cnv_bed="shuffled_cnv/$meta.${freq_code}.pheno_shuf.bed.gz"

        # Count CNVs
        /opt/rCNV2/analysis/sliding_windows/count_cnvs_per_window.py \
          --fraction ${bin_overlap} \
          --pad-controls ${pad_controls} \
          -t $CNV \
          --hpo ${hpo} \
          -z \
          -o "$meta.${prefix}.${freq_code}.$CNV.sliding_window.counts.bed.gz" \
          $cnv_bed \
          ${binned_genome}

        # Perform burden test
        /opt/rCNV2/analysis/sliding_windows/window_burden_test.R \
          --pheno-table ${metacohort_sample_table} \
          --cohort-name $meta \
          --case-hpo ${hpo} \
          --bgzip \
          "$meta.${prefix}.${freq_code}.$CNV.sliding_window.counts.bed.gz" \
          "$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz"
          tabix -f "$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz"

      done < <( fgrep -v "mega" ${metacohort_list} )

      # Perform meta-analysis
      while read meta cohorts; do
        echo -e "$meta\t$meta.${prefix}.${freq_code}.$CNV.sliding_window.stats.bed.gz"
      done < <( fgrep -v mega ${metacohort_list} ) \
      > ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.input.txt
      /opt/rCNV2/analysis/sliding_windows/window_meta_analysis.R \
        --model ${meta_model_prefix} \
        --p-is-phred \
        ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.input.txt \
        ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.perm_$i.bed
      bgzip -f ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.perm_$i.bed
      tabix -f ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.perm_$i.bed.gz

      # Copy results to output bucket
      gsutil cp \
        ${prefix}.${freq_code}.$CNV.sliding_window.meta_analysis.stats.perm_$i.bed.gz \
        "${rCNV_bucket}/analysis/sliding_windows/${prefix}/${freq_code}/permutations/"

    done
    echo "Done" > complete.txt
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:cf77a0b11c44d82d132c30778952de73da10d20216558718237737e9c8956331"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  # Note: must delocalize _something_ otherwise Cromwell will bypass this step
  output {
    File completion_marker = "complete.txt"
  }
}