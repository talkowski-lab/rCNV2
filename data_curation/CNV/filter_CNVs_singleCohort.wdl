######################
#    rCNV Project    #
######################

# Copyright (c) 2019-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Filter raw CNV data to rare, very rare, and ultra rare subsets for a single cohort


workflow filter_cnvs_singleCohort {
  String cohort
  File raw_CNVs
  File contiglist
  String rCNV_bucket
  String rCNV_docker

  Array[Array[String]] contigs = read_tsv(contiglist)

  # Scatter over contigs
  scatter ( contig in contigs ) {
    # Generate rare, very rare, and ultra-rare CNVs for each cohort per contig
    call filter_cnvs_singleChrom as filter_rare {
      input:
        cohort=cohort,
        raw_CNVs=raw_CNVs,
        contig=contig[0],
        rCNV_bucket=rCNV_bucket,
        max_freq="0.01",
        CNV_suffix="rCNV",
        rCNV_docker=rCNV_docker
    }
    call filter_cnvs_singleChrom as filter_veryrare {
      input:
        cohort=cohort,
        raw_CNVs=raw_CNVs,
        contig=contig[0],
        rCNV_bucket=rCNV_bucket,
        max_freq="0.001",
        CNV_suffix="vCNV",
        rCNV_docker=rCNV_docker
    }
    call filter_cnvs_singleChrom as filter_ultrarare {
      input:
        cohort=cohort,
        raw_CNVs=raw_CNVs,
        contig=contig[0],
        rCNV_bucket=rCNV_bucket,
        max_freq="0.0001",
        CNV_suffix="uCNV",
        rCNV_docker=rCNV_docker
    }
  }

  # Merge results across chromosomes
  call merge_beds as merge_rare {
    input:
      beds=filter_rare.filtered_cnvs,
      prefix="${cohort}.rCNV",
        output_bucket="${rCNV_bucket}/cleaned_data/cnv",
        rCNV_docker=rCNV_docker
  }
  call merge_beds as merge_veryrare {
    input:
      beds=filter_veryrare.filtered_cnvs,
      prefix="${cohort}.vCNV",
        output_bucket="${rCNV_bucket}/cleaned_data/cnv",
        rCNV_docker=rCNV_docker
  }
  call merge_beds as merge_ultrarare {
    input:
      beds=filter_ultrarare.filtered_cnvs,
      prefix="${cohort}.uCNV",
        output_bucket="${rCNV_bucket}/cleaned_data/cnv",
        rCNV_docker=rCNV_docker
  }

  output {
    File rCNVs = merge_rare.merged_bed
    File rCNVs_idx = merge_rare.merged_bed_idx
    File vCNVs = merge_veryrare.merged_bed
    File vCNVs_idx = merge_veryrare.merged_bed_idx
    File uCNVs = merge_ultrarare.merged_bed
    File uCNVs_idx = merge_ultrarare.merged_bed_idx
  }
}


# Task to filter CNVs for all cohorts on a single chromosome
task filter_cnvs_singleChrom {
  String cohort
  File raw_CNVs
  String contig
  String rCNV_bucket
  Float max_freq
  String CNV_suffix
  String rCNV_docker

  command <<<
    # Copy all raw CNV data
    gsutil -m cp -r ${rCNV_bucket}/raw_data/cnv ./
    gsutil -m cp -r ${rCNV_bucket}/refs ./
    gsutil -m cp ${rCNV_bucket}/analysis/paper/data/large_segments/lit_GDs.*.bed.gz ./refs/

    # Make BED files for each cohort to be used during filtering
    # Restrict to >= 50kb to reduce size of file
    # Note: it's impossible to get >50% RO with a <50kb call given minimum size of 100kb
    while read cohort N; do
      if [ -s cnv/$cohort.raw.bed.gz ]; then
        tabix cnv/$cohort.raw.bed.gz ${contig} \
        | awk -v FS="\t" -v OFS="\t" '{ if ($3-$2>=50000) print $0 }' \
        | bgzip -c \
        > $cohort.filtered.bed.gz
        echo "$cohort"
        echo "$N"
        echo "$cohort.filtered.bed.gz"
      fi \
      | paste -s
    done < <( fgrep -v "#" /opt/rCNV2/refs/rCNV_sample_counts.txt | cut -f1-2 ) \
    | sort -nk2,2 \
    > raw_CNVs.per_cohort.txt


    # Make list of all populations to filter for allele frequency
    for pop in AFR AMR CSA EAS EUR MID OCN OTH SAS; do
      echo $pop | awk -v OFS="\n" '{ print $1"_AF", $1"_NONREF_FREQ" }'
    done > all_pop_af_fields.txt


    # Reassign location of TMPDIR to local disk, rather than boot disk
    # pybedtools can use a _ton_ of /tmp space when processing large BED files
    mkdir pbt_tmp
    mkdir pbt_tmp/tmp
    chmod 777 pbt_tmp/tmp
    ln -s pbt_tmp/tmp /tmp
    export TEMP=pbt_tmp/tmp
    export TMPDIR=pbt_tmp/tmp


    # Filter CNVs
    /opt/rCNV2/data_curation/CNV/filter_cnv_bed.py \
      --chr ${contig} \
      --minsize 100000 \
      --maxsize 20000000 \
      --maxfreq ${max_freq} \
      --recipoverlap 0.5 \
      --dist 100000 \
      --whitelist refs/lit_GDs.hc.bed.gz \
      --whitelist refs/lit_GDs.mc.bed.gz \
      --wrecip 0.75 \
      --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
      --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
      --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
      --xcov 0.5 \
      --cohorts-list raw_CNVs.per_cohort.txt \
      --vcf refs/gnomad_v2.1_sv.nonneuro.sites.vcf.gz \
      --vcf refs/CCDG_Abel_bioRxiv.sites.vcf.gz \
      --vcf refs/1000Genomes_HGSV_highCov.sites.vcf.gz \
      --vcf refs/HGDP.hg19.sites.vcf.gz \
      --vcf-af-fields "AF,$( paste -s -d, all_pop_af_fields.txt )" \
      --genome refs/GRCh37.autosomes.genome \
      --bgzip \
      ${raw_CNVs} \
      ${cohort}.${contig}.${CNV_suffix}.bed.gz
  >>>

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
    memory: "4 GB"
    disks: "local-disk 200 SSD"
    bootDiskSizeGb: "20"
  }

  output {
   File filtered_cnvs = "${cohort}.${contig}.${CNV_suffix}.bed.gz"
  }
}


# Merge beds across chromosomes
task merge_beds {
  Array[File] beds
  String prefix
  String output_bucket
  String rCNV_docker

  command <<<
    # Simple merge
    echo -e "#chr\tstart\tend\tname\tcnv\tpheno" > header.txt
    zcat ${sep=" " beds} \
    | fgrep -v "#" \
    | sort -Vk1,1 -k2,2n -k3,3n \
    | uniq \
    | cat header.txt - \
    | bgzip -c \
    > "${prefix}.bed.gz"
    tabix -f "${prefix}.bed.gz"
    # Copy to google bucket
    gsutil -m cp "${prefix}*bed.gz" ${output_bucket}/
    gsutil -m cp "${prefix}*bed.gz.tbi" ${output_bucket}/
  >>>

  runtime {
    docker: "${rCNV_docker}"
    preemptible: 1
  }

  output {
   File merged_bed = "${prefix}.bed.gz"
   File merged_bed_idx = "${prefix}.bed.gz.tbi"
   # File merged_case_bed = "${prefix}.CASE.bed.gz"
   # File merged_case_bed_idx = "${prefix}.CASE.bed.gz.tbi"
   # File merged_control_bed = "${prefix}.CTRL.bed.gz"
   # File merged_control_bed_idx = "${prefix}.CTRL.bed.gz.tbi"
  }
}
