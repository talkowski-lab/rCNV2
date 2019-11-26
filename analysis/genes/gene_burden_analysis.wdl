######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control CNV burdens per gene


workflow gene_burden_analysis {
  File phenotype_list
  File metacohort_list
  File metacohort_sample_table
  File gtf
  Int pad_controls
  String weight_mode
  Float min_cds_ovr
  Float p_cutoff
  String rCNV_bucket

  Array[Array[String]] phenotypes = read_tsv(phenotype_list)
  Array[String] cnvs = ['DEL', 'DUP', 'CNV']

  # Scatter over CNV types
  scatter ( cnv in cnvs ) {
    # Step 1.1: build null distributions of rCNV counts in all cases vs all controls
    call build_null_distrib as build_rCNV_nulls {
      input:
        hpo="HP:0000118",
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="rCNV",
        CNV=cnv,
        gtf=gtf,
        pad_controls=pad_controls,
        weight_mode=weight_mode,
        min_cds_ovr=min_cds_ovr,
        rCNV_bucket=rCNV_bucket,
        prefix="HP000018"
    }

    # Step 1.1: build null distributions of uCNV counts in all cases vs all controls
    call build_null_distrib as build_uCNV_nulls {
      input:
        hpo="HP:0000118",
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        freq_code="uCNV",
        CNV=cnv,
        gtf=gtf,
        pad_controls=pad_controls,
        weight_mode=weight_mode,
        min_cds_ovr=min_cds_ovr,
        rCNV_bucket=rCNV_bucket,
        prefix="HP000018"
    }
  }
  # Step 1.2: merge output of null distributions into master table
  call merge_null_tables as merge_rCNV_nulls {
    input:
      null_tables=build_rCNV_nulls.null_fits,
      freq_code="rCNV"
  }
  call merge_null_tables as merge_uCNV_nulls {
    input:
      null_tables=build_uCNV_nulls.null_fits,
      freq_code="uCNV"
  }


  # Scatter over phenotypes
  scatter ( pheno in phenotypes ) {
    # Step 2: run rCNV assocation tests per phenotype
    call burden_test as rCNV_burden_test {
      input:
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        null_table=merge_rCNV_nulls.null_fit_table,
        freq_code="rCNV",
        gtf=gtf,
        pad_controls=pad_controls,
        weight_mode=weight_mode,
        min_cds_ovr=min_cds_ovr,
        p_cutoff=p_cutoff,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0]
    }

    # Step 2: run uCNV assocation tests per phenotype
    call burden_test as uCNV_burden_test {
      input:
        hpo=pheno[1],
        metacohort_list=metacohort_list,
        metacohort_sample_table=metacohort_sample_table,
        null_table=merge_uCNV_nulls.null_fit_table,
        freq_code="uCNV",
        gtf=gtf,
        pad_controls=pad_controls,
        weight_mode=weight_mode,
        min_cds_ovr=min_cds_ovr,
        p_cutoff=p_cutoff,
        rCNV_bucket=rCNV_bucket,
        prefix=pheno[0]
    }

    # Step 3: run meta-analysis to combine association stats per metacohort
    # TBD
  }
}


# Build null distribution of expected case:control CNV differences
task build_null_distrib {
  String hpo
  File metacohort_list
  File metacohort_sample_table
  String freq_code
  String CNV
  File gtf
  Int pad_controls
  String weight_mode
  Float min_cds_ovr
  String rCNV_bucket
  String prefix

  command <<<
    set -e

    # Copy CNV data and constrained gene coordinates
    mkdir cleaned_cnv/
    gsutil -m cp ${rCNV_bucket}/cleaned_data/cnv/* cleaned_cnv/
    mkdir refs/
    gsutil cp ${rCNV_bucket}/analysis/analysis_refs/gencode.v19.canonical.constrained.bed.gz \
      refs/
    gsutil -m cp ${rCNV_bucket}/refs/GRCh37.*.bed.gz refs/

    # Iterate over metacohorts
    while read meta cohorts; do
      echo $meta

      # Set metacohort parameters
      cnv_bed="cleaned_cnv/$meta.${freq_code}.bed.gz"
      meta_idx=$( head -n1 "${metacohort_sample_table}" \
                  | sed 's/\t/\n/g' \
                  | awk -v meta="$meta" '{ if ($1==meta) print NR }' )
      ncase=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
               | awk -v FS="\t" -v meta_idx="$meta_idx" '{ print $meta_idx }' \
               | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
      nctrl=$( fgrep -w "HEALTHY_CONTROL" "${metacohort_sample_table}" \
               | awk -v FS="\t" -v meta_idx="$meta_idx" '{ print $meta_idx }' \
               | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )

      # Count CNVs
      /opt/rCNV2/analysis/genes/count_cnvs_per_gene.py \
        --pad-controls ${pad_controls} \
        --weight-mode ${weight_mode} \
        --min-cds-ovr ${min_cds_ovr} \
        -t ${CNV} \
        --hpo ${hpo} \
        -z \
        -o "$meta.${prefix}.${freq_code}.${CNV}.gene_burden.counts.bed.gz" \
        "$cnv_bed" \
        ${gtf}
      tabix -f "$meta.${prefix}.${freq_code}.${CNV}.gene_burden.counts.bed.gz"

      # Build null distribution
      /opt/rCNV2/analysis/genes/gene_burden_test.R \
        --pheno-table ${metacohort_sample_table} \
        --cohort-name $meta \
        --case-hpo ${hpo} \
        --build-null \
        --null-dist-plot $meta.${prefix}.${freq_code}.${CNV}.gene_burden.null_fit.jpg \
        --precision 8 \
        "$meta.${prefix}.${freq_code}.${CNV}.gene_burden.counts.bed.gz" \
        "$meta.${prefix}.${freq_code}.${CNV}.gene_burden.null_fit.txt"

    done < ${metacohort_list}

    # Merge null tables
    cat *.${prefix}.${freq_code}.${CNV}.gene_burden.null_fit.txt \
    | grep -ve '^cohort' \
    | awk -v OFS="\t" -v freq_code=${freq_code} -v CNV=${CNV} \
      '{ print freq_code, CNV, $0 }' \
    | sort -Vk3,3 \
    | cat <( head -n1 mega.${prefix}.${freq_code}.${CNV}.gene_burden.null_fit.txt \
             | sed 's/^/#freq_code\tCNV\t/g' ) \
          - \
    > ${prefix}.${freq_code}.${CNV}.gene_burden.all_null_fits.txt

    # Copy plots to rCNV bucket (note: requires permissions)
    gsutil -m cp *null_fit.jpg \
      gs://rcnv_project/analysis/gene_burden/null_fits/${freq_code}/plots/
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:99f315980611283e31e8c120f480fc477445b5250050b7f91fc26755868ab203"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    File null_fits = "${prefix}.${freq_code}.${CNV}.gene_burden.all_null_fits.txt"
  }
}


# Merge null distribution tables
task merge_null_tables {
  Array[File] null_tables
  String freq_code

  command <<<
    head -n1 ${null_tables[0]} > header.txt
    cat ${sep=' ' null_tables} \
    | fgrep -v "#" \
    | sort -Vk1,1 -k2,2V -k3,3V \
    | cat header.txt - \
    > "${freq_code}.gene_burden.all_null_fits.txt"
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:99f315980611283e31e8c120f480fc477445b5250050b7f91fc26755868ab203"
    preemptible: 1
  }

  output {
    File null_fit_table = "${freq_code}.gene_burden.all_null_fits.txt"
  }

}


# Run burden test for a single phenotype for all metacohorts
task burden_test {
  String hpo
  File metacohort_list
  File metacohort_sample_table
  File null_table
  String freq_code
  File gtf
  Int pad_controls
  String weight_mode
  Float min_cds_ovr
  Float p_cutoff
  String rCNV_bucket
  String prefix

  command <<<
    set -e

    # Copy CNV data and constrained gene coordinates
    mkdir cleaned_cnv/
    gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/* cleaned_cnv/
    gsutil -m cp -r gs://rcnv_project/cleaned_data/genes ./
    mkdir refs/
    gsutil -m cp ${rCNV_bucket}/refs/GRCh37.*.bed.gz refs/
    gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/

    # Set HPO-specific parameters
    descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
               | awk -v FS="\t" '{ print $2 }' )
    zcat refs/gencode.v19.canonical.constrained.bed.gz \
    | fgrep -wf genes/gene_lists/${prefix}.HPOdb.constrained.genes.list \
    | cat <( echo -e "#chr\tstart\tend\tgene" ) - \
    > ${prefix}.highlight_regions.bed

    # Iterate over metacohorts
    while read meta cohorts; do
      echo $meta

      # Set metacohort-specific parameters
      cnv_bed="cleaned_cnv/$meta.${freq_code}.bed.gz"
      meta_idx=$( head -n1 "${metacohort_sample_table}" \
                  | sed 's/\t/\n/g' \
                  | awk -v meta="$meta" '{ if ($1==meta) print NR }' )
      ncase=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
               | awk -v FS="\t" -v meta_idx="$meta_idx" '{ print $meta_idx }' \
               | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
      nctrl=$( fgrep -w "HEALTHY_CONTROL" "${metacohort_sample_table}" \
               | awk -v FS="\t" -v meta_idx="$meta_idx" '{ print $meta_idx }' \
               | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
      title="$descrip (${hpo})\n$ncase cases vs $nctrl controls in '$meta' cohort"

      # Iterate over CNV types
      for CNV in CNV DEL DUP; do
        echo $CNV
        # Count CNVs
        /opt/rCNV2/analysis/genes/count_cnvs_per_gene.py \
          --pad-controls ${pad_controls} \
          --weight-mode ${weight_mode} \
          --min-cds-ovr ${min_cds_ovr} \
          -t $CNV \
          --hpo ${hpo} \
          -z \
          -o "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.bed.gz" \
          "$cnv_bed" \
          ${gtf}
        tabix -f "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.bed.gz"

        # Perform burden test
        /opt/rCNV2/analysis/genes/gene_burden_test.R \
          --pheno-table ${metacohort_sample_table} \
          --cohort-name $meta \
          --cnv $CNV \
          --null-table-in ${null_table} \
          --null-model "gaussian" \
          --case-hpo ${hpo} \
          --bgzip \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden.counts.bed.gz" \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz"
        tabix -f "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz"

        # Generate Manhattan & QQ plots
        /opt/rCNV2/utils/plot_manhattan_qq.R \
          --p-col-name "phred_p" \
          --p-is-phred \
          --max-phred-p 100 \
          --cutoff ${p_cutoff} \
          --highlight-bed "${prefix}.highlight_regions.bed" \
          --highlight-name "Constrained genes associated with this phenotype" \
          --title "$title" \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz" \
          "$meta.${prefix}.${freq_code}.$CNV.gene_burden"
      done

      # Generate Miami & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --miami \
        --p-col-name "phred_p" \
        --p-is-phred \
        --max-phred-p 100 \
        --cutoff ${p_cutoff} \
        --highlight-bed "${prefix}.highlight_regions.bed" \
        --highlight-name "Constrained genes associated with this phenotype" \
        --label-prefix "DUP" \
        --highlight-bed-2 "${prefix}.highlight_regions.bed" \
        --highlight-name-2 "Constrained genes associated with this phenotype" \
        --label-prefix-2 "DEL" \
        --title "$title" \
        "$meta.${prefix}.${freq_code}.DUP.gene_burden.stats.bed.gz" \
        "$meta.${prefix}.${freq_code}.DEL.gene_burden.stats.bed.gz" \
        "$meta.${prefix}.${freq_code}.gene_burden"
    done < ${metacohort_list}

    # Copy results to output bucket
    gsutil -m cp *.gene_burden.counts.bed.gz* \
      "${rCNV_bucket}/analysis/gene_burden/${prefix}/${freq_code}/counts/"
    gsutil -m cp *.gene_burden.stats.bed.gz* \
      "${rCNV_bucket}/analysis/gene_burden/${prefix}/${freq_code}/stats/"
    gsutil -m cp *.gene_burden.*.png \
      "${rCNV_bucket}/analysis/gene_burden/${prefix}/${freq_code}/plots/"
  >>>

  runtime {
    docker: "talkowski/rcnv@sha256:99f315980611283e31e8c120f480fc477445b5250050b7f91fc26755868ab203"
    preemptible: 1
    memory: "4 GB"
    bootDiskSizeGb: "20"
  }

  output {
    Array[File] stats_beds = glob("*.gene_burden.stats.bed.gz")
    Array[File] stats_bed_idxs = glob("*.gene_burden.stats.bed.gz.tbi")
  }
}
