#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Analysis of case-control CNV burdens for canonical protein-coding genes


# Launch docker image
docker run --rm -it talkowski/rcnv


# Copy all filtered CNV data, gene coordinates, and other references 
# from the project Google Bucket (note: requires permissions)
gcloud auth login
mkdir cleaned_cnv/
gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/* cleaned_cnv/
gsutil -m cp -r gs://rcnv_project/cleaned_data/genes ./
mkdir refs/
gsutil -m cp gs://rcnv_project/refs/GRCh37.*.bed.gz refs/
gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/


# Test/dev parameters for building null distribution
hpo="HP:0012759"
meta="meta1"
freq_code="uCNV"
CNV="DEL"
metacohort_list="refs/rCNV_metacohort_list.txt"
metacohort_sample_table="refs/HPOs_by_metacohort.table.tsv"
gtf="genes/gencode.v19.canonical.gtf.gz"
weight_mode="weak"
pad_controls=50000
min_cds_ovr=0.1
prefix="HP0012759"


# Iterate over CNV types
for CNV in DEL DUP CNV; do
  echo $CNV

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
      --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
      --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
      --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
      -z \
      --verbose \
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
done

# Merge null tables
#TBD: produce uCNV.gene_burden.all_null.fits.txt




# Test/dev parameters (seizures)
hpo="HP:0001250"
prefix="HP0001250"
meta="meta2"
freq_code="uCNV"
CNV="DEL"
metacohort_list="refs/rCNV_metacohort_list.txt"
metacohort_sample_table="refs/HPOs_by_metacohort.table.tsv"
gtf="genes/gencode.v19.canonical.gtf.gz"
null_table="uCNV.gene_burden.all_null.fits.txt"
pad_controls=50000
weight_mode="weak"
min_cds_ovr=0.1
p_cutoff=0.000002587992

# Test/dev parameters (NDD)
hpo="HP:0012759"
prefix="HP0012759"
meta="meta1"
freq_code="uCNV"
CNV="DEL"
metacohort_list="refs/rCNV_metacohort_list.txt"
metacohort_sample_table="refs/HPOs_by_metacohort.table.tsv"
gtf="genes/gencode.v19.canonical.gtf.gz"
null_table="HP0000118.uCNV.DEL.gene_burden.all_null_fits.txt"
pad_controls=50000
weight_mode="weak"
min_cds_ovr=0.1
p_cutoff=0.000002587992


# Count CNVs in cases and controls per phenotype, split by metacohort and CNV type
# Iterate over phenotypes
while read prefix hpo; do

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
        --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
        --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
        --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
        -z \
        --verbose \
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
done < refs/test_phenotypes.list

