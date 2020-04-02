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
gcloud auth login


# Copy all filtered CNV data, gene coordinates, and other references 
# from the project Google Bucket (note: requires permissions)
mkdir cleaned_cnv/
gsutil -m cp -r gs://rcnv_project/cleaned_data/cnv/* cleaned_cnv/
gsutil -m cp -r gs://rcnv_project/cleaned_data/genes ./
mkdir refs/
gsutil -m cp gs://rcnv_project/refs/GRCh37.*.bed.gz refs/
gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/


# Test/dev parameters
# Params for seizures
hpo="HP:0001250"
prefix="HP0001250"
meta="meta2"
# Params for NDDs
hpo="HP:0012759"
prefix="HP0012759"
meta="meta1"
# General params
freq_code="rCNV"
CNV="DEL"
phenotype_list="refs/test_phenotypes.list"
metacohort_list="refs/rCNV_metacohort_list.txt"
metacohort_sample_table="refs/HPOs_by_metacohort.table.tsv"
gtf="genes/gencode.v19.canonical.gtf.gz"
rCNV_bucket="gs://rcnv_project"
pad_controls=0
weight_mode="weak"
min_cds_ovr_del=0.1
min_cds_ovr_dup=0.5
max_genes_per_cnv=20000
p_cutoff=0.000002896368
max_manhattan_phred_p=30
n_pheno_perms=50
meta_model_prefix="fe"
i=1


# Count CNVs in cases and controls per phenotype, split by metacohort and CNV type
# Iterate over phenotypes
while read pheno hpo; do
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
    for CNV in DEL DUP; do
      echo $CNV

      # Set CNV-specific parameters
      case "$CNV" in
        DEL)
          min_cds_ovr=${min_cds_ovr_del}
          ;;
        DUP)
          min_cds_ovr=${min_cds_ovr_dup}
          ;;
      esac

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
        "$cnv_bed" \
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

      # Generate Manhattan & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --p-col-name "fisher_phred_p" \
        --p-is-phred \
        --max-phred-p ${max_manhattan_phred_p} \
        --cutoff ${p_cutoff} \
        --highlight-bed "${prefix}.highlight_regions.bed" \
        --highlight-name "Constrained genes associated with this phenotype" \
        --label-prefix "$CNV" \
        --title "$title" \
        "$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz" \
        "$meta.${prefix}.${freq_code}.$CNV.gene_burden"
    done

    # Generate Miami & QQ plots
    /opt/rCNV2/utils/plot_manhattan_qq.R \
      --miami \
      --p-col-name "fisher_phred_p" \
      --p-is-phred \
      --max-phred-p ${max_manhattan_phred_p} \
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
done < ${phenotype_list}




# Permute CNVs in cases and controls per phenotype, split by metacohort and CNV type
# Iterate over phenotypes
while read prefix hpo; do
  for i in $( seq 1 ${n_pheno_perms} ); do

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

  done
done < ${phenotype_list}

# Gather all permutation results and compute FDR CDFs
mkdir perm_res/
while read prefix hpo; do
  gsutil -m cp \
    "${rCNV_bucket}/analysis/gene_burden/$prefix/${freq_code}/permutations/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.perm_*.bed.gz" \
    perm_res/
  for i in $( seq 1 ${n_pheno_perms} ); do
    p_idx=$( zcat perm_res/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.perm_$i.bed.gz \
             | sed -n '1p' | sed 's/\t/\n/g' | awk -v OFS="\t" '{ print $1, NR }' \
             | fgrep -w ${p_val_column_name} | cut -f2 )
    zcat perm_res/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.perm_$i.bed.gz \
    | grep -ve '^#' \
    | awk -v p_idx=$p_idx '{ print $(p_idx) }' \
    | cat <( echo "$prefix.${CNV}.$i" ) - \
    > perm_res/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.permuted_p_values.$i.txt
  done
  rm perm_res/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.perm_*.bed.gz
done < ${phenotype_list}
paste perm_res/*.gene_burden.meta_analysis.permuted_p_values.*.txt \
| gzip -c \
> ${freq_code}.${CNV}.permuted_pval_matrix.txt.gz

# Analyze p-values and compute FDR
/opt/rCNV2/analysis/sliding_windows/calc_empirical_fdr.R \
  --cnv ${CNV} \
  --fdr-target ${fdr_target} \
  --linear-fit \
  --flat-ladder \
  --plot gene_burden.${freq_code}.${CNV}.${fdr_table_suffix}_permutation_results.png \
  ${freq_code}.${CNV}.permuted_pval_matrix.txt.gz \
  ${metacohort_sample_table} \
  gene_burden.${freq_code}.${CNV}.${fdr_table_suffix}




# Run meta-analysis for each phenotype
while read prefix hpo; do

    # Get metadata for meta-analysis
    mega_idx=$( head -n1 "${metacohort_sample_table}" \
                | sed 's/\t/\n/g' \
                | awk '{ if ($1=="mega") print NR }' )
    ncase=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" -v mega_idx="$mega_idx" '{ print $mega_idx }' \
             | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
    nctrl=$( fgrep -w "HEALTHY_CONTROL" "${metacohort_sample_table}" \
             | awk -v FS="\t" -v mega_idx="$mega_idx" '{ print $mega_idx }' \
             | sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta' )
    descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
               | awk -v FS="\t" '{ print $2 }' )
    title="$descrip (${hpo})\nMeta-analysis of $ncase cases and $nctrl controls"
    DEL_p_cutoff=$( awk -v hpo=${prefix} '{ if ($1==hpo) print $2 }' \
                    gene_burden.${freq_code}.DEL.empirical_genome_wide_pval.hpo_cutoffs.tsv )
    DUP_p_cutoff=$( awk -v hpo=${prefix} '{ if ($1==hpo) print $2 }' \
                    gene_burden.${freq_code}.DUP.empirical_genome_wide_pval.hpo_cutoffs.tsv )

    # Set HPO-specific parameters
    descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
               | awk -v FS="\t" '{ print $2 }' )
    zcat refs/gencode.v19.canonical.constrained.bed.gz \
    | fgrep -wf genes/gene_lists/${prefix}.HPOdb.constrained.genes.list \
    | cat <( echo -e "#chr\tstart\tend\tgene" ) - \
    > ${prefix}.highlight_regions.bed


    # Run meta-analysis for each CNV type
    for CNV in DEL DUP; do
      # Set CNV-specific parameters
      case "$CNV" in
        DEL)
          meta_p_cutoff=$DEL_p_cutoff
          ;;
        DUP)
          meta_p_cutoff=$DUP_p_cutoff
          ;;
      esac

      # Perform meta-analysis for unweighted CNVs
      while read meta cohorts; do
        echo -e "$meta\t$meta.${prefix}.${freq_code}.$CNV.gene_burden.stats.bed.gz"
      done < <( fgrep -v mega ${metacohort_list} ) \
      > ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.input.txt
      /opt/rCNV2/analysis/genes/gene_meta_analysis.R \
        --or-corplot ${prefix}.${freq_code}.$CNV.gene_burden.or_corplot_grid.jpg \
        --model ${meta_model_prefix} \
        --p-is-phred \
        --spa \
        ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.input.txt \
        ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed
      bgzip -f ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed
      tabix -f ${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed.gz

      # Generate Manhattan & QQ plots
      /opt/rCNV2/utils/plot_manhattan_qq.R \
        --p-col-name "meta_phred_p" \
        --p-is-phred \
        --max-phred-p ${max_manhattan_phred_p} \
        --cutoff $meta_p_cutoff \
        --highlight-bed "${prefix}.highlight_regions.bed" \
        --highlight-name "Constrained genes associated with this phenotype" \
        --label-prefix "$CNV" \
        --title "$title" \
        "${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed.gz" \
        "${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis"
    done

    # Generate Miami & QQ plots
    /opt/rCNV2/utils/plot_manhattan_qq.R \
      --miami \
      --p-col-name "meta_phred_p" \
      --p-is-phred \
      --max-phred-p ${max_manhattan_phred_p} \
      --cutoff $DUP_p_cutoff \
      --highlight-bed "${prefix}.highlight_regions.bed" \
      --highlight-name "Constrained genes associated with this phenotype" \
      --label-prefix "DUP" \
      --cutoff-2 $DEL_p_cutoff \
      --highlight-bed-2 "${prefix}.highlight_regions.bed" \
      --highlight-name-2 "Constrained genes associated with this phenotype" \
      --label-prefix-2 "DEL" \
      --title "$title" \
      "${prefix}.${freq_code}.DUP.gene_burden.meta_analysis.stats.bed.gz" \
      "${prefix}.${freq_code}.DEL.gene_burden.meta_analysis.stats.bed.gz" \
      "${prefix}.${freq_code}.gene_burden.meta_analysis"

  done
done < refs/test_phenotypes.list

# Collapse all meta-analysis p-values into single matrix for visualizing calibration
mkdir meta_res/
# Download data
while read prefix hpo; do
  for CNV in DEL DUP; do
    echo -e "${rCNV_bucket}/analysis/gene_burden/${prefix}/${freq_code}/stats/${prefix}.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed.gz"
  done
done < ${phenotype_list} \
| gsutil -m cp -I meta_res/
# Primary p-values
p_val_column_name="meta_phred_p"
while read prefix hpo; do
  echo -e "$prefix\n\n"
  for CNV in DEL DUP; do
      stats=meta_res/$prefix.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed.gz
      if [ -e $stats ]; then
        p_idx=$( zcat $stats | sed -n '1p' | sed 's/\t/\n/g' \
                 | awk -v OFS="\t" '{ print $1, NR }' \
                 | fgrep -w ${p_val_column_name} | cut -f2 )
        zcat $stats | grep -ve '^#' \
        | awk -v p_idx=$p_idx '{ print $(p_idx) }' \
        | cat <( echo "$prefix.$CNV" ) - \
        > meta_res/$prefix.${freq_code}.$CNV.gene_burden.meta_analysis.p_values.txt
      fi
  done
done < ${phenotype_list}
paste meta_res/*.${freq_code}.$CNV.gene_burden.meta_analysis.p_values.txt \
| gzip -c \
> ${freq_code}.observed_pval_matrix.txt.gz
# Secondary p-values
p_val_column_name="meta_phred_p_secondary"
while read prefix hpo; do
  echo -e "$prefix\n\n"
  for CNV in DEL DUP; do
      stats=meta_res/$prefix.${freq_code}.$CNV.gene_burden.meta_analysis.stats.bed.gz
      if [ -e $stats ]; then
        p_idx=$( zcat $stats | sed -n '1p' | sed 's/\t/\n/g' \
                 | awk -v OFS="\t" '{ print $1, NR }' \
                 | fgrep -w ${p_val_column_name} | cut -f2 )
        zcat $stats | grep -ve '^#' \
        | awk -v p_idx=$p_idx '{ print $(p_idx) }' \
        | cat <( echo "$prefix.$CNV" ) - \
        > meta_res/$prefix.${freq_code}.$CNV.gene_burden.meta_analysis.secondary_p_values.txt
      fi
  done
done < ${phenotype_list}
paste meta_res/*.${freq_code}.$CNV.gene_burden.meta_analysis.secondary_p_values.txt \
| gzip -c \
> ${freq_code}.observed_pval_matrix_secondary.txt.gz
# DEV NOTE: these p-values can be visualized with plot_gene_burden_meta_analysis_p_values.R,
#           which is currently just a code snippet referencing local filepaths




# Fine-map significant gene blocks
# Test/dev parameters
freq_code="rCNV"
CNV="DEL"
rCNV_bucket="gs://rcnv_project"
phenotype_list="test_phenotypes.list"
metacohort_sample_table="HPOs_by_metacohort.table.tsv"
meta_p_cutoffs_tsv="gene_burden.rCNV.DEL.empirical_genome_wide_pval.hpo_cutoffs.tsv"
gene_features="gencode.v19.canonical.pext_filtered.all_features.eigenfeatures.bed.gz"
finemap_output_label="all_features"
meta_secondary_p_cutoff=0.05
meta_nominal_cohorts_cutoff=2
finemap_elnet_alpha=0.1
finemap_elnet_l1_l2_mix=1
finemap_distance=1000000

# Copy association stats & gene lists from the project Google Bucket (note: requires permissions)
mkdir stats
gsutil -m cp \
  ${rCNV_bucket}/analysis/gene_burden/**.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz \
  stats/
gsutil -m cp -r \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists \
  ./

# Write tsv input
while read prefix hpo; do
  for wrapper in 1; do
    echo "$hpo"
    echo "stats/$prefix.${freq_code}.${CNV}.gene_burden.meta_analysis.stats.bed.gz"
    awk -v x=$prefix -v FS="\t" '{ if ($1==x) print $2 }' \
      gene_burden.${freq_code}.${CNV}.empirical_genome_wide_pval.hpo_cutoffs.tsv
  done | paste -s
done < ${phenotype_list} \
> ${freq_code}.${CNV}.gene_fine_mapping.stats_input.tsv

# Make lists of a priori "true causal" genes for null variance estimation
case ${CNV} in
  "DEL")
    for wrapper in 1; do
      echo "gene_lists/DDG2P.hmc_lof.genes.list"
      echo "gene_lists/DDG2P.hmc_other.genes.list"
      echo "gene_lists/ClinGen.hmc_haploinsufficient.genes.list"
      echo "gene_lists/HP0000118.HPOdb.constrained.genes.list"
      echo "gene_lists/gnomad.v2.1.1.lof_constrained.genes.list"
    done > known_causal_gene_lists.tsv
    ;;
  "DUP")
    for wrapper in 1; do
      echo "gene_lists/DDG2P.all_gof.genes.list"
      echo "gene_lists/DDG2P.hmc_other.genes.list"
      echo "gene_lists/ClinGen.all_triplosensitive.genes.list"
      echo "gene_lists/HP0000118.HPOdb.constrained.genes.list"
      echo "gene_lists/gnomad.v2.1.1.lof_constrained.genes.list"
    done > known_causal_gene_lists.tsv
    ;;
esac

# Run functional fine-mapping procedure
/opt/rCNV2/analysis/genes/finemap_genes.py \
  --secondary-p-cutoff ${meta_secondary_p_cutoff} \
  --min-nominal ${meta_nominal_cohorts_cutoff} \
  --secondary-or-nominal \
  --regularization-alpha ${finemap_elnet_alpha} \
  --regularization-l1-l2-mix ${finemap_elnet_l1_l2_mix} \
  --distance ${finemap_distance} \
  --known-causal-gene-lists known_causal_gene_lists.tsv \
  --outfile ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.${finemap_output_label}.tsv \
  --all-genes-outfile ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.${finemap_output_label}.all_genes_from_blocks.tsv \
  --naive-outfile ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.naive_priors.${finemap_output_label}.tsv \
  --genetic-outfile ${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genetics_only.${finemap_output_label}.tsv \
  --coeffs-out ${freq_code}.${CNV}.gene_fine_mapping.logit_coeffs.${finemap_output_label}.tsv \
  ${freq_code}.${CNV}.gene_fine_mapping.stats_input.tsv \
  ${gene_features}




# Plot results of fine-mapping
# Test/dev parameters
freq_code="rCNV"
CNV="DEL"
phenotype_list="test_phenotypes.list"
raw_features_genomic="gencode.v19.canonical.pext_filtered.genomic_features.bed.gz"
raw_features_expression="gencode.v19.canonical.pext_filtered.expression_features.bed.gz"
raw_features_constraint="gencode.v19.canonical.pext_filtered.constraint_features.bed.gz"
raw_features_merged="gencode.v19.canonical.pext_filtered.all_features.bed.gz"
rCNV_bucket="gs://rcnv_project"

# Copy all fine-mapped gene lists
mkdir finemap_stats/
gsutil -m cp \
  ${rCNV_bucket}/analysis/gene_burden/fine_mapping/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.*.tsv \
  finemap_stats/

# Make input tsvs
for wrapper in 1; do
  echo -e "Prior\tgrey70\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.naive_priors.genomic_features.tsv"
  echo -e "Posterior\t'#264653'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genetics_only.genomic_features.tsv"
  echo -e "Genomic features\t'#E76F51'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genomic_features.tsv"
  echo -e "Gene expression\t'#E9C46A'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.expression_features.tsv"
  echo -e "Gene constraint\t'#F4A261'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.constraint_features.tsv"
  echo -e "Full model\t'#2A9D8F'\t1\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.merged_features.tsv"
done > finemap_roc_input.tsv
for wrapper in 1; do
  echo -e "Genomic features\t'#E76F51'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.genomic_features.all_genes_from_blocks.tsv\t${raw_features_genomic}"
  echo -e "Gene expression\t'#E9C46A'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.expression_features.all_genes_from_blocks.tsv\t${raw_features_expression}"
  echo -e "Gene constraint\t'#F4A261'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.constraint_features.all_genes_from_blocks.tsv\t${raw_features_constraint}"
  echo -e "Full model\t'#2A9D8F'\tfinemap_stats/${freq_code}.${CNV}.gene_fine_mapping.gene_stats.merged_features.all_genes_from_blocks.tsv\t${raw_features_merged}"
done > finemap_feature_cor_input.tsv

# Make all gene truth sets
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/genes/gene_lists ./

# HPO-associated
while read prefix hpo; do
  awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }' \
    gene_lists/$prefix.HPOdb.genes.list
done < ${phenotype_list} \
| sort -Vk1,1 -k2,2V | uniq \
| cat <( echo -e "#HPO\tgene" ) - \
> hpo_truth_set.tsv

# Make CNV type-dependent truth sets
case ${CNV} in
  "DEL")
    # Union (ClinGen HI + DDG2P dominant LoF)
    while read prefix hpo; do
      cat gene_lists/ClinGen.hmc_haploinsufficient.genes.list \
          gene_lists/DDG2P.hmc_lof.genes.list \
      | awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }'
    done < ${phenotype_list} \
    | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > union_truth_set.lof.tsv

    # Intersection (ClinGen HI + DDG2P dominant LoF)
    while read prefix hpo; do
      fgrep -wf \
        gene_lists/ClinGen.hmc_haploinsufficient.genes.list \
        gene_lists/DDG2P.hmc_lof.genes.list \
      | awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }'
    done < ${phenotype_list} \
    | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > intersection_truth_set.lof.tsv

    # ClinGen HI alone
    while read prefix hpo; do
      awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }' \
        gene_lists/ClinGen.all_haploinsufficient.genes.list
    done < ${phenotype_list} \
    | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > clingen_truth_set.lof.tsv

    # DDG2P lof + other alone
    while read prefix hpo; do
      cat gene_lists/DDG2P.all_lof.genes.list \
          gene_lists/DDG2P.all_other.genes.list \
      | awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }'
    done < ${phenotype_list} \
    | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > ddg2p_truth_set.lof.tsv

    # Combine all truth sets into master union
    cat union_truth_set.lof.tsv \
        clingen_truth_set.lof.tsv \
        ddg2p_truth_set.lof.tsv \
        hpo_truth_set.tsv \
    | fgrep -v "#" | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > master_union_truth_set.lof.tsv

    # Write truth set input tsv
    for wrapper in 1; do
      echo -e "Union of all truth sets\tmaster_union_truth_set.lof.tsv"
      echo -e "ClinGen HI & DECIPHER LoF (union)\tunion_truth_set.lof.tsv"
      echo -e "ClinGen HI & DECIPHER LoF (intersection)\tintersection_truth_set.lof.tsv"
      echo -e "ClinGen HI (any confidence)\tclingen_truth_set.lof.tsv"
      echo -e "DECIPHER dominant LoF/unk. (any confidence)\tddg2p_truth_set.lof.tsv"
      echo -e "HPO-matched disease genes\thpo_truth_set.tsv"
    done > finemap_roc_truth_sets.tsv
    ;;

  "DUP")
    # Union (ClinGen HI + DDG2P dominant CG)
    while read prefix hpo; do
      cat gene_lists/ClinGen.hmc_triplosensitive.genes.list \
          gene_lists/DDG2P.hmc_gof.genes.list \
      | awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }'
    done < ${phenotype_list} \
    | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > union_truth_set.gof.tsv

    # ClinGen triplo alone
    while read prefix hpo; do
      awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }' \
        gene_lists/ClinGen.all_triplosensitive.genes.list
    done < ${phenotype_list} \
    | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > clingen_truth_set.triplo.tsv

    # DDG2P gof + other alone
    while read prefix hpo; do
      cat gene_lists/DDG2P.all_gof.genes.list \
          gene_lists/DDG2P.all_other.genes.list \
      | awk -v OFS="\t" -v hpo=$hpo '{ print hpo, $1 }'
    done < ${phenotype_list} \
    | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > ddg2p_truth_set.gof.tsv

    # Combine all truth sets into master union
    cat union_truth_set.gof.tsv \
        clingen_truth_set.gof.tsv \
        ddg2p_truth_set.gof.tsv \
        hpo_truth_set.tsv \
    | fgrep -v "#" | sort -Vk1,1 -k2,2V | uniq \
    | cat <( echo -e "#HPO\tgene" ) - \
    > master_union_truth_set.gof.tsv

    # Write truth set input tsv
    for wrapper in 1; do
      echo -e "Union of all truth sets\tmaster_union_truth_set.gof.tsv"
      echo -e "ClinGen TS & DECIPHER GoF (union)\tunion_truth_set.gof.tsv"
      echo -e "ClinGen TS (any confidence)\tclingen_truth_set.triplo.tsv"
      echo -e "DECIPHER dominant GoF/unk. (any confidence)\tddg2p_truth_set.gof.tsv"
      echo -e "HPO-matched disease genes\thpo_truth_set.tsv"
    done > finemap_roc_truth_sets.tsv
    ;;

esac

# Plot finemapping QC & feature correlations
mkdir ${freq_code}_${CNV}_finemap_plots/
/opt/rCNV2/analysis/genes/plot_finemap_results.R \
  finemap_roc_input.tsv \
  finemap_roc_truth_sets.tsv \
  ${freq_code}_${CNV}_finemap_plots/${freq_code}.${CNV}.finemap_results
/opt/rCNV2/analysis/genes/plot_finemap_coefficients.R \
  finemap_feature_cor_input.tsv \
  ${freq_code}_${CNV}_finemap_plots/${freq_code}.${CNV}.finemap_feature_cors

# Compress results
tar -czvf ${freq_code}_${CNV}_finemap_plots.tgz ${freq_code}_${CNV}_finemap_plots

