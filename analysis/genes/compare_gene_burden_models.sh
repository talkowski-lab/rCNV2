#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Determine best gene meta-analysis parameters based on distribution of lambda ratios

# Note: assumes mixture of parameters considered during Feb 2020 methods development
#       at a later date, some of these parameters/frequencies/models may be obsolete


# Launch docker image
docker run --rm -it talkowski/rcnv


# Download basic refs
gcloud auth login
gsutil -m cp -r gs://rcnv_project/cleaned_data/genes ./
mkdir refs/
gsutil -m cp gs://rcnv_project/analysis/analysis_refs/* refs/


# Set general params
metacohort_list="refs/rCNV_metacohort_list.txt"
metacohort_sample_table="refs/HPOs_by_metacohort.table.tsv"
gtf="genes/gencode.v19.canonical.gtf.gz"
meta_p_cutoff=0.000002896368


# Iterate over each phenotype, CNV type, frequency bin, and analysis model, and copy data
while read prefix hpo; do
  for CNV in DEL DUP; do
    for freq in rCNV vCNV uCNV; do
      for model in meta_analysis weighted_meta_analysis; do
        # Echo stats file
        echo -e "gs://rcnv_project/analysis/gene_burden/$prefix/$freq/stats/$prefix.$freq.$CNV.gene_burden.$model.stats.bed.gz"
      done
    done
  done
done < refs/test_phenotypes.list \
| gsutil -m cp -I ./


# Iterate over each phenotype, CNV type, frequency bin, and analysis model, and compute lambdas
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

  # Set HPO-specific parameters
  descrip=$( fgrep -w "${hpo}" "${metacohort_sample_table}" \
             | awk -v FS="\t" '{ print $2 }' )
  zcat refs/gencode.v19.canonical.constrained.bed.gz \
  | fgrep -wf genes/gene_lists/${prefix}.HPOdb.constrained.genes.list \
  | cat <( echo -e "#chr\tstart\tend\tgene" ) - \
  > ${prefix}.highlight_regions.bed

  for CNV in DEL DUP; do
    for freq in rCNV vCNV uCNV; do
      for model in meta_analysis weighted_meta_analysis; do
        # Compute lambdas
        /opt/rCNV2/utils/plot_manhattan_qq.R \
          --p-col-name "meta_phred_p" \
          --p-is-phred \
          --cutoff ${meta_p_cutoff} \
          --highlight-bed "${prefix}.highlight_regions.bed" \
          --highlight-name "Constrained genes associated with this phenotype" \
          --label-prefix "$CNV" \
          --title "$title" \
          --echo-lambdas \
          --qq-only \
          $prefix.$freq.$CNV.gene_burden.$model.stats.bed.gz \
          qq_out \
        | fgrep -v "Printing" \
        | paste <( echo -e "$prefix\t$CNV\t$freq\t$model" ) -
      done
    done
  done
done < refs/test_phenotypes.list \
| cat <( echo -e "hpo\tcnv\tfreq\tmodel\tla\tlc\tratio" ) - \
> rCNV2.gene_burden_meta_analysis.lambda_ratios.tsv


