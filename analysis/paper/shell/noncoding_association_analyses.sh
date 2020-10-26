#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Master code block for secondary analyses of noncoding association for rCNV2 manuscript


# Launch docker image & authenticate GCP credentials
docker run --rm -it talkowski/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d1"


# Download necessary data (note: requires permissions)
gsutil -m cp -r ${rCNV_bucket}/cleaned_data/cnv ./
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/genome_annotations/rCNV.burden_stats.tsv.gz \
  ${rCNV_bucket}/cleaned_data/genome_annotations/rCNV.crbs.bed.gz \
  ./
mkdir refs/
gsutil -m cp -r \
  ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* \
  ${rCNV_bucket}/refs/GRCh37.* \
  ${rCNV_bucket}/analysis/analysis_refs/test_phenotypes.list \
  ${rCNV_bucket}/analysis/analysis_refs/HPOs_by_metacohort.table.tsv \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists \
  ${rCNV_bucket}/refs/REP_state_manifest.tsv \
  refs/


# Plot observed effect size distributions split by gene set membership
# Used to justify inclusion of some coding effects in noncoding association test
fgrep -wvf \
  refs/gene_lists/HP0000118.HPOdb.genes.list \
  refs/gene_lists/gnomad.v2.1.1.likely_unconstrained.genes.list \
| sort -V | uniq \
> loose_noncoding_whitelist.genes.list
mkdir genes_per_cnv
for CNV in DEL DUP; do
  # Annotate all CNVs based on any exon overlap
  for i in $( seq 1 4 ); do
    /opt/rCNV2/analysis/genes/count_cnvs_per_gene.py \
      cnv/meta${i}.rCNV.bed.gz \
      refs/gencode.v19.canonical.pext_filtered.gtf.gz \
      --min-cds-ovr "10e-10" \
      -t ${CNV} \
      --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
      --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
      --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
      -o /tmp/junk.bed.gz \
      --bgzip \
      --cnvs-out /dev/stdout \
    | cut -f4,6-7,9 | gzip -c \
    > genes_per_cnv/rCNV.${CNV}.meta${i}.genes_per_cnv.tsv.gz
    echo -e "meta${i}\tgenes_per_cnv/rCNV.${CNV}.meta${i}.genes_per_cnv.tsv.gz"
  done \
  > genes_per_cnv.${CNV}.input.tsv
  # Summarize CNV counts by gene list per phenotype
  /opt/rCNV2/analysis/paper/scripts/noncoding_association/summarize_ncCNV_counts.py \
    --genes-per-cnv genes_per_cnv.${CNV}.input.tsv \
    --hpos <( cut -f2 refs/test_phenotypes.list ) \
    --unconstrained-genes loose_noncoding_whitelist.genes.list \
    --summary-counts unconstrained_cnv_counts.${CNV}.tsv.gz \
    --gzip
done
/opt/rCNV2/analysis/paper/plot/noncoding_association/plot_unconstrained_effect_sizes.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  unconstrained_cnv_counts.DEL.tsv.gz \
  unconstrained_cnv_counts.DUP.tsv.gz \
  refs/HPOs_by_metacohort.table.tsv \
  ${prefix}



# Plot annotation track distributions & stats
/opt/rCNV2/analysis/paper/plot/noncoding_association/plot_track_stats.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.burden_stats.tsv.gz \
  ${prefix}


# Volcano plots of DEL & DUP track-level burden tests
/opt/rCNV2/analysis/paper/plot/noncoding_association/plot_volcanos.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.burden_stats.tsv.gz \
  ${prefix}


# Plot ChromHMM enrichments as positive controls
/opt/rCNV2/analysis/paper/plot/noncoding_association/chromhmm_enrichments.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.burden_stats.tsv.gz \
  refs/REP_state_manifest.tsv \
  ${prefix}


# Plot CRB distributions & stats
/opt/rCNV2/analysis/paper/plot/noncoding_association/plot_crb_stats.R \
  --rcnv-config /opt/rCNV2/config/rCNV2_rscript_config.R \
  rCNV.crbs.bed.gz \
  ${prefix}


# Copy all plots to final gs:// directory
gsutil -m cp -r \
  ${prefix}.cnv_lnORs_by_genic_context.pdf \
  ${prefix}.track_stats.*.pdf \
  ${prefix}.*volcano*pdf \
  ${prefix}.chromhmm_enrichment.*.pdf \
  ${prefix}.crb_stats.*.pdf \
  ${rCNV_bucket}/analysis/paper/plots/noncoding_association/
