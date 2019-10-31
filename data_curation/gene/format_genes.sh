#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Curate & annotate genes used in subsequent analyses


# Launch docker image
docker run --rm -it talkowski/rcnv
gcloud auth login


# Download all gene data from Gencode v19
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.pc_translations.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.pc_transcripts.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.metadata.Transcript_source.gz


# Build table of canonical protein-coding transcripts from Gencode
/opt/rCNV2/data_curation/gene/get_canonical_transcripts.py \
  gencode.v19.annotation.gtf.gz \
  gencode.v19.pc_translations.fa.gz \
  gencode.v19.pc_transcripts.fa.gz \
  gencode.v19.metadata.Transcript_source.gz \
  gencode.v19.canonical.tsv
gzip -f gencode.v19.canonical.tsv


# Subset GTF to autosomal canonical protein-coding transcripts
cat \
  <( zcat gencode.v19.canonical.tsv.gz \
     | cut -f2 | sed '1d' \
     | fgrep -wf - <( zcat gencode.v19.annotation.gtf.gz ) ) \
  <( zcat gencode.v19.canonical.tsv.gz \
     | cut -f1 | sed '1d' \
     | fgrep -wf - <( zcat gencode.v19.annotation.gtf.gz ) \
     | awk '($3=="gene")' ) \
| sed -e 's/^chr//' \
| grep -e '^[1-9]' \
| sort -k10,10 \
| /opt/rCNV2/data_curation/gene/filter_UTRs.py stdin stdout \
| sort -k1,1V -k4,4n \
| bgzip -c \
> gencode.v19.canonical.gtf.gz


# Copy canonical gene metadata to rCNV bucket (note: requires permissions)
gsutil cp gencode.v19.canonical.*.gz gs://rcnv_project/cleaned_data/genes/


# Generate BED file of pLoF constrained genes from gnomAD
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
pli_idx=$( zcat gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
           | head -n1 | sed 's/\t/\n/g' \
           | awk '{ if ($1=="pLI") print NR }' )
zcat gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
| awk -v col_idx=${pli_idx} -v FS="\t" '{ if ($col_idx >= 0.95) print $1 }' \
| sort -Vk1,1 \
> gnomad.v2.1.1.lof_constrained.genes.list
awk '{ print "gene_name \""$1"\";" }' \
gnomad.v2.1.1.lof_constrained.genes.list \
| fgrep -wf - <( zcat gencode.v19.canonical.gtf.gz ) \
| awk -v FS="\t" '{ if ($3=="transcript") print }' \
> gencode.v19.canonical.constrained.gtf
while read gene; do
  fgrep -w $gene gencode.v19.canonical.constrained.gtf \
  | awk -v OFS="\t" -v gene=$gene '{ print $1, $4, $5, gene }'
done < gnomad.v2.1.1.lof_constrained.genes.list \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
| cat <( echo -e "#chr\tstart\tend\tgene" ) - \
| bgzip -c \
> gencode.v19.canonical.constrained.bed.gz


# Copy constrained gene BED file to rCNV bucket (note: requires permissions)
gsutil cp gencode.v19.canonical.constrained.bed.gz \
  gs://rcnv_project/analysis/analysis_refs/
gsutil cp 

