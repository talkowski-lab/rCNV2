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
gzip gencode.v19.canonical.tsv


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
