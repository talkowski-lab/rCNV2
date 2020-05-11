#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Helper script to build comparison table of predicted NAHR-mediated CNVs
# Note: assumes running in Docker image with local files per curate_known_gds.sh


# Download segdups file and GTF, if needed
if ! [ -e genomicSuperDups.txt.gz ]; then
  wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/genomicSuperDups.txt.gz
fi


# Find pairs of segdups with the following criteria:
# 1. Same chromosome
# 2. Distance ≥ 100kb & ≤ 10Mb
# 3. Both segdups ≥ 1kb
# 4. Strict homology ≥ 95%
# 5. Direct orentation of repeats (i.e., same strand)
# 6. Total intervening sequence ≤30% covered by rCNV blacklist 
#    (including segdup/simple repeat/satellites/Nmask/somatic hypermutable)
# 7. At least 100kb of intervening sequence after subtracting blacklist regions
zcat genomicSuperDups.txt.gz \
| awk -v FS="\t" -v OFS="\t" \
  '{ if ($2==$8 && $7=="+" && $4-$3>=1000 && $10-$9>=1000 && \
         $9-$4>=100000 && $9-$4<=10000000 && $27>=0.95) \
     print $2, $4, $9, $2, $3, $4, $8, $9, $10 }' \
| sed 's/chr//g' \
| grep -e '^[0-9]' \
| bedtools coverage -a - \
  -b <( zcat refs/GRCh37.Nmask.autosomes.bed.gz \
             refs/GRCh37.somatic_hypermutable_sites.bed.gz \
             refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
        | cut -f1-3 | sort -Vk1,1 -k2,2n -k3,3n | bedtools merge -i - ) \
| awk -v FS="\t" -v OFS="\t" \
  '{ if ($NF<=0.3 && $12-$11>=100000) print $4, $5, $6, $7, $8, $9 }' \
| sort -Vk1,1 -k2,2n -k3,3n -k5,5n -k6,6n \
| awk -v FS="\t" -v OFS="\t" '{ print $0, "candidate_segdup_pair_"NR }' \
> candidate_segdup_pairs.bedpe
# Cluster candidate segdup pairs with bedtools pairtopair
# To cluster, candidate pairs must:
# 1. Have both ends within 1Mb of each other, and
# 2. Have >50% reciprocal overlap of intervening sequence
# Also pre-formats for gene annotation
/opt/rCNV2/data_curation/other/collapse_segdup_pairs.py \
  --distance 1000000 \
  --recip 0.5 \
  candidate_segdup_pairs.bedpe \
| sort -Vk1,1 -k2,2n -k3,3n \
| awk -v FS="\t" -v OFS="\t" '{ print $0, "nahr_region_"NR, $1":"$2"-"$3, "." }' \
| cat <( echo -e "#chr\tstart\tend\tnahr_id\tcoords\tdummy" ) - \
> clustered_nahr_regions.no_genes.bed


# Annotate vs genes
/opt/rCNV2/analysis/sliding_windows/get_genes_per_region.py \
    -o clustered_nahr_regions.w_genes.bed \
    clustered_nahr_regions.no_genes.bed \
    refs/gencode.v19.canonical.pext_filtered.gtf.gz
cut -f1-4,7-8 clustered_nahr_regions.w_genes.bed \
| bgzip -c > clustered_nahr_regions.bed.gz

