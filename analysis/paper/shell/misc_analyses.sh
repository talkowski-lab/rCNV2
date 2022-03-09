#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Miscellaneous helper analyses for rCNV2 manuscript


# Launch docker image & authenticate GCP credentials
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv
gcloud auth login


# Set global parameters
export rCNV_bucket="gs://rcnv_project"
export prefix="rCNV2_analysis_d2"
export finemap_dist=1000000


# Download necessary data (note: requires permissions)
gsutil -m cp -r \
  ${rCNV_bucket}/cleaned_data/cnv \
  ${rCNV_bucket}/analysis/analysis_refs/HPOs_by_metacohort.table.tsv \
  ${rCNV_bucket}/results/gene_association/* \
  ${rCNV_bucket}/results/segment_association/* \
  ${rCNV_bucket}/analysis/gene_burden/fine_mapping/rCNV.*.gene_fine_mapping.gene_stats.merged_no_variation_features.tsv \
  ${rCNV_bucket}/analysis/gene_burden/fine_mapping/rCNV.*.gene_fine_mapping.credible_sets_per_hpo.merged_no_variation_features.bed \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists \
  ${rCNV_bucket}/analysis/gene_scoring/gene_lists \
  ${rCNV_bucket}/analysis/paper/data/global/${prefix}.global_burden_stats.tsv.gz \
  ./
mkdir refs
gsutil -m cp \
  ${rCNV_bucket}/analysis/analysis_refs/* \
  refs/

# Count total number of CNVs
for CNV in DEL DUP; do
  zcat cnv/mega.rCNV.bed.gz | fgrep -w $CNV | wc -l
done


# Get total number of samples, cases. and controls
n_control=$( fgrep -w HEALTHY_CONTROL HPOs_by_metacohort.table.tsv | cut -f3 )
n_case=$( fgrep -w "HP:0000118" HPOs_by_metacohort.table.tsv | cut -f3 )
echo -e "$n_case\t$n_control" | awk -v OFS="\n" '{ print $1+$2, $1, $2 }'


# Compute total number of significant loci
for CNV in DEL DUP; do
  zcat rCNV.final_segments.loci.bed.gz \
       rCNV.final_genes.credible_sets.bed.gz \
  | fgrep -w $CNV \
  | cut -f1-3 \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  | wc -l
done
# (As above, but split by significance)
for CNV in DEL DUP; do
  # Genome- or exome-wide
  zcat rCNV.final_segments.loci.bed.gz \
       rCNV.final_genes.credible_sets.bed.gz \
  | fgrep -w $CNV | grep -e '[genome|exome]_wide' \
  | cut -f1-3 \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  > gw_ew.$CNV.bed
  cat gw_ew.$CNV.bed | wc -l
  # FDR
  zcat rCNV.final_segments.loci.bed.gz \
       rCNV.final_genes.credible_sets.bed.gz \
  | fgrep -w $CNV | fgrep -w FDR \
  | cut -f1-3 \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  | bedtools intersect -v -a - -b gw_ew.$CNV.bed \
  | wc -l
  rm gw_ew.$CNV.bed
done


# Preformat genic associations to split out per HPO
while read chrom start end rid cnv sig hpos; do
  while read hpo; do
    echo -e "$chrom\t$start\t$end\t${rid}_$hpo\t$cnv\t$sig\t$hpo"
  done < <( echo $hpos | sed 's/;/\n/g' )
done < <( zcat rCNV.final_genes.credible_sets.bed.gz | cut -f1-6,21 | fgrep -v "#" ) \
| cat <( zcat rCNV.final_genes.credible_sets.bed.gz | head -n1 | cut -f1-6,21 ) - \
| bgzip -c > rCNV.final_genes.credible_sets.split_by_hpo.bed.gz
cat rCNV.*.gene_fine_mapping.credible_sets_per_hpo.merged_no_variation_features.bed \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V -k5,5V | fgrep -v "#" | uniq | bgzip -c \
> rCNV.final_genes.credible_sets_hpo_specific.bed.gz


# Compute total number of significant associations
for CNV in DEL DUP; do
  cat \
    <( zcat rCNV.final_segments.associations.bed.gz | fgrep -w $CNV | awk -v FS="\t" -v OFS="\t" '{ print $6"_"$1, $2, $3 }' ) \
    <( zcat rCNV.final_genes.credible_sets.split_by_hpo.bed.gz | fgrep -w $CNV | awk -v FS="\t" -v OFS="\t" '{ print $7"_"$1, $2, $3 }' ) \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  | wc -l
done
# (As above, but split by significance)
for CNV in DEL DUP; do
  # Genome- or exome-wide
  zcat rCNV.final_segments.associations.bed.gz \
       rCNV.final_genes.credible_sets_hpo_specific.bed.gz \
  | fgrep -w $CNV | grep -e '[genome|exome]_wide' \
  | awk -v FS="\t" -v OFS="\t" '{ print $6"_"$1, $2, $3 }' \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  > gw_ew.$CNV.bed
  cat gw_ew.$CNV.bed | wc -l
  # FDR
  zcat rCNV.final_segments.associations.bed.gz \
       rCNV.final_genes.credible_sets_hpo_specific.bed.gz \
  | fgrep -w $CNV | fgrep -w FDR \
  | awk -v FS="\t" -v OFS="\t" '{ print $6"_"$1, $2, $3 }' \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bedtools merge -i - \
  | bedtools intersect -v -a - -b gw_ew.$CNV.bed \
  | wc -l
  rm gw_ew.$CNV.bed
done


# Calculate number of genic associations not overlapping significant large segment
# Note that this is slightly imperfect because the gene credible sets get their coordinates refined
# during joint fine mapping, so there may be a few situations where the overlap is incorrect
# It may be more accurate to just subtract the number of large segment associations from the total
# number of associations calculated above (as that is more accurate)
for CNV in DEL DUP; do
  zcat rCNV.final_genes.credible_sets_hpo_specific.bed.gz | fgrep -w $CNV \
  | awk -v OFS="\t" '{ print $6"_"$1, $2, $3 }' \
  | bedtools intersect -v -a - \
    -b <( zcat rCNV.final_segments.associations.bed.gz | fgrep -w $CNV \
          | awk -v OFS="\t" '{ print $6"_"$1, $2, $3 }' ) \
  | wc -l
done


# Get number of candidate novel TS disease genes
zcat rCNV.final_genes.genes.bed.gz | fgrep -w DUP | cut -f4 \
| fgrep -wvf gene_lists/DDG2P.all_gof.genes.list \
| fgrep -wvf gene_lists/ClinGen.all_triplosensitive.genes.list \
| sort | uniq | wc -l


# Get total number of fine-mapped candidate driver genes by CNV type
zcat rCNV.final_genes.genes.bed.gz \
| fgrep -v "#" | cut -f4 | sort | uniq | wc -l 


# Get gold-standard dosage-sensitive genes & PIPs to higlight in bottom of graphical abstract
cat \
  <( zcat rCNV.final_genes.genes.bed.gz | fgrep -w DEL \
     | fgrep -wf  gene_lists/gold_standard.haploinsufficient.genes.list ) \
  <( zcat rCNV.final_genes.genes.bed.gz | fgrep -w DUP \
     | fgrep -wf  gene_lists/gold_standard.triplosensitive.genes.list ) \
| cut -f4,5,12 | sort -nrk3,3 | head -n10


# Get constrained genes & PIPs to higlight in bottom of graphical abstract
cat gene_lists/gnomad.v2.1.1.*_constrained.genes.list | sort | uniq \
| fgrep -wf - <( zcat rCNV.final_genes.genes.bed.gz ) | cut -f4,5,12 \
| sort -nrk3,3 | head -n10


# Get median & IQR of GD effect sizes across phenotypes
# (Note: category 3 = GDs & 7 = constrained genes outside of GDs)
for categ in 3 7; do
  cat << EOF > get_OR_distrib.cat$categ.R
x <- read.table("${prefix}.global_burden_stats.tsv.gz", header=T, comment.char="")
round(summary(exp(x[which(x\$category==$categ & x\$CNV=="CNV"), "meta_lnOR"])), 2)
EOF
Rscript get_OR_distrib.cat$categ.R
rm get_OR_distrib.cat$categ.R
done

# Get interquartile range of CNV sizes
Rscript -e "summary(abs(apply(read.table('cnv/mega.rCNV.bed.gz', \
            header=T, sep=\"\\\t\", comment.char=\"\")[, 2:3], 1, diff)))"


######################################
# FREQUENCY COMPARISONS TO GNOMAD-SV #
######################################
# Note: everything below this header is related to frequency comparisons
# of CNVs per cohort vs gnomAD-SV v2.1

# Install svtk for clustering CNVs
conda install -y Cython && \
cd /opt && \
git clone https://github.com/broadinstitute/gatk-sv.git && \
cd gatk-sv/src/svtk && \
pip install -e . && \
cd /

# Cluster CNVs per cohort
while read meta cohorts; do
  echo $meta
  bedtools intersect -wa -wb -r -f 0.5 \
    -a cnv/$meta.rCNV.bed.gz -b cnv/$meta.rCNV.bed.gz \
  | awk -v OFS="\t" -v dist=100000 \
    '{ if (sqrt(($8-$2)^2) <= dist && sqrt(($9-$3)^2) <= dist) \
       print $1, $2, $3, $4, $4, $5, $7, $8, $9, $10, $10, $11 }' \
  | bgzip -c > cnv/$meta.preclustered.bed.gz
  tabix -f cnv/$meta.preclustered.bed.gz
  zcat cnv/$meta.preclustered.bed.gz | cut -f1-6 | sort -Vk1,1 -k2,2n -k3,3n | uniq \
  | bgzip -c > cnv/$meta.reformatted.bed.gz
  tabix -f cnv/$meta.reformatted.bed.gz
  cidx=$( head -n1 refs/HPOs_by_metacohort.table.tsv | sed 's/\t/\n/g' \
          | awk -v meta=$meta '{ if ($1==meta) print NR }' )
  n_ctrl=$( fgrep -w HEALTHY_CONTROL refs/HPOs_by_metacohort.table.tsv | cut -f$cidx )
  n_case=$( fgrep -w "HP:0000118" refs/HPOs_by_metacohort.table.tsv | cut -f$cidx )
  # Cluster one chromosome at a time for memory purposes
  for contig in $( seq 1 22 ); do
    tabix -h cnv/$meta.preclustered.bed.gz $contig > preclust.bed
    tabix -h cnv/$meta.reformatted.bed.gz $contig > cnvs.bed
    svtk bedcluster -f 0.5 -p "${meta}_clustered" -m -s preclust.bed cnvs.bed
  done \
  | awk -v n_case=$n_case -v n_ctrl=$n_ctrl -v OFS="\t" \
    '{ print $1, $2, $3, $4, $5, $9/(n_case+n_ctrl) }' \
  | sort -Vk1,1 -k2,2n -k3,3n -k5,5V | uniq | fgrep -v "#" \
  | cat <( echo -e "#chr\tstart\tend\tvid\tcnv\tfreq" ) - \
  | bgzip -c > cnv/$meta.clustered.bed.gz
  rm preclust.bed cnvs.bed
done < <( fgrep -v mega refs/rCNV_metacohort_list.txt )

# Filter gnomAD-SV in prior to comparison
wget https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.bed.gz
wget https://storage.googleapis.com/gcp-public-data--gnomad/papers/2019-sv/gnomad_v2.1_sv.sites.bed.gz.tbi
svtype_idx=$( tabix -H gnomad_v2.1_sv.sites.bed.gz | sed 's/\t/\n/g' | awk '{ if ($1=="SVTYPE") print NR }' )
fref_idx=$( tabix -H gnomad_v2.1_sv.sites.bed.gz | sed 's/\t/\n/g' | awk '{ if ($1=="FREQ_HOMREF") print NR }' )
cpx_idx=$( tabix -H gnomad_v2.1_sv.sites.bed.gz | sed 's/\t/\n/g' | awk '{ if ($1=="CPX_INTERVALS") print NR }' )
af_idx=$( tabix -H gnomad_v2.1_sv.sites.bed.gz | sed 's/\t/\n/g' | awk '{ if ($1=="AF") print NR }' )
# Biallelic CNVs
zcat gnomad_v2.1_sv.sites.bed.gz | fgrep -w PASS | grep -e '^[0-9]' \
| awk -v FS="\t" -v OFS="\t" -v svtype_idx=$svtype_idx -v fref_idx=$fref_idx \
  '{ if ($svtype_idx ~ /DEL|DUP/) print $1, $2, $3, $4, $svtype_idx, 1-$fref_idx }' \
> gnomAD_SV.cnvs.bed
# Complex SVs
while read vid ints freq; do
  echo $ints | sed -e 's/,/\n/g' -e 's/[_\:-]/\t/g' \
  | awk -v OFS="\t" -v vid=$vid -v freq="$freq" \
    '{ if ($1 ~ /DEL|DUP/) print $2, $3, $4, vid, $1, freq }'
done < <( zcat gnomad_v2.1_sv.sites.bed.gz | fgrep -w PASS | grep -e '^[0-9]' | fgrep -w CPX \
          | awk -v FS="\t" -v OFS="\t" -v cpx_idx=$cpx_idx -v fref_idx=$fref_idx \
            '{ print $4, $cpx_idx, 1-$fref_idx }' ) \
>> gnomAD_SV.cnvs.bed
# Multiallelic CNVs
while read chrom start end vid afs; do
  del_freq=$( echo $afs | sed 's/,/\n/g' | head -n2 | awk '{ sum+=$1 }END{ print sum }' )
  dup_freq=$( echo $afs | sed 's/,/\n/g' | sed '1,3d' | awk '{ sum+=$1 }END{ print sum }' )
  echo -e "$chrom\t$start\t$end\t$vid" \
  | awk -v OFS="\t" -v del_freq="$del_freq" -v dup_freq="$dup_freq" \
    '{ print $0, "DUP", dup_freq"\n"$0, "DEL", del_freq }' \
  | awk '{ if ($6>0) print $0 }'
done < <( zcat gnomad_v2.1_sv.sites.bed.gz | grep -e '^[0-9]' | fgrep -w MCNV \
          | awk -v FS="\t" -v OFS="\t" -v af_idx=$af_idx \
            '{ if ($NF=="MULTIALLELIC") print $1, $2, $3, $4, $af_idx }' ) \
>> gnomAD_SV.cnvs.bed
# Sort & compress
sort -Vk1,1 -k2,2n -k3,3n -k5,5V -k6,6n gnomAD_SV.cnvs.bed \
| cat <( echo -e "#chr\tstart\tend\tvid\tcnv\tfreq" ) - \
| bgzip -c > gnomAD_SV.cnvs.bed.gz

# Intersect clustered CNVs vs. gnomAD-SV
while read meta cohorts; do
  bedtools intersect -wao -r -f 0.5 \
    -a cnv/$meta.clustered.bed.gz \
    -b gnomAD_SV.cnvs.bed.gz \
  | awk -v OFS="\t" '{ if ($12==".") $12="NA"; if ($5==$11 || $11==".") print $1, $2, $3, $4, $5, $6, $12, $13 }' \
  | sort -Vk1,1 -k2,2n -k3,3n -k5,5V -k4,4V -k8,8nr \
  | cat <( echo -e "#chrom\tstart\tend\tvid\tcnv\tfreq\tgnomad\toverlap" ) - \
  | bgzip -c > $meta.cnv_vs_gnomad.bed.gz
done < <( fgrep -v mega refs/rCNV_metacohort_list.txt )

# Compare frequencies
if [ -e gnomad_frequency_comparisons ]; then
  rm -rf gnomad_frequency_comparisons
fi
mkdir gnomad_frequency_comparisons
while read meta cohorts; do
  /opt/rCNV2/analysis/paper/plot/misc/plot_gnomad_freq_comparisons.R \
    $meta.cnv_vs_gnomad.bed.gz \
    gnomad_frequency_comparisons/$meta
done < <( fgrep -v mega refs/rCNV_metacohort_list.txt )

# Copy plots to rCNV bucket (note: requires permissions)
gsutil -m cp -r \
  gnomad_frequency_comparisons \
  ${rCNV_bucket}/analysis/paper/plots/misc/
