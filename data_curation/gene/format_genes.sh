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


# Add helper alias for formatting long integers
alias addcom="sed -e :a -e 's/\(.*[0-9]\)\([0-9]\{3\}\)/\1,\2/;ta'"


# Download analysis references
mkdir refs/
gsutil -m cp gs://rcnv_project/analysis/analysis_refs/** \
  refs/


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


# Make gene list of all canonical autosomal genes used in analysis
zcat gencode.v19.canonical.gtf.gz \
| awk -v FS="\t" '{ if ($3=="gene") print $9 }' \
| sed 's/;/\n/g' | fgrep "gene_name" \
| sed 's/gene_name\ //g' | tr -d '"' \
| awk '{ print $1 }' | sort -Vk1,1 | uniq \
> gencode.v19.canonical.genes.list


# Copy canonical gene metadata to rCNV bucket (note: requires permissions)
gsutil cp gencode.v19.canonical.*.gz gs://rcnv_project/cleaned_data/genes/
gsutil cp gencode.v19.canonical.genes.list \
  gs://rcnv_project/cleaned_data/genes/gene_lists/


# Generate BED file of pLoF constrained genes from gnomAD
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
pli_idx=$( zcat gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
           | head -n1 | sed 's/\t/\n/g' \
           | awk '{ if ($1=="pLI") print NR }' )
loeuf_idx=$( zcat gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
             | head -n1 | sed 's/\t/\n/g' \
             | awk '{ if ($1=="oe_lof_upper_bin_6") print NR }' )
zcat gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
| fgrep -wf gencode.v19.canonical.genes.list \
| awk -v pli_idx=${pli_idx} -v loeuf_idx=${loeuf_idx} -v FS="\t" \
  '{ if ($pli_idx >= 0.9 || $loeuf_idx == 0 ) print $1 }' \
| sort -Vk1,1 \
> gnomad.v2.1.1.lof_constrained.genes.list
awk '{ print "gene_name \""$1"\";" }' \
gnomad.v2.1.1.lof_constrained.genes.list \
| fgrep -wf - <( zcat gencode.v19.canonical.gtf.gz ) \
| awk -v FS="\t" '{ if ($3=="transcript") print }' \
> gencode.v19.canonical.constrained.gtf
while read gene; do
  fgrep -w "\"$gene\";" gencode.v19.canonical.constrained.gtf \
  | awk -v OFS="\t" -v gene=$gene '{ print $1, $4, $5, gene }'
done < gnomad.v2.1.1.lof_constrained.genes.list \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
| cat <( echo -e "#chr\tstart\tend\tgene" ) - \
| bgzip -c \
> gencode.v19.canonical.constrained.bed.gz


# Copy constrained gene files to rCNV bucket (note: requires permissions)
gsutil cp gencode.v19.canonical.constrained.bed.gz \
  gs://rcnv_project/analysis/analysis_refs/
gsutil cp gnomad.v2.1.1.lof_constrained.genes.list \
  gs://rcnv_project/cleaned_data/genes/gene_lists/


# Generate gene list of all genes associated with each HPO term
wget http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt
while read pheno hpo; do
  awk -v hpo=${hpo} '{ if ($1==hpo) print $NF }' \
    ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt \
  | fgrep -wf gencode.v19.canonical.genes.list \
  | sort -Vk1,1 | uniq \
  > $pheno.HPOdb.genes.list
done < refs/test_phenotypes.list
# For HP0000118 & UNKNOWN, take union of all genes associated with any HPO term
sed '1d' ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt \
| awk '{ print $NF }' \
| sort -Vk1,1 | uniq \
> HP0000118.HPOdb.genes.list
cp HP0000118.HPOdb.genes.list UNKNOWN.HPOdb.genes.list


# Generate gene list of all constrained genes associated with each HPO term
while read pheno hpo; do
  fgrep -wf gnomad.v2.1.1.lof_constrained.genes.list \
    $pheno.HPOdb.genes.list \
  > $pheno.HPOdb.constrained.genes.list
done < refs/test_phenotypes.list


# Copy HPO-based gene lists to rCNV bucket (note: requires permissions)
gsutil -m cp *.HPOdb.*genes.list gs://rcnv_project/cleaned_data/genes/gene_lists/


# Generate HTML table of gene lists for README
# Header
echo "| Gene set | Genes | Filename prefix | Source | Description |" \
> genelist_table.html.txt
echo "| :--- | ---: | :--- | :--- | :---- | " \
>> genelist_table.html.txt
# All genes
for wrapper in 1; do
  echo "All genes"
  wc -l gencode.v19.canonical.genes.list \
  | awk -v OFS="\t" '{ print $1, "`"$2"`" }' \
  | sed 's/\.genes\.list//g' | addcom
  echo "Gencode v19 [Harrow _et al._, _Genome Res._, 2012](https://www.ncbi.nlm.nih.gov/pubmed/22955987)"
  echo "Canonical transcripts from autosomal, protein-coding genes"
done | paste -s \
| sed 's/\t/\ \|\ /g' \
| sed -e 's/^/\|\ /g' -e 's/$/\ \|/g' \
>> genelist_table.html.txt
# Constrained genes
for wrapper in 1; do
  echo "Constrained genes"
  wc -l gnomad.v2.1.1.lof_constrained.genes.list \
  | awk -v OFS="\t" '{ print $1, "`"$2"`" }' \
  | sed 's/\.genes\.list//g' | addcom
  echo "gnomAD v2.1.1 [Karczewski _et al._, _bioRxiv_, 2019](https://www.biorxiv.org/content/10.1101/531210v3)"
  echo "pLI ≥ 0.9 or in the first LOEUF sextile"
done | paste -s \
| sed 's/\t/\ \|\ /g' \
| sed -e 's/^/\|\ /g' -e 's/$/\ \|/g' \
>> genelist_table.html.txt
# HPO-associated genes
while read pheno hpo; do
  for wrapper in 1; do
    descrip=$( fgrep -w $hpo refs/HPOs_by_metacohort.table.tsv | cut -f2 )
    echo "$descrip (${hpo})-associated genes"
    cat $pheno.HPOdb.genes.list | wc -l | addcom
    echo "\`$pheno.HPOdb\`"
    echo "HPO database (accessed $( date | awk '{ print $2, $NF }' )) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478)"
    echo "Genes linked to $hpo"
  done | paste -s \
  | sed 's/\t/\ \|\ /g' \
  | sed -e 's/^/\|\ /g' -e 's/$/\ \|/g'
done < refs/test_phenotypes.list \
>> genelist_table.html.txt
cat genelist_table.html.txt



