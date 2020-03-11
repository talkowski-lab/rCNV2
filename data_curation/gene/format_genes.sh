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
gsutil -m cp gs://rcnv_project/analysis/analysis_refs/** refs/
gsutil -m cp gs://rcnv_project/refs/** refs/


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


# Apply pext filter to exons from canonical transcripts
# Note: pext pre-formatting takes several hours locally; was parallelized in cloud 
# using process_pext.wdl, and can be downloaded (with permissions) at 
# gs://rcnv_project/cleaned_data/genes/annotations/gnomad.v2.1.1.pext.bed.gz
gsutil -m cp \
  gs://gnomad-public/papers/2019-tx-annotation/pre_computed/all.possible.snvs.tx_annotated.022719.tsv.bgz \
  ./
tabix -s 1 -b 2 -e 2 -S 1 \
  all.possible.snvs.tx_annotated.022719.tsv.bgz
/opt/rCNV2/data_curation/gene/process_pext.py \
  --pan-tissue \
  all.possible.snvs.tx_annotated.022719.tsv.bgz \
| sort -Vk1,1 -k2,2n -k3,3n \
| bgzip -c \
> gnomad.v2.1.1.pext.bed.gz
tabix -f gnomad.v2.1.1.pext.bed.gz
/opt/rCNV2/data_curation/gene/apply_pext_filter.py \
  --min-pext 0.2 \
  -o gencode.v19.canonical.pext_filtered.gtf.gz \
  --bgzip \
  --lost-genes genes_lost_during_pext_filtering.genes.list \
  gencode.v19.canonical.gtf.gz \
  gnomad.v2.1.1.pext.bed.gz


# Make gene list of all canonical autosomal genes used in analysis
for gtf in gencode.v19.canonical gencode.v19.canonical.pext_filtered; do
  zcat ${gtf}.gtf.gz \
  | awk -v FS="\t" '{ if ($3=="transcript") print $9 }' \
  | sed 's/;/\n/g' | fgrep "gene_name" \
  | sed 's/gene_name\ //g' | tr -d '"' \
  | awk '{ print $1 }' | sort -Vk1,1 | uniq \
  > ${gtf}.genes.list
done
zcat gencode.v19.canonical.tsv.gz \
| fgrep -wf gencode.v19.canonical.pext_filtered.genes.list \
| gzip -c \
> gencode.v19.canonical.pext_filtered.tsv.gz


# Preprocess GTEx expression matrix to compute summaries
wget https://storage.googleapis.com/gtex_analysis_v7/rna_seq_data/GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz
wget https://storage.googleapis.com/gtex_analysis_v7/annotations/GTEx_v7_Annotations_SampleAttributesDS.txt
mkdir gtex_stats/
/opt/rCNV2/data_curation/gene/preprocess_GTEx.py \
  --gzip \
  --prefix gtex_stats/gencode.v19.canonical.pext_filtered.GTEx_v7_expression_stats \
  gencode.v19.canonical.pext_filtered.tsv.gz \
  GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_tpm.gct.gz \
  GTEx_v7_Annotations_SampleAttributesDS.txt
# Copy precomputed GTEx summary data to rCNV bucket (note: requires permissions)
gsutil -m cp -r gtex_stats \
  gs://rcnv_project/cleaned_data/genes/annotations/


# Gather per-gene metadata (genomic)
# Note: this is parallelized in FireCloud with get_gene_metadata.wdl
wget ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
bgzip Homo_sapiens.GRCh37.dna.primary_assembly.fa
samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
for wrapper in 1; do 
  echo -e "refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz\tcoverage\trepeat_cov"
  echo -e "refs/GRCh37.somatic_hypermutable_sites.bed.gz\tcoverage\thypermutable_cov"
  echo -e "https://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/wgEncodeCrgMapabilityAlign100mer.bw\tmap-mean\talign_100mer"
done > gene_features.athena_tracklist.tsv
/opt/rCNV2/data_curation/gene/get_gene_features.py \
  --get-genomic \
  --centro-telo-bed refs/GRCh37.centromeres_telomeres.bed.gz \
  --ref-fasta Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz \
  --athena-tracks gene_features.athena_tracklist.tsv \
  --outbed gencode.v19.canonical.pext_filtered.genomic_features.bed.gz \
  --bgzip \
  gencode.v19.canonical.pext_filtered.gtf.gz


# Gather per-gene metadata (expression)
/opt/rCNV2/data_curation/gene/get_gene_features.py \
  --get-expression \
  --gtex-medians gtex_stats/gencode.v19.canonical.pext_filtered.GTEx_v7_expression_stats.median.tsv.gz \
  --gtex-mads gtex_stats/gencode.v19.canonical.pext_filtered.GTEx_v7_expression_stats.mad.tsv.gz \
  --outbed gencode.v19.canonical.pext_filtered.expression_features.bed.gz \
  --bgzip \
  gencode.v19.canonical.pext_filtered.gtf.gz


# Gather per-gene constraint metadata
wget https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
wget http://genic-intolerance.org/data/RVIS_Unpublished_ExACv2_March2017.txt
gsutil -m cp gs://rcnv_project/cleaned_data/genes/annotations/EDS.Wang_2018.tsv.gz ./
# Add: promoter conservation, average exon conservation, EDS
/opt/rCNV2/data_curation/gene/get_gene_features.py \
  --get-constraint \
  --ref-fasta Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz \
  --gnomad-constraint gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
  --rvis-tsv RVIS_Unpublished_ExACv2_March2017.txt \
  --eds-tsv EDS.Wang_2018.tsv.gz \
  --outbed gencode.v19.canonical.pext_filtered.constraint_features.bed.gz \
  --bgzip \
  gencode.v19.canonical.pext_filtered.gtf.gz


# Join all per-gene metadata
/opt/rCNV2/data_curation/gene/join_gene_metadata.R \
  gencode.v19.canonical.pext_filtered.genomic_features.bed.gz \
  gencode.v19.canonical.pext_filtered.expression_features.bed.gz \
  gencode.v19.canonical.pext_filtered.constraint_features.bed.gz \
| bgzip -c \
> gencode.v19.canonical.pext_filtered.all_features.bed.gz


# Copy canonical gene metadata to rCNV bucket (note: requires permissions)
gsutil -m cp gencode.v19.canonical*.gz gs://rcnv_project/cleaned_data/genes/
gsutil -m cp gencode.v19.canonical*genes.list \
  gs://rcnv_project/cleaned_data/genes/gene_lists/
gsutil -m cp genes_lost_during_pext_filtering.genes.list \
  gs://rcnv_project/cleaned_data/genes/gene_lists/


# Generate BED file of pLoF constrained genes from gnomAD
# wget https://storage.googleapis.com/gnomad-public/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz
pli_idx=$( zcat gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
           | head -n1 | sed 's/\t/\n/g' \
           | awk '{ if ($1=="pLI") print NR }' )
loeuf_idx=$( zcat gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
             | head -n1 | sed 's/\t/\n/g' \
             | awk '{ if ($1=="oe_lof_upper_bin_6") print NR }' )
zcat gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz \
| fgrep -wf gencode.v19.canonical.pext_filtered.genes.list \
| awk -v pli_idx=${pli_idx} -v loeuf_idx=${loeuf_idx} -v FS="\t" \
  '{ if ($pli_idx >= 0.9 || $loeuf_idx == 0 ) print $1 }' \
| sort -Vk1,1 \
> gnomad.v2.1.1.lof_constrained.genes.list
awk '{ print "gene_name \""$1"\"" }' \
gnomad.v2.1.1.lof_constrained.genes.list \
| fgrep -wf - <( zcat gencode.v19.canonical.pext_filtered.gtf.gz ) \
| awk -v FS="\t" '{ if ($3=="transcript") print }' \
> gencode.v19.canonical.pext_filtered.constrained.gtf
while read gene; do
  fgrep -w "\"$gene\"" gencode.v19.canonical.pext_filtered.constrained.gtf \
  | awk -v OFS="\t" -v gene=$gene '{ print $1, $4, $5, gene }'
done < gnomad.v2.1.1.lof_constrained.genes.list \
| sort -Vk1,1 -k2,2n -k3,3n -k4,4V \
| cat <( echo -e "#chr\tstart\tend\tgene" ) - \
| bgzip -c \
> gencode.v19.canonical.pext_filtered.constrained.bed.gz


# Copy constrained gene files to rCNV bucket (note: requires permissions)
gsutil -m cp gencode.v19.canonical.pext_filtered.constrained.bed.gz \
  gs://rcnv_project/analysis/analysis_refs/
gsutil -m cp gnomad.v2.1.1.lof_constrained.genes.list \
  gs://rcnv_project/cleaned_data/genes/gene_lists/


# Generate gene list of all genes associated with each HPO term
wget http://compbio.charite.de/jenkins/job/hpo.annotations.monthly/lastSuccessfulBuild/artifact/annotation/ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt
while read pheno hpo; do
  awk -v hpo=${hpo} '{ if ($1==hpo) print $NF }' \
    ALL_SOURCES_ALL_FREQUENCIES_phenotype_to_genes.txt \
  | fgrep -wf gencode.v19.canonical.pext_filtered.genes.list \
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


# Generate ClinGen dosage sensitive gene lists
wget ftp://ftp.clinicalgenome.org/ClinGen_gene_curation_list_GRCh37.tsv



# Generate HTML table of gene lists for README
# Header
echo "| Gene set | Genes | Filename prefix | Source | Description |" \
> genelist_table.html.txt
echo "| :--- | ---: | :--- | :--- | :---- | " \
>> genelist_table.html.txt
# All genes
for wrapper in 1; do
  echo "All genes"
  wc -l gencode.v19.canonical.pext_filtered.genes.list \
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

