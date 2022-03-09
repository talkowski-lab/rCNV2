#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Curate all genome annotation tracks used for noncoding association testing


# Launch docker image
docker run --rm -it gcr.io/gnomad-wgs-v2-sv/rcnv
gcloud auth login


# Set global paramters
export rCNV_bucket="gs://rcnv_project"


# Download necessary reference files
mkdir refs/
gsutil -m cp \
  ${rCNV_bucket}/refs/** \
  ${rCNV_bucket}/analysis/analysis_refs/* \
  ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists \
  refs/
mkdir cnvs/
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/cnv/*bed.gz* \
  ${rCNV_bucket}/cleaned_data/cnv/noncoding/** \
  cnvs/


# Prepare clustered somatic hypermutability reference blacklist
bedtools merge -d 200000 -i refs/GRCh37.somatic_hypermutable_sites.bed.gz \
| bgzip -c \
> GRCh37.somatic_hypermutable_sites.200kb_clustered.bed.gz
gsutil -m cp \
  GRCh37.somatic_hypermutable_sites.200kb_clustered.bed.gz \
  ${rCNV_bucket}/refs/


# Preprocess Roadmap Epigenomics ChromHMM states
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/all.mnemonics.bedFiles.tgz
mkdir roadmap_raw/
tar -xzvf all.mnemonics.bedFiles.tgz -C roadmap_raw/
wget https://egg2.wustl.edu/roadmap/data/byFileType/chromhmmSegmentations/ChmmModels/core_K27ac/jointModel/final/labelmap_18_core_K27ac.tab
gsutil -m cp ${rCNV_bucket}/refs/REP_sample_manifest.tsv ./
mkdir chromhmm_beds/
/opt/rCNV2/data_curation/genome_annotations/preprocess_chromhmm.py \
  --state-manifest labelmap_18_core_K27ac.tab \
  --sample-manifest REP_sample_manifest.tsv \
  --prefix chromhmm_beds/roadmap_chromhmm \
  roadmap_raw/
# Copy ChromHMM beds (and tracklist) to gs:// bucket (note: requires permissions)
gsutil -m cp -r \
  chromhmm_beds \
  ${rCNV_bucket}/cleaned_data/genome_annotations/
find chromhmm_beds/ -name "*.bed.gz" \
| awk -v gs=${rCNV_bucket} '{ print gs"/cleaned_data/genome_annotations/"$1 }' \
> chromhmm_tracks.gs_paths.list
gsutil -m cp \
  chromhmm_tracks.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Preprocess ENCODE TAD boundaries
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/genome_annotations/manifests/encode.hic_tads.manifest.tsv.gz \
  ./
mkdir encode_tad_boundaries/
while read trackname path; do
  wget -O $trackname.raw.bed.gz $path
  bedtools flank \
    -i $trackname.raw.bed.gz \
    -g <( awk '{ print "chr"$0 }' refs/GRCh37.genome ) \
    -b 5000 \
  | cut -f1-3 \
  | sort -Vk1,1 -k2,2n -k3,3n \
  | bgzip -c \
  > encode_tad_boundaries/$trackname.tad_boundaries.bed.gz
done < <( zcat encode.hic_tads.manifest.tsv.gz \
          | awk -v FS="\t" -v OFS="\t" '{ print $1, $NF }' \
          | sed '1d' )
# Copy ENCODE beds (and tracklist) to gs:// bucket (note: requires permissions)
gsutil -m cp -r \
  encode_tad_boundaries \
  ${rCNV_bucket}/cleaned_data/genome_annotations/
find encode_tad_boundaries/ -name "*.bed.gz" \
| awk -v gs=${rCNV_bucket} -v OFS="\t" \
  '{ print gs"/cleaned_data/genome_annotations/"$1, "encode_tad_boundaries" }' \
> encode.tad_boundaries.gs_paths.list
gsutil -m cp \
  encode.tad_boundaries.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Curate HACER enhancers
mkdir hacer_enhancers
gsutil -m cp \
  ${rCNV_bucket}/raw_data/genome_annotations/HACER.GROseq.manifest.tsv \
  ${rCNV_bucket}/raw_data/genome_annotations/HACER.PROseq.manifest.tsv \
  ${rCNV_bucket}/raw_data/genome_annotations/HACER.CAGE.manifest.tsv \
  ./
while read sample url; do
  echo $sample
  if [ $( echo $url | grep -e '\.zip$' | wc -l ) -gt 0 ]; then
    wget $url
    if [ -e $sample/ ]; then
      rm -r $sample/
    fi
    mkdir $sample/
    unzip $( basename $url )  -d $sample
    for file in $sample/*; do
      tname=$( basename $file | sed 's/.txt//g' | cut -f2- -d\- | paste -s -d \- )
      sed '1d' $file | cut -f2-4 | bgzip -c > hacer_enhancers/$tname.bed.gz
    done
  else
    wget $url -q -O - | sed '1d' | cut -f2-4 | bgzip -c > hacer_enhancers/$sample.bed.gz
  fi
done < <( cat HACER.GROseq.manifest.tsv HACER.PROseq.manifest.tsv )
gsutil -m cp -r \
  hacer_enhancers \
  ${rCNV_bucket}/cleaned_data/genome_annotations/
find hacer_enhancers/ -name "*.bed.gz" \
| awk -v gs=${rCNV_bucket} '{ print gs"/cleaned_data/genome_annotations/"$1 }' \
| cat - <( awk -v FS="\t" '{ print $2 }' HACER.CAGE.manifest.tsv ) \
| awk -v OFS="\t" '{ print $1, "hacer_enh" }' \
> HACER.track_urls.list
gsutil -m cp \
  HACER.track_urls.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Curate super enhancers and typical enhancers per sample from SEdb
gsutil cp ${rCNV_bucket}/raw_data/genome_annotations/SEdb.sample_manifest.tsv ./
mkdir SEdb_super_enhancers/
/opt/rCNV2/data_curation/genome_annotations/preprocess_SEdb.py \
  --format 'SE' \
  --outdir SEdb_super_enhancers/ \
  http://www.licpathway.net/sedb/download/package/SE_package.bed \
  SEdb.sample_manifest.tsv
mkdir SEdb_enhancers
/opt/rCNV2/data_curation/genome_annotations/preprocess_SEdb.py \
  --format 'TE' \
  --outdir SEdb_enhancers/ \
  http://www.licpathway.net/sedb/download/package/TE_package.bed \
  SEdb.sample_manifest.tsv
gsutil -m cp -r \
  SEdb_super_enhancers \
  SEdb_enhancers \
  ${rCNV_bucket}/cleaned_data/genome_annotations/
find SEdb_super_enhancers/ -name "*bed.gz" \
|  awk -v gs=${rCNV_bucket} '{ print gs"/cleaned_data/genome_annotations/"$1 }' \
> SEdb.track_gspaths.list
find SEdb_enhancers/ -name "*bed.gz" \
|  awk -v gs=${rCNV_bucket} '{ print gs"/cleaned_data/genome_annotations/"$1 }' \
>> SEdb.track_gspaths.list
gsutil -m cp \
  SEdb.track_gspaths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Curate super enhancers from dbSUPER
wget http://asntech.org/dbsuper/data/bed/hg19/all_hg19_bed.zip
mkdir dbSUPER_super_enhancers/
unzip all_hg19_bed.zip
for file in all_hg19_bed/*.bed; do
  newname="$( basename "$file" | sed -e 's/\.bed//g' -e 's/\ /_/g' ).bed"
  mv "$file" dbSUPER_super_enhancers/$newname
  bgzip -f dbSUPER_super_enhancers/$newname
done
gsutil -m cp -r \
  dbSUPER_super_enhancers \
  ${rCNV_bucket}/cleaned_data/genome_annotations/
find dbSUPER_super_enhancers/ -name "*.bed.gz" \
| awk -v gs=${rCNV_bucket} -v OFS="\t" \
  '{ print gs"/cleaned_data/genome_annotations/"$1, "dbSUPER" }' \
> dbSUPER.gs_paths.list
gsutil -m cp \
  dbSUPER.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Curate VISTA validated enhancers
gsutil -m cp ${rCNV_bucket}/raw_data/genome_annotations/vista_raw.fa.gz ./
zcat vista_raw.fa.gz \
| fgrep positive | fgrep Human \
| awk -v FS="|" '{ print $2 }' \
| sed -e 's/:\|-/\t/g' \
| awk -v OFS="\t" '{ print $1, $2, $3 }' \
| sort -Vk1,1 -k2,2n -k3,3n \
| bgzip -c \
> vista_enhancers.bed.gz
gsutil -m cp vista_enhancers.bed.gz \
  ${rCNV_bucket}/cleaned_data/genome_annotations/misc_tracks/


# Curate DENdb enhancers split by cell type
wget https://www.cbrc.kaust.edu.sa/dendb/src/enhancers.csv.zip
unzip enhancers.csv.zip
mkdir DENdb_enhancers
/opt/rCNV2/data_curation/genome_annotations/preprocess_DENdb.py \
  enhancers.csv \
  DENdb_enhancers/
gsutil -m cp -r \
  DENdb_enhancers \
  ${rCNV_bucket}/cleaned_data/genome_annotations/  
find DENdb_enhancers/ -name "*.bed.gz" \
| awk -v gs=${rCNV_bucket} -v OFS="\t" \
  '{ print gs"/cleaned_data/genome_annotations/"$1 }' \
> DENdb.gs_paths.list
gsutil -m cp \
  DENdb.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Curate SEA super enhancers split by cell type
wget http://113.6.251.5:8080/SEA3/download/SEA00101.bed
awk -v FS="\t" -v OFS="\t" '{ print $2, $3, $4, $7 }' SEA00101.bed \
| sed 's/\ /_/g' | bgzip -c \
> SEA.hg38.gz
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
liftOver -minMatch=0.75 -bedPlus=4 \
  SEA.hg38.gz \
  hg38ToHg19.over.chain.gz \
  SEA.hg19.bed \
  SEA.hg38_liftFailed.txt
mkdir SEA_super_enhancers/
/opt/rCNV2/data_curation/genome_annotations/preprocess_SEA.py \
  SEA.hg19.bed \
  SEA_super_enhancers/
gsutil -m cp -r \
  SEA_super_enhancers \
  ${rCNV_bucket}/cleaned_data/genome_annotations/  
find SEA_super_enhancers/ -name "*.bed.gz" \
| awk -v gs=${rCNV_bucket} -v OFS="\t" \
  '{ print gs"/cleaned_data/genome_annotations/"$1 }' \
> SEA.gs_paths.list
gsutil -m cp \
  SEA.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Curate various PsychENCODE datasets
wget http://resource.psychencode.org/Datasets/Derived/DER-03a_hg19_PEC_enhancers.bed
sed '1d' DER-03a_hg19_PEC_enhancers.bed \
| awk -v OFS="\t" '{ print $1, $2, $3 }' \
| sort -Vk1,1 -k2,2n -k3,3n \
| bgzip -c \
> psychencode.enhancers.bed.gz
curl -s http://resource.psychencode.org/Datasets/Derived/DER-18_TAD_adultbrain.bed \
| awk -v OFS="\t" -v buffer=5000 \
  '{ print $1, $2-buffer, $2+buffer"\n"$1, $3-buffer, $3+buffer }' \
  DER-18_TAD_adultbrain.bed \
| awk -v OFS="\t" '{ if ($2<0) $2=0; print $1, $2, $3 }' \
| sort -Vk1,1 -k2,2n -k3,3n | uniq \
| bgzip -c \
> psychencode.tad_boundaries.bed.gz
gsutil -m cp \
  psychencode.enhancers.bed.gz \
  psychencode.tad_boundaries.bed.gz \
  ${rCNV_bucket}/cleaned_data/genome_annotations/misc_tracks/


# Curate human TFs from JASPAR 2020
gsutil -m cp \
  ${rCNV_bucket}/raw_data/genome_annotations/jaspar_2020.human.tf_ids.manifest.tsv.gz \
  ./
if ! [ -e hg38ToHg19.over.chain.gz ]; then
  wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz
fi
mkdir jaspar_tfbs
while read id tfname; do
  echo $tfname
  curl -s http://jaspar.genereg.net/download/bed_files/${id}.bed \
  | fgrep hg19 | cut -f1-3 \
  | sort -Vk1,1 -k2,2 -k3,3n \
  > jaspar_tfbs/${tfname}.bed
  if [ $( cat jaspar_tfbs/${tfname}.bed | wc -l ) -gt 0 ]; then
    bgzip -f jaspar_tfbs/${tfname}.bed
  else
    rm jaspar_tfbs/${tfname}.bed
    curl -s http://jaspar.genereg.net/download/bed_files/${id}.bed \
    | fgrep hg38 | cut -f1-3 \
    | sort -Vk1,1 -k2,2 -k3,3n \
    > jaspar_tfbs/${tfname}.hg38.bed
    if [ $( cat jaspar_tfbs/${tfname}.hg38.bed | wc -l ) -gt 0 ]; then
      liftOver \
        -minMatch=1.0 \
        -bedPlus=5 \
        jaspar_tfbs/${tfname}.hg38.bed \
        hg38ToHg19.over.chain.gz \
        jaspar_tfbs/${tfname}.bed \
        jaspar_tfbs/${tfname}.hg38_fail.bed
        rm jaspar_tfbs/${tfname}.hg38.bed jaspar_tfbs/${tfname}.hg38_fail.bed
        bgzip -f jaspar_tfbs/${tfname}.bed
    else
      rm jaspar_tfbs/${tfname}.hg38.bed
    fi
  fi
done < <( zcat jaspar_2020.human.tf_ids.manifest.tsv.gz | fgrep -v "#" )
gsutil -m cp -r \
  jaspar_tfbs \
  ${rCNV_bucket}/cleaned_data/genome_annotations/  
find jaspar_tfbs/ -name "*.bed.gz" \
| awk -v gs=${rCNV_bucket} -v OFS="\t" \
  '{ print gs"/cleaned_data/genome_annotations/"$1, "jaspar" }' \
> jaspar_tfbs.gs_paths.list
gsutil -m cp \
  jaspar_tfbs.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Curate BOCA brain chromatin accessibility data
wget https://bendlj01.u.hpc.mssm.edu//multireg/resources/boca_peaks.zip
mkdir boca_beds_tmp boca_beds
unzip boca_peaks.zip -d boca_beds_tmp/
wget \
  https://bendlj01.u.hpc.mssm.edu/multireg/resources/boca_peaks_consensus_no_blacklisted_regions.bed \
  -O boca_beds_tmp/BOCA_consensus.bed
for file in boca_beds_tmp/*bed; do
  awk -v FS="\t" '{ printf "%s\t%.0f\t%.0f\n", $1, $2, $3 }' $file \
  | bgzip -c \
  > boca_beds/$( basename $file ).gz
done
gsutil -m cp -r \
  boca_beds \
  ${rCNV_bucket}/cleaned_data/genome_annotations/
find boca_beds/ -name "*.bed.gz" \
| awk -v gs=${rCNV_bucket} -v OFS="\t" \
  '{ print gs"/cleaned_data/genome_annotations/"$1, "boca" }' \
> boca.gs_paths.list
gsutil -m cp \
  boca.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Curate misc. TAD datasets
mkdir misc_tad_boundaries
gsutil -m cp -r ${rCNV_bucket}/raw_data/genome_annotations/dixon_tbrs ./
if ! [ -e hg18ToHg19.over.chain.gz ]; then
  wget http://hgdownload.cse.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
fi
for ctype in hESC IMR90; do
  liftOver \
    -minMatch=0.5 \
    dixon_tbrs/dixon_2012.tad_boundaries.$ctype.hg18.bed.gz \
    hg18ToHg19.over.chain.gz \
    misc_tad_boundaries/dixon_2012.$ctype.bed \
    dixon_2012.liftFail.$ctype.bed
  bgzip -f misc_tad_boundaries/dixon_2012.$ctype.bed
done
for ctype in GM12878_primary+replicate K562 IMR90 HMEC KBM7 NHEK HUVEC HeLa; do
  echo $ctype
  curl -s https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_${ctype}_Arrowhead_domainlist.txt.gz \
  | gunzip -c | sed '1d' \
  | awk -v OFS="\t" -v buffer=5000 '{ print $1, $2-buffer, $2+buffer"\n"$1, $3-buffer, $3+buffer }' \
  | awk -v OFS="\t" '{ if ($2<0) $2=0; print $1, $2, $3 }' \
  | sort -Vk1,1 -k2,2n -k3,3n | uniq \
  | bgzip -c \
  > misc_tad_boundaries/rao_2014.$ctype.arrowhead_domain_boundaries.bed.gz
  curl -s https://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_${ctype}_HiCCUPS_looplist.txt.gz \
  | gunzip -c | sed '1d' \
  | awk -v OFS="\t" '{ print $1, $2, $3"\n"$4, $5, $6 }' \
  | awk -v OFS="\t" '{ if ($2<0) $2=0; print $1, $2, $3 }' \
  | sort -Vk1,1 -k2,2n -k3,3n | uniq \
  | bgzip -c \
  > misc_tad_boundaries/rao_2014.$ctype.loop_domain_boundaries.bed.gz
done
for ctype in CP GZ; do
  echo $ctype
  curl -s https://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77565/suppl/GSE77565_${ctype}_TAD.bed.gz \
  | gunzip -c | sed '1d' \
  | awk -v OFS="\t" -v buffer=5000 '{ print $1, $2-buffer, $2+buffer"\n"$1, $3-buffer, $3+buffer }' \
  | awk -v OFS="\t" '{ if ($2<0) $2=0; print $1, $2, $3 }' \
  | sort -Vk1,1 -k2,2n -k3,3n | uniq \
  | bgzip -c \
  > misc_tad_boundaries/gao_2016.$ctype.tad_boundaries.bed.gz
done
gsutil -m cp -r \
  misc_tad_boundaries \
  ${rCNV_bucket}/cleaned_data/genome_annotations/
find misc_tad_boundaries/ -name "*.bed.gz" \
| awk -v gs=${rCNV_bucket} -v OFS="\t" \
  '{ print gs"/cleaned_data/genome_annotations/"$1 }' \
> misc_tad_boundaries.gs_paths.list
gsutil -m cp \
  misc_tad_boundaries.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Curate noncoding RNAs from Gencode GTF
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.pc_transcripts.fa.gz
/opt/rCNV2/data_curation/gene/get_noncoding_rnas.py \
  --outdir gencode_noncoding_rnas/ \
  --prefix gencode.v19 \
  gencode.v19.annotation.gtf.gz \
  gencode.v19.pc_transcripts.fa.gz
gsutil -m cp -r \
  gencode_noncoding_rnas \
  ${rCNV_bucket}/cleaned_data/genome_annotations/
find gencode_noncoding_rnas/ -name "*.bed.gz" \
| awk -v gs=${rCNV_bucket} -v OFS="\t" \
  '{ print gs"/cleaned_data/genome_annotations/"$1 }' \
> gencode_noncoding_rnas.gs_paths.list
gsutil -m cp \
  gencode_noncoding_rnas.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Combine all enhancer tracks into single tracklist for sharding on FireCloud
gsutil cat \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/EnhancerAtlas2.track_urls.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/HACER.track_urls.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/SEdb.track_gspaths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/dbSUPER.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/DENdb.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/SEA.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/fantom_enh.track_urls.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/plac_brain_enhancer.gs_paths.list \
> enhancer_databases.track_urls.list
echo -e \
"${rCNV_bucket}/cleaned_data/genome_annotations/misc_tracks/vista_enhancers.bed.gz
${rCNV_bucket}/cleaned_data/genome_annotations/misc_tracks/psychencode.enhancers.bed.gz
http://resource.psychencode.org/Datasets/Derived/DER-03b_hg19_high_confidence_PEC_enhancers.bed\tpsychencode" \
>> enhancer_databases.track_urls.list
gsutil -m cp \
  enhancer_databases.track_urls.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/


# Combine all miscellaneous tracks into single tracklist for sharding on FireCloud
gsutil cat \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/encode.tad_boundaries.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/jaspar_tfbs.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/boca.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/misc_tad_boundaries.gs_paths.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/gencode_noncoding_rnas.gs_paths.list \
> misc_genome_annotations.track_urls.list
echo -e \
"http://resource.psychencode.org/Datasets/Derived/DER-05_PFC_H3K27ac_peaks.bed\tpsychencode
http://resource.psychencode.org/Datasets/Derived/DER-06_TC_H3K27ac_peaks.bed\tpsychencode
http://resource.psychencode.org/Datasets/Derived/DER-07_CBC_H3K27ac_peaks.bed\tpsychencode
http://resource.psychencode.org/Datasets/Pipeline/TARs/PIP-04_all_TARs.70pc.active.bed\tpsychencode
https://ccg.epfl.ch/UCNEbase/data/download/ucnes/hg19_UCNE_coord.bed\tUCNEbase
${rCNV_bucket}/cleaned_data/genome_annotations/dev_cerebrum_beds/dev_cerebrum_pREs.bed.gz
${rCNV_bucket}/cleaned_data/genome_annotations/dev_cerebrum_beds/dev_cerebrum_OCRs.bed.gz
${rCNV_bucket}/cleaned_data/genome_annotations/misc_tracks/psychencode.tad_boundaries.bed.gz" \
>> misc_genome_annotations.track_urls.list
gsutil -m cp \
  misc_genome_annotations.track_urls.list \
  ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/



# Count number of CNVs per cohort by pheno based on genic overlap
fgrep -wvf \
  refs/gene_lists/HP0000118.HPOdb.genes.list \
  refs/gene_lists/gnomad.v2.1.1.likely_unconstrained.genes.list \
| sort -V | uniq \
> loose_noncoding_whitelist.genes.list
mkdir genes_per_cnv
for CNV in DEL DUP; do
  echo "Starting $CNV"
  while read meta cohorts; do
    # Get conditional exclusion genes for that cohort as BED
    zcat refs/gencode.v19.canonical.pext_filtered.cohort_exclusion.bed.gz \
    | fgrep -w $meta | cut -f1-3 | bedtools merge -i - \
    | bgzip -c > $meta.cond_excl_genes.bed.gz

    # Annotate all CNVs based on any exon overlap
    /opt/rCNV2/analysis/genes/count_cnvs_per_gene.py \
      cnvs/$meta.rCNV.bed.gz \
      refs/gencode.v19.canonical.pext_filtered.gtf.gz \
      --min-cds-ovr "10e-10" \
      -t $CNV \
      --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
      --blacklist refs/GRCh37.somatic_hypermutable_sites.bed.gz \
      --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
      --blacklist $meta.cond_excl_genes.bed.gz \
      -o .junk.bed.gz \
      --bgzip \
      --cnvs-out /dev/stdout \
    | cut -f4,6-7,9 | gzip -c \
    > genes_per_cnv/rCNV.$CNV.$meta.genes_per_cnv.tsv.gz
    echo -e "$meta\tgenes_per_cnv/rCNV.$CNV.$meta.genes_per_cnv.tsv.gz"
  done < <( fgrep -v mega refs/rCNV_metacohort_list.txt ) \
  > genes_per_cnv.$CNV.input.tsv

  # Summarize CNV counts by gene list per phenotype
  /opt/rCNV2/analysis/paper/scripts/noncoding_association/summarize_ncCNV_counts.py \
    --genes-per-cnv genes_per_cnv.$CNV.input.tsv \
    --hpos <( cut -f2 refs/test_phenotypes.list ) \
    --unconstrained-genes loose_noncoding_whitelist.genes.list \
    --summary-counts unconstrained_cnv_counts.$CNV.tsv.gz \
    --developmental-hpos refs/rCNV2.hpos_by_severity.developmental.list \
    --gzip
  echo "Finished $CNV"
done
# Copy counts to Google bucket for future use
gsutil -m cp \
  unconstrained_cnv_counts.*.tsv.gz \
  ${rCNV_bucket}/analysis/crb_burden/other_data/





# Development parameters for curate_annotations.wdl
prefix="all_tracks"
tracklist="test.annotations.list"
min_element_size=5
max_element_size=200000
min_element_overlap=1.0
p_cutoff=0.05
track_prefix="test_annotations"

# Local dev only: localize CNV counts by gene context
# gsutil -m cp \
#   ${rCNV_bucket}/analysis/crb_burden/other_data/unconstrained_cnv_counts.*.tsv.gz \
#   ./

# Combine CNV counts by gene context into a single file
find / -name "unconstrained_cnv_counts.*.tsv.gz" | xargs -I {} mv {} ./ && \
cat \
  <( zcat unconstrained_cnv_counts.DEL.tsv.gz | head -n1 \
     | awk -v OFS="\t" '{ print $0, "cnv" }' ) \
  <( zcat unconstrained_cnv_counts.DEL.tsv.gz | grep -ve '^#' \
     | awk -v OFS="\t" '{ print $0, "DEL" }' ) \
  <( zcat unconstrained_cnv_counts.DUP.tsv.gz | grep -ve '^#' \
     | awk -v OFS="\t" '{ print $0, "DUP" }' ) \
| gzip -c > unconstrained_cnv_counts.tsv.gz

# Make dummy file of 10 annotations for development purposes
gsutil -m cat ${rCNV_bucket}/cleaned_data/genome_annotations/tracklists/misc_genome_annotations.track_urls.list \
| shuf | head -n10 > test.annotations.list

# Curate all annotations in an arbitrary input list of paths
while IFS=$'\t' read path tprefix; do
  if [ -z $tprefix ]; then
    tprefix=${track_prefix}
  fi
  echo -e "Curating $path"
  /opt/rCNV2/data_curation/genome_annotations/curate_track.py \
    --genome refs/GRCh37.genome \
    --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
    --blacklist refs/GRCh37.somatic_hypermutable_sites.200kb_clustered.bed.gz \
    --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
    --min-size ${min_element_size} \
    --max-size ${max_element_size} \
    --stats \
    --prefix "$tprefix" \
    "$path"
done < ${tracklist}
export IFS=$' \t\n'

# Compute CNV counts per annotation track per metacohort per CNV type
cat *.curated.stats.tsv | fgrep -v "#" | sort -Vk1,1 | uniq \
| awk -v FS="\t" -v OFS="\t" '{ print $0, $1".curated.bed.gz " }' \
| cat <( cat *.curated.stats.tsv | grep -e '^#' | sed -n '1p' \
         | awk -v OFS="\t" '{ print $0, "local_path" }' ) \
      - \
> ${prefix}.stats.tsv
while read cohort; do
  for dummy in 1; do
    echo $cohort
    fgrep -w $cohort refs/rCNV2.hpos_by_severity.developmental.counts.tsv | cut -f2
    cidx=$( sed -n '1p' refs/HPOs_by_metacohort.table.tsv \
            | sed 's/\t/\n/g' \
            | awk -v cohort=$cohort '{ if ($1==cohort) print NR }' )
    fgrep -w "HEALTHY_CONTROL" refs/HPOs_by_metacohort.table.tsv | cut -f$cidx
    echo -e "cnvs/$cohort.rCNV.strict_noncoding.bed.gz"
  done | paste -s
done < <( fgrep -v mega refs/rCNV_metacohort_list.txt | cut -f1 ) \
> metacohorts.input.tsv
/opt/rCNV2/data_curation/genome_annotations/count_cnvs_per_track.py \
  --cohorts metacohorts.input.tsv \
  --hpo-list refs/rCNV2.hpos_by_severity.developmental.list \
  --track-stats ${prefix}.stats.tsv \
  --frac-overlap ${min_element_overlap} \
  --norm-by-samplesize \
  --conditional-exclusion refs/GRCh37.200kb_bins_10kb_steps.raw.cohort_exclusion.bed.gz \
  --counts-by-gene-context unconstrained_cnv_counts.tsv.gz \
  --outfile ${prefix}.stats.with_counts.tsv.gz \
  --gzip 


# # Dev code:
# stats=${prefix}.stats.with_counts.tsv.gz
# merged_tracklist=test.annotations.list


# Burden meta-analysis of cohort CNV counts for each track
/opt/rCNV2/data_curation/genome_annotations/trackwise_cnv_burden_meta_analysis.R \
  --cutoff ${p_cutoff} \
  --signif-tracks ${prefix}.signif_paths_and_tracks.list \
  --volcano ${prefix}.track_volcanos.png \
  ${stats} \
  ${prefix}.burden_stats.tsv
gzip -f ${prefix}.burden_stats.tsv


# # Dev code:
# prefix="crb_clustering_test"
# min_prop_tracks_per_crb=0.01304348
# min_prop_track_representation=0.01304348
# clustering_neighborhood_dist=10000
# min_crb_separation=10000
# max_crb_size=500000
# blacklist_buffer=100000
# blacklist_buffer_min_size=100000
# contig=3


# Copy & index all final curated significant tracks locally
mkdir sig_tracks/
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/genome_annotations/significant_tracks/*.curated.bed.gz \
  sig_tracks/
find sig_tracks/ -name "*.curated.bed.gz" \
| xargs -I {} tabix -f {}

# Subset genome file to contig of interest
awk -v FS="\t" -v OFS="\t" -v contig=${contig} \
  '{ if ($1==contig) print $0 }' refs/GRCh37.genome \
> ${contig}.genome

# Create whitelist for CRB clustering
gsutil -m cp \
  ${rCNV_bucket}/cleaned_data/genes/gencode.v19.canonical.pext_filtered.gtf.gz* \
  ${rCNV_bucket}/cleaned_data/genes/gene_lists/gnomad.v2.1.1.likely_unconstrained.genes.list \
  ./
tabix gencode.v19.canonical.pext_filtered.gtf.gz ${contig} \
| gzip -c \
> ${contig}.gtf.gz
/opt/rCNV2/data_curation/genome_annotations/make_crb_whitelist.py \
  --bgzip \
  ${contig}.gtf.gz \
  gnomad.v2.1.1.likely_unconstrained.genes.list \
  ${contig}.genome \
  ${contig}.crb_whitelist.bed.gz


# Cluster significant tracks into CRBs
/opt/rCNV2/data_curation/genome_annotations/build_crbs.py \
  --genome ${contig}.genome \
  --blacklist refs/GRCh37.segDups_satellites_simpleRepeats_lowComplexityRepeats.bed.gz \
  --blacklist refs/GRCh37.somatic_hypermutable_sites.200kb_clustered.bed.gz \
  --blacklist refs/GRCh37.Nmask.autosomes.bed.gz \
  --blacklist-buffer ${blacklist_buffer} \
  --blacklist-buffer-min-size ${blacklist_buffer_min_size} \
  --whitelist ${contig}.crb_whitelist.bed.gz \
  --prop-min-elements ${min_prop_tracks_per_crb} \
  --prop-min-tracks ${min_prop_track_representation} \
  --neighborhood-dist ${clustering_neighborhood_dist} \
  --min-crb-separation ${min_crb_separation} \
  --max-crb-size ${max_crb_size} \
  --crb-prefix "${prefix}_CRB" \
  --crb-outbed ${prefix}.${contig}.crbs.bed.gz \
  --element-outbed ${prefix}.${contig}.crb_elements.bed.gz \
  --bgzip \
  sig_tracks/*.curated.bed.gz

