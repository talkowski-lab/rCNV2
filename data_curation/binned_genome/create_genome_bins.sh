#!/usr/bin/env bash

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Segment the reference into sliding windows and add key annotations


# Launch docker image
docker run --rm -it talkowski/rcnv


# Download necessary references (note: requires permissions)
gcloud auth login
mkdir refs/
gsutil cp -r gs://rcnv_project/refs/GRCh37.Nmask.bed.gz refs/
gsutil cp -r gs://rcnv_project/refs/GRCh37.autosomes.genome refs/
gsutil cp -r gs://rcnv_project/refs/Affy_UKBB_axiom_probes.bed.gz* refs/


# Prep reference fasta
# wget ftp://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
# gunzip Homo_sapiens.GRCh37.dna.primary_assembly.fa.gz
# samtools faidx Homo_sapiens.GRCh37.dna.primary_assembly.fa
# gsutil cp Homo_sapiens.GRCh37.dna.primary_assembly.fa gs://rcnv_project/GRCh37_ref_build/GRCh37.primary_assembly.fa
# gsutil cp Homo_sapiens.GRCh37.dna.primary_assembly.fa.fai gs://rcnv_project/GRCh37_ref_build/GRCh37.primary_assembly.fa.fai
gsutil cp -r gs://rcnv_project/GRCh37_ref_build/* refs/


# Set desired bin size & step size resolution (in kb)
binsize=100
stepsize=10


# Create sliding windows
athena make-bins -z \
	-x refs/GRCh37.Nmask.bed.gz \
  -s "$stepsize"000 \
	--buffer "$binsize"000 \
	refs/GRCh37.autosomes.genome \
	"$binsize"000 \
	GRCh37."$binsize"kb_bins_"$stepsize"kb_steps.raw.bed.gz
gsutil cp GRCh37."$binsize"kb_bins_"$stepsize"kb_steps.raw.bed.gz \
  gs://rcnv_project/cleaned_data/binned_genome/


# Annotate bins
athena annotate-bins -z \
  --chroms 20 \
  -t refs/Affy_UKBB_axiom_probes.bed.gz -a count -n ukbbAxiom_probes \
  -u decodeSexAveraged -a map-mean -n mean_recomb_rate \
  -u decodeSexAveraged -a map-max -n max_recomb_rate \
  -u wgEncodeDukeMapabilityUniqueness20bp -a map-mean -n uniqueness_20mer_mean \
  -u wgEncodeDukeMapabilityUniqueness35bp -a map-mean -n uniqueness_35mer_mean \
  -u snp151Common -a count-unique -n dbsnp151_common \
  -u snpArrayAffy6 -a count-unique -n affy6_probes \
  -u snpArrayAffy6SV -a count-unique -n affy6sv_probes \
  -u snpArrayAffy5 -a count-unique -n affy5_probes \
  -u snpArrayIlluminaHumanOmni1_Quad -a count-unique -n illOmni_probes \
  -u snpArrayIllumina1M -a count-unique -n ill1M_probes \
  -u snpArrayIllumina650 -a count-unique -n ill650_probes \
  -u snpArrayIllumina550 -a count-unique -n ill550_probes \
  -u agilentCgh1x244k -a count -n agilent244k_probes \
  -u agilentCgh4x180k -a count -n agilent180k_probes \
  -u agilentCgh2x105k -a count -n agilent105k_probes \
  -u rmsk:genoName,genoStart,genoEnd -a coverage -n rmsk_all \
  -u rmsk:genoName,genoStart,genoEnd,"repClass=LINE" -a coverage -n rmsk_LINE \
  -u rmsk:genoName,genoStart,genoEnd,"repClass=SINE" -a coverage -n rmsk_SINE \
  -u rmsk:genoName,genoStart,genoEnd,"repClass=LTR" -a coverage -n rmsk_LTR \
  -u rmsk:genoName,genoStart,genoEnd,"repClass=DNA" -a coverage -n rmsk_DNA_repeats \
  -u rmsk:genoName,genoStart,genoEnd,"repClass=Simple_repeat" -a coverage -n rmsk_simple_repeats \
  -u rmsk:genoName,genoStart,genoEnd,"repClass=Low_complexity" -a coverage -n rmsk_low_complexity_repeats \
  -u rmsk:genoName,genoStart,genoEnd,"repClass=Satellite" -a coverage -n rmsk_satellites \
  -u nestedRepeats -a coverage -n nested_repeats \
  -u microsat -a coverage -n microsatellites \
  -u genomicSuperDups -a coverage -n segdups \
  -u genomicSuperDups:fracMatch -a map-max -n max_segdup_identity \
  -u simpleRepeat -a coverage -n simple_repeats \
  -u chainSelf:tName,tStart,tEnd -a coverage -n self_chain \
  -u chainSelf:tName,tStart,tEnd,normScore -a map-max -n max_self_chain_score \
  --fasta refs/GRCh37.primary_assembly.fa \
  --ucsc-ref hg19 \
  GRCh37."$binsize"kb_bins_"$stepsize"kb_steps.raw.bed.gz \
  GRCh37."$binsize"kb_bins_"$stepsize"kb_steps.annotated.bed.gz


# Decompose annotated bins to top 10 PCs
athena eigen-bins -z \
  -e 10 \
  --sqrt-transform max_recomb_rate \
  --sqrt-transform mean_recomb_rate \
  --log-transform affy6_probes \
  --log-transform affy5_probes \
  --log-transform illOmni_probes \
  --log-transform ill1M_probes \
  --log-transform ill650_probes \
  --log-transform ill550_probes \
  --log-transform ukbbAxiom_probes \
  --log-transform rmsk_LINE \
  --log-transform rmsk_SINE \
  --log-transform rmsk_LTR \
  --log-transform rmsk_DNA_repeats \
  --log-transform rmsk_simple_repeats \
  --log-transform rmsk_low_complexity_repeats \
  --log-transform segdups \
  --log-transform max_segdup_identity \
  --log-transform simple_repeats \
  --log-transform self_chain \
  --stats GRCh37."$binsize"kb_bins_"$stepsize"kb_steps.eigenfeature_stats.txt \
  GRCh37."$binsize"kb_bins_"$stepsize"kb_steps.annotated.bed.gz \
  GRCh37."$binsize"kb_bins_"$stepsize"kb_steps.annotated.eigen.bed.gz


