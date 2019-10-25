# Genome-Wide Bin Curation  

## Genome-Wide Bin Curation  

We segmented the genome into a set of regularized sliding windows for association testing. This process is described below.  

All commands executed to filter the CNV data as descrbed below are contained in `create_genome_bins.sh`.  

In practice, the commands in `create_genome_bins.sh` were parallelized in [FireCloud/Terra](https://portal.firecloud.org) with `create_genome_bins.wdl`.  

### Bin files  

All genome-wide bin files discussed below are stored in a protected Google Cloud bucket:  
```
$ gsutil ls gs://rcnv_project/cleaned_data/binned_genome/

gs://rcnv_project/cleaned_data/binned_genome/GRCh37.200kb_bins_10kb_steps.raw.bed.gz
gs://rcnv_project/cleaned_data/binned_genome/GRCh37.200kb_bins_10kb_steps.annotated.bed.gz
gs://rcnv_project/cleaned_data/binned_genome/GRCh37.200kb_bins_10kb_steps.annotated.eigen.bed.gz
```

### Bin creation & annotation

We created sliding windows for all autosomes at 200kb resolution and 10kb step size, and excluded any bins within ±200kb of any N-masked sequence or known somatically hypermutable site (as applied in [Collins\*, Brand\*, _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/biorxiv/early/2019/03/14/578674)).  

The window size of 200kb was selected to approximately match the median size of rare CNVs for most cohorts following [our CNV filtering protocol](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV/).  

After filtering, we retained a final set of 256,329 bins for analysis.  

To control for technical and genomic covariates, we annotated all bins against a suite of features, then performed Eigendecomposition to control for the inherent correlation structure of most genomic annotation tracks.  

We included the following annotations:  

| Track | Source | File / Table | Athena function(s) | Transformation |  
| :--- | :---- | :--- | :--- | :--- |  
| GC Content | GRCh37 reference | `GRCh37.primary_assembly.fa` | `fasta` | None |  
| Sex-averaged recombination rate | UCSC | `decodeSexAveraged` | `map-mean`, `map-max` | `sqrt(x)` |  
| Mean 20mer uniqueness | UCSC | `wgEncodeDukeMapabilityUniqueness20bp` | `map-mean` | None |  
| Mean 35mer uniqueness | UCSC | `wgEncodeDukeMapabilityUniqueness35bp` | `map-mean` | None |  
| Common SNP density | UCSC | `snp151Common` | `count-unique` | None |  
| Affy 6.0 probe density | UCSC | `snpArrayAffy6` | `count-unique` | `log(x+0.01max(x))` |  
| Affy 6.0 SV probe density | UCSC | `snpArrayAffy6SV` | `count-unique` | None |  
| Affy 5.0 probe density | UCSC | `snpArrayAffy5` | `count-unique` | `log(x+0.01max(x))` |  
| Illumina Omni probe density | UCSC | `snpArrayIlluminaHumanOmni1_Quad` | `count-unique` | `log(x+0.01max(x))` |  
| Illumina 1M probe density | UCSC | `snpArrayIllumina1M` | `count-unique` | `log(x+0.01max(x))` |  
| Illumina 650 probe density | UCSC | `snpArrayIllumina650` | `count-unique` | `log(x+0.01max(x))` |  
| Illumina 550 probe density | UCSC | `snpArrayIllumina550` | `count-unique` | `log(x+0.01max(x))` |  
| Agilent 244k aCGH probe density | UCSC | `agilentCgh1x244k` | `count` | None |  
| Agilent 180k aCGH probe density | UCSC | `agilentCgh4x180k` | `count` | None |  
| Agilent 105k aCGH probe density | UCSC | `agilentCgh2x105k` | `count` | None |  
| Affy UKBB Axiom probe density | [Thermo-Fisher](https://www.thermofisher.com/order/catalog/product/902502) | `Axiom_UKB_WCSG.na35.annot.csv` | `count` | `log(x+0.01max(x))` |  
| RepeatMasker (all) | UCSC | `rmsk` | `coverage` | None |  
| RepeatMasker (LINEs) | UCSC | `rmsk` | `coverage` | `log(x+0.01max(x))` |  
| RepeatMasker (SINEs) | UCSC | `rmsk` | `coverage` | `log(x+0.01max(x))` |  
| RepeatMasker (LTRs) | UCSC | `rmsk` | `coverage` | `log(x+0.01max(x))` |  
| RepeatMasker (DNA repeats) | UCSC | `rmsk` | `coverage` | `log(x+0.01max(x))` |  
| RepeatMasker (Simple repeats) | UCSC | `rmsk` | `coverage` | `log(x+0.01max(x))` |  
| RepeatMasker (Low-complexity repeats) | UCSC | `rmsk` | `coverage` | `log(x+0.01max(x))` |  
| RepeatMasker (Satellites) | UCSC | `rmsk` | `coverage` | None |  
| Nested repeats | UCSC | `nestedRepeats` | `coverage` | None |  
| Microsatellites | UCSC | `microsat` | `coverage` | None |  
| Segmental duplications | UCSC | `genomicSuperDups` | `coverage` | `log(x+0.01max(x))` |  
| Max seg. dup. identity | UCSC | `genomicSuperDups:fracMatch` | `map-max` | `log(x+0.01max(x))` |  
| Simple repeats | UCSC | `simpleRepeat` | `coverage` | `log(x+0.01max(x))` |  
| Self-chain | UCSC | `chainSelf` | `coverage` | `log(x+0.01max(x))` |  
| Max self-chain score | UCSC | `chainSelf:normScore` | `coverage` | None |  

Following annotation, we collapsed the feature correlation structure with Eigendecomposition, and retained the top 10 components (ranked by % variance explained), labeled `Eigenfeatures`.  

### Dependencies  

#### Athena  
Several commands used in the creation and processing of binned genome files require Athena, a toolkit currently in development by the Talkowski Lab.  

Athena is not currently publicly available. For now, the current Athena build is installed into the rCNV Docker, and can be accessed as follows:  
```
$ docker pull talkowski/rcnv

$ docker run --rm -it talkowski/rcnv

(rcnv) root@hash:/# athena --help
Usage: athena [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  annotate-bins  Annotate bins
  count-sv       Intersect SV and bins
  eigen-bins     Eigendecomposition of annotations
  make-bins      Create sequential bins
  query          Mutation rate lookup
  vcf-filter     Filter an input VCF
  vcf-stats      Get SV size & spacing
```
