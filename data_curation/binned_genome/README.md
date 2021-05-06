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
```

### Bin creation

We created sliding windows for all autosomes at 200kb resolution and 10kb step size, and excluded any bins with ≥30% coverage by N-masked sequences or known somatically hypermutable site (as applied in [Collins\*, Brand\*, _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/biorxiv/early/2019/03/14/578674)).  

The window size of 200kb was selected to approximately match the median size of rare CNVs for most cohorts following [our CNV filtering protocol](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV/).  

After filtering, we retained a final set of 267,237 bins for analysis.  


### Control array probe counts  

To control for technical differences between cohorts introduced by different microarray platforms, we annotated all 200kb bins with the count of probes for every major array used to genotype control samples in any cohort in this study. The specifics of these probesets are [described here](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/other). 

