# Noncoding Association Analyses  

We compared rare CNV (rCNV) counts between cases and controls for clusters of _cis_-regulatory elements (_cis_-regulatory blocks, or CRBs) across all 22 autosomes. This process is described below.  

The code to reproduce these analyses is contained in `noncoding_association_analysis.sh`.  

In practice, this analysis was parallelized on [FireCloud/Terra](https://terra.bio) using `noncoding_association_analysis.wdl`.  

---  

## Preprocessing  

There were two preprocessing steps used exclusively for these noncoding association analyses:

1. Noncoding CNV filtering ([described here](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#noncoding-subsets))  
2. CRB definition ([described here](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/genome_annotations))  

---  

## Noncoding association test procedure  

We executed a standardized procedure to conduct CRB-based burden tests for each [phenotype group](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/).  

The steps of this procedure are described below:  

### 1. Count rCNVs per CRB  

### 2. Calculate burden statistics between cases & controls  

#### Output files  

### 3. Combine association statistics across metacohorts  

#### 4. Determining & calibrating exome-wide significance threshold  

#### Output files  

#### 5. Refine correlated associations  

#### 6. Reporting of final association statistics  

---  

#### _A note on data curation_  

The information presented on this page references various curated datasets.  

The curation of these datasets is documented elsewhere in this repository.  

Please see the README available in [the `data_curation/` subdirectory](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/).  
