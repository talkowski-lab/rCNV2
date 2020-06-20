# Genome Annotation Database Curation  

We curated a large database of genome annotation tracks for association testing. This process is described below.  

All commands required to curate the datasets as descrbed below are contained in `curate_annotations.sh`.  

In practice, this analysis was parallelized on [FireCloud/Terra](https://terra.bio) using `curate_annotations.wdl`.  

## Genome annotation manifest  

We aimed to comprehensively assess all major repositories of genome annotations.

In total, we evaluated TBD genome annotation tracks across all sources.  

The following table outlines major sources considered in this analysis:

| Source | Number of tracks | Description | Citation | Website |  
| :--- | ---: | :--- | :--- | :--- |  
| Roadmap Epigenomics (ChromHMM only) | 1,764 | Chromatin states inferred by the extended 18-way ChromHMM model across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) | [Link](https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html) |  
| ENCODE Data Portal | TBD | All hg19-aligned BED files from any human tissue or cell line available for download from the ENCODE Data Portal (accessed June 2020) | ENCODE Data Portal [(Davis _et al._, _Nucleic Acids Res._, 2018)](https://pubmed.ncbi.nlm.nih.gov/29126249/) | [Link](https://www.encodeproject.org/matrix/) |  


## Track curation procedure  

We subjected each annotation track to an identical set of curation steps. 

This procedure was performed using `curate_track.py`, and is described below.

For each track, we:
1. Excluded any elements covered at least 50% by the same set of blacklists used during [CNV curation](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#curation-steps-rare-cnvs);  
2. Standardized contig nomenclature to be consistent with GRCh37; 
3. Merged overlapping elements; 
4. Excluded merged elements <5bp or >200kb in size; and
5. Restricted to autosomes.  

After curation, we computed the following statistics for each track:
* Number of elements in track
* Minimum, mean, median, and maximum element size
* Total nonredundant nucleotides covered by all elements


## Identifying annotation classes with burdens of noncoding rCNVs  

To reduce our search space for association testing, we next wanted to determine which annotation tracks had evidence of dosage sensitivity in disease.

We restricted subsequent analyses to a subset of tracks following the procedure described below:

For each track, we counted the number of strictly noncoding rare CNVs that completely overlapped at least one element per track.  

These counts of CNVs per track were split by [metacohort](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#case-control-metacohorts), CNV type, and [case/control status](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype#hpo-terms-per-metacohort).  

After tabulating counts of CNVs per track, we conducted a fixed-effects meta-analysis using the inverse-variance weighted Z-score method for each track with saddlepoint approximation applied across all tracks per cohort.  

We considered any track with a meta-analysis P < 0.05 to have sufficient evidence for possible disease relevance to be included for subsequent association testing.  

#### Output files  

From this process, we produced the following:  
1. `rCNV.burden_stats.tsv.gz`: a tab-delimited file containing all statistics for each annotation track evaluated, including number & distribution of elements, number of CNVs per meta-cohort, burden statistics per meta-cohort, and meta-analysis statistics
2. One BED file for each significant track.  

The data are available from the below Google storage bucket (note: requires permissions):
```
$ gsutil ls gs://rcnv_project/cleaned_data/genome_annotations/
gs://rcnv_project/cleaned_data/genome_annotations/rCNV.burden_stats.tsv.gz
gs://rcnv_project/cleaned_data/genome_annotations/significant_tracks/
```

### Cluster significant annotations into CRBs  

Given that many genome annotations are correlated, we next clustered elements into CRBs for association testing.  

This process was restricted to the subset of genome annotation tracks with statistical evidence for relevance in disease, as [described above](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/genome_annotations#identifying-annotation-classes-with-burdens-of-noncoding-rcnvs).  

All elements from significant annotations were clustered into CRBs per-chromosome using DBSCAN (a density-based clustering algorithm) with the following parameters:
* Neighborhood distance < 10kb
* Minimum number of elements per cluster proportional to 10% of the total number of annotation tracks being considered (rounding down)

After clustering, CRBs were assigned the minimum start and maximum end coordinate across all their constituent elements.  

Finally, pairs of CRBs within ±10kb were merged, and CRBs with ≥50% coverage by our [CNV blacklist regions](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#curation-steps-rare-cnvs) were excluded.

#### Output files  

This process produced the following files:  
1. `rCNV.crbs.bed.gz`: a BED file with one entry per final CRB
2. `rCNV.crb_elements.bed.gz`: a BED file with one entry per element belonging to one of the final CRBs

The data are available from the below Google storage bucket (note: requires permissions):
```
$ gsutil ls gs://rcnv_project/cleaned_data/genome_annotations/*bed.gz*
gs://rcnv_project/cleaned_data/genome_annotations/rCNV.crb_elements.bed.gz
gs://rcnv_project/cleaned_data/genome_annotations/rCNV.crb_elements.bed.gz.tbi
gs://rcnv_project/cleaned_data/genome_annotations/rCNV.crbs.bed.gz
gs://rcnv_project/cleaned_data/genome_annotations/rCNV.crbs.bed.gz.tbi
```
