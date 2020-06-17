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
4. Excluded merged elements <10bp or >100kb in size; and
5. Restricted to autosomes.  

After curation, we computed the following statistics for each track:
* Number of elements in track
* Minimum, mean, median, and maximum element size
* Total nonredundant nucleotides covered by all elements


## Identifying annotation classes with burdens of noncoding rCNVs  

To reduce our search space for association testing, we next wanted to determine which annotation tracks had evidence of dosage sensitivity in disease.

We restricted subsequent analyses to a subset of tracks following the procedure described below:

For each track, we counted the number of strictly noncoding rare CNVs that overlapped at least one element per track.  

These counts of CNVs per track were split by [metacohort](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#case-control-metacohorts), CNV type, and [case/control status](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype#hpo-terms-per-metacohort).  

After tabulating counts of CNVs per track, we conducted a fixed-effects meta-analysis for each track using the same approach as for our [sliding window](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/sliding_windows#3-combine-association-statistics-across-metacohorts) or [gene-based association tests](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/genes#3-combine-association-statistics-across-metacohorts).  

We applied a Benjamini-Hochberg correction to the meta-analysis P-values across all tracks & CNV types, and considered any track with FDR q < 0.05 to have sufficient evidence for possible disease relevance to be included for subsequent association testing.  

### Cluster significant annotations into CRBs  

Given that many genome annotations are correlated, we next clustered elements into CRBs for association testing.  

This process was restricted to the subset of genome annotation tracks with statistical evidence for relevance in disease, as [described above](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/genome_annotations#identifying-annotation-classes-with-burdens-of-noncoding-rcnvs).  

TBD

