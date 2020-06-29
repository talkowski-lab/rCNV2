# Genome Annotation Database Curation  

We curated a large database of genome annotation tracks for association testing. This process is described below.  

All commands required to curate the datasets as descrbed below are contained in `curate_annotations.sh`.  

In practice, this analysis was parallelized on [FireCloud/Terra](https://terra.bio) using `curate_annotations.wdl`.  

## Genome annotation manifest  

We aimed to comprehensively assess all major repositories of genome annotations.

In total, we evaluated _TBD_ genome annotation tracks across all sources.  

The following table outlines all sources included in this analysis:

| Source | Number of tracks | Description | Citation | Website |  
| :--- | ---: | :--- | :--- | :--- |  
| Roadmap ChromHMM | 1,764 | Chromatin states inferred by the extended 18-way ChromHMM model across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) | [Link](https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html) |  
| ENCODE DNA Accessibility | 1,228 | All DNA accessiblity assays meeting our [ENCODE data inclusion criteria](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/genome_annotations#encode-data-inclusion-criteria) available for download in broad- or narrow-peak format from the ENCODE Data Portal | ENCODE Data Portal [(Davis _et al._, _Nucleic Acids Res._, 2018)](https://pubmed.ncbi.nlm.nih.gov/29126249/) (accessed June 2020) | [Link](https://www.encodeproject.org/matrix/) |  
| ENCODE Histone Modifications | 2,977 | All histone modification ChIP-seq tracks meeting our [ENCODE data inclusion criteria](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/genome_annotations#encode-data-inclusion-criteria) available for download in broad- or narrow-peak format from the ENCODE Data Portal | ENCODE Data Portal [(Davis _et al._, _Nucleic Acids Res._, 2018)](https://pubmed.ncbi.nlm.nih.gov/29126249/) (accessed June 2020) | [Link](https://www.encodeproject.org/matrix/) |  
| ENCODE Transcription Factor Binding Sites | 2,859 | All transcription factor (or otherwise DNA-binding protein) ChIP-seq tracks meeting our [ENCODE data inclusion criteria](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/genome_annotations#encode-data-inclusion-criteria) available for download in broad- or narrow-peak format from the ENCODE Data Portal | ENCODE Data Portal [(Davis _et al._, _Nucleic Acids Res._, 2018)](https://pubmed.ncbi.nlm.nih.gov/29126249/) (accessed June 2020) | [Link](https://www.encodeproject.org/matrix/) |  
| ENCODE Transcription | 458 | A subset of transcription assays (CAGE, RAMPAGE, small RNAseq, total RNAseq, polyA RNAseq) meeting our [ENCODE data inclusion criteria](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/genome_annotations#encode-data-inclusion-criteria) available for download in BED format from the ENCODE Data Portal | ENCODE Data Portal [(Davis _et al._, _Nucleic Acids Res._, 2018)](https://pubmed.ncbi.nlm.nih.gov/29126249/) (accessed June 2020) | [Link](https://www.encodeproject.org/matrix/) |  
| ENCODE TAD boundaries | 30 | TAD boundaries (defined as the start and end coordinates for each TAD ± 5kb) from 30 samples meeting our [ENCODE data inclusion criteria](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/genome_annotations#encode-data-inclusion-criteria) available for download from the ENCODE Data Portal | ENCODE Data Portal [(Davis _et al._, _Nucleic Acids Res._, 2018)](https://pubmed.ncbi.nlm.nih.gov/29126249/) (accessed June 2020) | [Link](https://www.encodeproject.org/matrix/) |  
| EnhancerAtlas 2.0 | 197 | Enhancer predictions in 197 human cell lines & tissues | EnhancerAtlas 2.0 [(Gao _et al._, _Nucleic Acids Res._, 2019)](https://doi.org/10.1093/nar/gkz980) | [Link](http://www.enhanceratlas.org/indexv2.php) |  
| HACER | 289 | Active enhancer predictions in human cell lines & tissues based on PRO-seq, GRO-seq, or CAGE data | HACER [(Wang _et al._, _Nucleic Acids Res._, 2019)](https://doi.org/10.1093/nar/gky864) | [Link](http://bioinfo.vanderbilt.edu/AE/HACER/index.html) |  
| SEdb | 1,082 | Super enhancer and typical enhancer predictions from 541 human cell lines and tissues | SEdb v1.03 [(Jiang _et al._, _Nucleic Acids Res._, 2018)](https://doi.org/10.1093/nar/gky1025) | [Link](http://www.licpathway.net/sedb/index.php) |  
| dbSUPER | 99 | Super enhancers from 99 human cell lines and tissues | dbSUPER [(Khan _et al._, _Nucleic Acids Res._, 2016)](https://academic.oup.com/nar/article/44/D1/D164/2502575) | [Link](http://asntech.org/dbsuper/index.php) |  
| VISTA | 1 | Experimentally-validated mammalian enhancers | VISTA Browser [(Visel _et al._, _Nucleic Acids Res._, 2007)](https://dx.doi.org/10.1093%2Fnar%2Fgkl822) (accessed June 2020) | [Link](https://enhancer.lbl.gov/frnt_page_n.shtml) |  
| DENdb | 15 | Enhancer predictions from 15 human cell lines | DENdb [(Ashoor _et al._, _Database_, 2015)](doi:10.1093/database/bav085) | [Link](https://www.cbrc.kaust.edu.sa/dendb/index.php) |  
| SEA | 143 | Super enhancer predictions from 143 human cell lines and tissues (mapped back to hg19 using liftOver with minimum 75% match) | SEA v3.0 [(Ashoor _et al._, _Database_, 2015)](doi:10.1093/database/bav085) | [Link](http://sea.edbc.org/) |  
| FANTOM Enhancers | 115 | Enhancer predictions for human tissues and cell types from the FANTOM5 consortium | FANTOM5 [(Andersson _et al._, _Nature_, 2014)](https://www.nature.com/articles/nature12787) | [Link](http://enhancer.binf.ku.dk/presets/#download_view_div) |  
| PsychENCODE | 7 | Selected "derived" datasets from PsychENCODE Integrated Analysis Package, including cortex enhancers, transcriptionally active regions, TAD boundaries, and H3k27ac peaks | PsychEncode [(Wang _et al._, _Science_, 2018)](https://science.sciencemag.org/content/362/6420/eaat8464) | [Link](http://resource.psychencode.org/#Derived) |  


#### ENCODE data inclusion criteria  

To be considered in this analysis, ENCODE tracks had to meet all of the following criteria:  
1. Available for download in BED format
2. Aligned to hg19
3. Sample from an unperturbed human cell type, tissue, or cell line
4. No ENCODE data audit errors (color code: red) or non-compliances (color code: orange)
5. Marked as "released" status (_i.e._, not archived and/or retracted)  
6. Must be peaks-only (no background sample)  

Manifests for all ENCODE tracks considered in this analysis are stored in the following Google Bucket (note: requires permissions):  
```
$ gsutil ls gs://rcnv_project/cleaned_data/genome_annotations/manifests/

gs://rcnv_project/cleaned_data/genome_annotations/manifests/encode.dna_accessibility.manifest.tsv.gz
gs://rcnv_project/cleaned_data/genome_annotations/manifests/encode.hic_tads.manifest.tsv.gz
gs://rcnv_project/cleaned_data/genome_annotations/manifests/encode.histone_mods.manifest.tsv.gz
gs://rcnv_project/cleaned_data/genome_annotations/manifests/encode.tfbs.manifest.tsv.gz
gs://rcnv_project/cleaned_data/genome_annotations/manifests/encode.transcription.manifest.tsv.gz
```

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

For each track, we counted the number of [loose noncoding rare CNVs](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#noncoding-subsets) that completely overlapped at least one element per track.  

These counts of CNVs per track were split by [metacohort](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#case-control-metacohorts), CNV type, and [case/control status](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype#hpo-terms-per-metacohort).  

After tabulating counts of CNVs per track, we conducted a fixed-effects meta-analysis using the inverse-variance weighted Z-score method for each track with saddlepoint approximation applied across all tracks per cohort.  

We considered any track with a Benjamini-Hochberg FDR q < 0.1 to have sufficient evidence for possible disease relevance to be included for subsequent association testing.  

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

Finally, pairs of CRBs within ±10kb were merged, and CRBs with ≥30% coverage by our [CNV blacklist regions](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#curation-steps-rare-cnvs) were excluded.

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
