# Gene-Based Analyses  

We compared CNV counts between cases and controls for each of 19,346 canonical, protein-coding genes across all 22 autosomes. This process is described below.  

The code to reproduce these analyses is contained in `gene_burden_analysis.sh`.  

In practice, this analysis was parallelized on [FireCloud/Terra](https://terra.bio) using `gene_burden_analysis.wdl`.  

## Gene-based burden test procedure

We executed a standardized procedure to conduct gene-based burden tests for each [phenotype group](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/).  

The steps of this procedure are described below:  

### 1. Counting weighted CNVs per gene  

For each [metacohort](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/phenotype/), we intersected rare CNVs against [all canonical, autosomal, protein-coding genes](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/gene/) separately for cases and controls.  

We conducted this procedure a total of three times per phenotype group & metacohort: once each for deletions, duplications, and all CNVs (deletions + duplications).  

The code to perform this step is contained in `count_cnvs_per_gene.py`.  

We considered a CNV to overlap a gene if it overlapped any exons from that gene.  

For each overlapping CNV-gene pair, we computed the total fraction of exonic bases (CDS) from the gene overlapped by the CNV.  

Finally, for each CNV, we distributed weights across all overlapping genes proportionally to the sum of the CDS overlapped per gene.  

For clarity, two toy examples of this weighting process are highlighted below:  

#### Example #1  
| Gene | CNV overlap? | Pct of CDS overlapped by CNV | Weighted CNV count |  
| :--- | :--- | ---: | ---: |  
| Gene A | Yes | 100% | 0.5 |  
| Gene B | Yes | 100% | 0.5 |  

In Example #1 above, the CNV completely overlaps two genes and no others.  

Since both genes are fully overlapped by the CNV, the CNV count is divided evenly between them. Each gene receives a weighted count of 0.5.  

#### Example #2  
| Gene | CNV overlap? | Pct of CDS overlapped by CNV | Weighted CNV count |  
| :--- | :--- | ---: | ---: |  
| Gene A | No | 0 | 0 |  
| Gene B | Yes | 50% | 0.25 |  
| Gene C | Yes | 100% | 0.5 |  
| Gene D | Yes | 50% | 0.25 |  

In Example #2 above, the CNV overlaps three of four genes.  

Gene A is not overlapped, so is not assigned any weighted CNV count.  

Genes B, C, and D are all at least partially overlapped, but only Gene C is completely overlapped.  

Thus, genes B, C, and D are each assigned a weight proportional to the fraction of their CDS overlapped by the CNV. This is computed as the fraction of CDS overlapped divided by the sum of CDS overlap fractions across all genes.  

---  

#### _A note on data curation_  

The information presented on this page references various curated datasets.  

The curation of these datasets is documented elsewhere in this repository.  

Please see the README available in [the `data_curation/` subdirectory](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/).  
