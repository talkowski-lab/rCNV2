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

Weights for each gene _i_ were computed as w<sub>_i_</sub> = CDS<sub>_i_</sub> / 2 <sup>&Sigma;(CDS) - 1</sup>  

Where w<sub>_i_</sub> is the per-gene weight, CDS<sub>_i_</sub> is the fraction of CDS from gene _i_ overlapped by the CNV, and &Sigma;(CDS) is the sum of the fractions of CDS for all genes overlapped by the CNV.  

For clarity, two toy examples of this weighting process are highlighted below:  

#### Example #1  
| Gene | CNV overlap? | Pct of CDS overlapped by CNV | Weighted CNV count |  
| :--- | :--- | ---: | ---: |  
| Gene A | Yes | 100% | 0.5 |  
| Gene B | Yes | 100% | 0.5 |  

In Example #1 above, the CNV completely overlaps two genes and no others.  

Since both genes are fully overlapped by the CNV, the CNV count is divided evenly by 2<sup>(1.0 + 1.0)</sup>, since both genes are 100% overlapped by the CNV. Therefore, each gene receives a weighted count of 0.5.  

#### Example #2  
| Gene | CNV overlap? | Pct of CDS overlapped by CNV | Weighted CNV count |  
| :--- | :--- | ---: | ---: |  
| Gene A | No | 0 | 0 |  
| Gene B | Yes | 100% | 1 / 2<sup>2.5</sup> = 0.18 |  
| Gene C | Yes | 100% | 1 / 2<sup>2.5</sup> = 0.18 |  
| Gene D | Yes | 50% | 0.5 / 2<sup>2.5</sup> = 0.09 |  

In Example #2 above, the CNV overlaps three of four genes.  

Gene A is not overlapped, so is not assigned any weighted CNV count.  

Genes B, C, and D are all at least partially overlapped, but only Genes B and C are completely overlapped.  

Thus, genes B, C, and D are each assigned a weight proportional to the fraction of their CDS overlapped by the CNV (CDS<sub>_i_</sub> / 2<sup>(1.0 + 1.0 + 0.5)</sup> =  CDS<sub>_i_</sub> / 2<sup>2.5</sup> = CDS<sub>_i_</sub> / 5.66)


### 2. Calculate burden statistics between cases & controls  

Once wighted CNVs were tallied in cases and controls for each window, we next compared the ratios of CNV carriers between cases and controls using the following procedure:  



The code to perform this step is contained in `window_burden_test.R`.  

#### Output files  

For each combination of phenotype group, metacohort, and CNV type, the following files are generated:  

1. `$metacohort.$hpo.rCNV.$CNV_type.sliding_window.stats.bed.gz`: a bgzipped BED file containing CNV counts and association statistics for each window  
2. `$metacohort.$hpo.rCNV.$CNV_type.sliding_window.manhattan.png`: a Manhattan plot of association statistics for each window
3. `$metacohort.$hpo.rCNV.$CNV_type.sliding_window.qq.png`: a QQ plot of observed P-values compared to expected P-values under a uniform null  
4. `$metacohort.$hpo.rCNV.$CNV_type.sliding_window.manhattan_with_qq.png`: a two-panel composite plot combining plots `2` and `3`, above  

Furthermore, for each pair of phenotype group & metacohort, two additional files are generated:  
5. `$metacohort.$hpo.rCNV.sliding_window.miami.png`: a Miami plot of association statistics for each window, with duplications above and deletions below the x-axis  
6. `$metacohort.$hpo.rCNV.sliding_window.miami_with_qq.png`: a multi-panel composite plot of association statistics, combining a Miami plot with QQ plots for deletions and duplications separately  

All `.png` plots are annotated with 54 known genomic disorder loci (adapted from [Owen _et al._, _BMC Genomics_, 2018](https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-018-5292-7)), for reference.  



---  

#### _A note on data curation_  

The information presented on this page references various curated datasets.  

The curation of these datasets is documented elsewhere in this repository.  

Please see the README available in [the `data_curation/` subdirectory](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/).  
