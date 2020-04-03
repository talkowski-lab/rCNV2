# CNV Data Curation

We aggregated existing microarray-based CNV calls from numerous sources. This process is described below.  

All commands executed to filter the CNV data as descrbed below are contained in `filter_CNV_data.sh`.  

In practice, the commands in `filter_CNV_data.sh` were parallelized in [FireCloud/Terra](https://portal.firecloud.org) with `filter_CNV_data.wdl`.  

## CNV data sources  

We aggregated CNV data from multiple sources, listed alphabetically below:  

| Dataset | Citation | PMID | Platform(s) | Build | Phenos | N Cases | N Ctrls |
| --- | :--- | :--- | :--- | :--- | :--- | ---: | ---: |
| BCH | [Talkowski _et al._, _Cell_ (2012)](https://www.sciencedirect.com/science/article/pii/S0092867412004114) | [22521361](https://www.ncbi.nlm.nih.gov/pubmed/22521361) | aCGH (?) | hg18 | Mixed | 3,591 | 0 |  
| CHOP | - | - | Mixed Illumina SNP genotyping platforms | hg19 | Mixed | 153,870<sup>1</sup> | 24,161<sup>1</sup> |
| Coe | [Coe _et al._, _Nat. Genet._ (2014)](https://www.nature.com/articles/ng.3092) | [25217958](https://www.ncbi.nlm.nih.gov/pubmed/25217958) | Cases: SignatureChip OS v2.0 (58%), SignatureChip OS v1.0 (34%), Other (8%); Controls: Affy 6.0 (100%) | hg19 | Developmental disorders | 29,083 | 11,256 |
| Cooper<sup>2</sup> | [Cooper _et al._, _Nat. Genet._ (2011)](https://www.nature.com/articles/ng.909) | [21841781](https://www.ncbi.nlm.nih.gov/pubmed/21841781) | Ill. 550k-610k (75%), Custom 1.2M (25%) | hg19 | Developmental disorders | 0<sup>2</sup> | 8,329 |
| Epi25k | [Niestroj _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/10.1101/651299v1) | - | Illumina GSA-MD v1.0 (100%) | hg19 | Epilepsy | 12,053<sup>1</sup> | 8,173<sup>1</sup> |  
| GDX | - | - | aCGH (?) | hg18 & hg19 | Mixed | 9,958 | 0 |
| IU<sup>3</sup> | - | - | CMA(?) | hg19(?) | Mixed | 1,576<sup>3</sup> | 0 |  
| PGC | [Marshall _et al._, _Nat. Genet._ (2017)](https://www.nature.com/articles/ng.3725) | [27869829](https://www.ncbi.nlm.nih.gov/pubmed/27869829) | Affy 6.0 (37%), Omni Express (31%), Omni Express Plus (12%), Other (20%) | hg18 | Schizophrenia | 21,094 | 20,277 |
| SSC<sup>4</sup> | [Sanders _et at._, _Neuron_ (2015)](https://www.sciencedirect.com/science/article/pii/S0896627315007734?) | [26402605](https://www.ncbi.nlm.nih.gov/pubmed/26402605) | Omni 1Mv3 (46%), Omni 2.5 (41%), Omni 1Mv1 (13%) | hg18 | Autism | 2,795 | 0<sup>4</sup> |
| SickKids | [Zarrei _et al._, _NPJ Genomic Medicine_ (2019)](https://www.nature.com/articles/s41525-019-0098-3) | [31602316](https://www.ncbi.nlm.nih.gov/pubmed/31602316) | Affy 6.0 (100%) | hg19 | NDDs | 2,689<sup>1</sup> | 0<sup>4</sup> |  
| TCGA | [Zack _et al._, _Nat. Genet._ (2013)](https://www.nature.com/articles/ng.2760) | [24071852](https://www.ncbi.nlm.nih.gov/pubmed/24071852) | Affy 6.0 (100%) | hg19 | Cancer | 0<sup>5</sup> | 8,670<sup>5</sup> |  
| TSAICG | [Huang _et al._, _Neuron_ (2017)](https://www.sciencedirect.com/science/article/pii/S0896627317305081) | [28641109](https://www.ncbi.nlm.nih.gov/pubmed/28641109) | OmniExpress (100%) | hg19 | Tourette Syndrome | 2,434 | 4,093 |
| UKBB | [Macé _et al._, _Nat. Comms._ (2017)](https://www.nature.com/articles/s41467-017-00556-x) | [28963451](https://www.ncbi.nlm.nih.gov/pubmed/28963451) | UKBB Affy Axiom (100%) | hg19 | Mixed | 54,071<sup>1</sup> | 375,800<sup>1</sup> |

#### Notes on raw CNV data   
1. Counts represent the number of samples retained after filtering outliers, described below.  
2. Only retained control samples from Cooper _et al._ All cases from Cooper _et al._ also appear in Coe _et al._  
3. Excluded samples from Indiana University (IU) cohort derived from buccal swab DNA, samples with known aneuploidies or large runs of homozygosity, and samples with no phenotypic indication specified.  
4. Only retained affected children from Sanders _et al._ and Zarrei _et al._, since all controls were first-degree relatives of affected cases.  
5. Only retained normal samples from tumor:normal pairs, and excluded any donors with known blood cancer.  

## Raw CNV data processing steps  

All CNV data native to hg18 was lifted over to hg19 using UCSC liftOver, requiring at least 50% of the original CNV to map successfully to hg19 in order to be retained.  

Some datasets required manual curation prior to inclusion. Where necessary, these steps are enumerated below:  

 * **CHOP**: CNVs were filtered on quality score ≥40 and CNV size ≥25kb while requiring at least 10 SNPs per CNV. After CNV filtering, samples with `LRR_SD` <0.25, >20 CNV calls, or SNP call rate <98% were excluded as outliers, as well as samples genotyped on arrays with <175k SNP probes or samples labeled as cancer or Down's Syndrome patients. Finally, we identified 19 loci with apparently platform-specific artifactual CNV pileups. CNVs covered ≥10% by any of these artifact regions were removed from the callset. Lists of CHOP-specific blacklisted loci for deletions and duplications are provided as [reference files](https://github.com/talkowski-lab/rCNV2/tree/master/refs).  
 * **Epi25k**: CNVs were filtered on ≥10 probes and ≥25kb. After CNV filtering, samples with >25 CNV calls were excluded as outliers.  
 * **SSC**: CNVs were filtered on pCNV ≤10<sup>-9</sup>, per recommendation of the authors.  
 * **SickKids**: CNVs were filtered on ≥25kb. After CNV filtering, samples with >80 CNV calls were excluded as outliers. Finally, we identified a single locus on chr12 that had CNVs only appearing in ADHD samples at 2.8% frequency; these CNVs were removed from the callset.  
 * **TCGA**: CNVs were filtered on ≥10 probes and ≥25kb. Deletions were required to have a mean log<sub>2</sub> intensity ≤-1 and duplications were required to have a mean log<sub>2</sub> intensity of ≥0.5849625.  
 * **UKBB**: CNVs were filtered on quality score ≥17 and CNV size ≥25kb. After CNV filtering, samples with >10 CNV calls were excluded as outliers as well as any samples with known malignant cancers or chromosomal disorders (e.g., Down's Syndrome or sex chromosome aneuploidies).  

#### CNV defragmentation  

Many array-based CNV calling algorithms infrequently result in fragmentation (i.e., over-segmentation) of large CNV calls.  

For the purposes of this study, fragmentation of CNV calls has the potential to bias association tests, as single individuals might be counted multiple times for a given locus or gene.  

Thus, we applied a standardized defragmentation step to raw CNV calls for all studies where this was possible.  

Defragmentation was performed with `defragment_cnvs.py` using `--max-dist 0.25`, which merges CNVs of the same type found in the same sample if their breakpoints are within ±25% of the size of their corresponding original CNV calls.  

For example, below is a 4.5Mb deletion reported as three deletion fragments in the `TCGA` cohort:

| chrom | start | end | CNV type |  
| ---: | ---: | ---: | :--- |
| 10 | 54,387,162 | 54,542,011 | DEL |  
| 10 | 54,544,685 | 55,877,827 | DEL |  
| 10 | 55,883,144 | 58,906,151 | DEL |  

Given that these three fragments are (i) all reported in the same sample and (ii) separated by just 2.7kb and 5.3kb, respectively, it is highly unlikely that these deletions are independent CNV events.  

The defragmentation algorithm as applied to these data would merge these three deletion fragments into a single CNV spanning from `54387162` to `58906151`.  

Four studies were unable to be defragmented due to inadequate sample-level information: `Coe` (_controls only_), `Cooper`, `PGC`, and `TSAICG`.

### Raw CNV callset properties  

The properties of each callset are listed below after initial data processing steps but prior to further filtering.  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| BCH | 3,591 | 5,159 | 1.44 | 204.0 kb | 1:1.27 | 0 | 0 | 0 | NA | NA |  
| CHOP | 153,870 | 641,938 | 4.17 | 87.0 kb | 1.05:1 | 24,161 | 108,559 | 4.49 | 87.0 kb | 1.86:1 |  
| Coe | 29,083 | 28,667 | 0.99 | 187.0 kb | 1:1.11 | 11,256 | 273,331 | 24.28 | 53.0 kb | 1.24:1 |  
| Cooper | 0 | 0 | 0 | NA | NA | 8,329 | 432,478 | 51.92 | 2.0 kb | 8.04:1 |  
| Epi25k | 12,053 | 92,412 | 7.67 | 97.0 kb | 1:1.54 | 8,173 | 64,036 | 7.84 | 100.0 kb | 1:1.33 |  
| GDX | 9,958 | 20,658 | 2.07 | 196.0 kb | 1:1.76 | 0 | 0 | 0 | NA | NA |  
| IU | 1,576 | 4,890 | 3.1 | 98.0 kb | 1.41:1 | 0 | 0 | 0 | NA | NA |  
| PGC | 21,094 | 42,096 | 2.0 | 80.0 kb | 1:1.04 | 20,277 | 40,464 | 2.0 | 78.0 kb | 1:1.03 |  
| SSC | 2,795 | 30,856 | 11.04 | 21.0 kb | 3.09:1 | 0 | 0 | 0 | NA | NA |  
| SickKids | 2,689 | 72,278 | 26.88 | 55.0 kb | 1.04:1 | 0 | 0 | 0 | NA | NA |  
| TCGA | 0 | 0 | 0 | NA | NA | 8,670 | 139,082 | 16.04 | 58.0 kb | 1:1.69 |  
| TSAICG | 2,434 | 3,541 | 1.45 | 91.0 kb | 1.01:1 | 4,093 | 5,834 | 1.43 | 91.0 kb | 1:1.08 |  
| UKBB | 54,071 | 140,082 | 2.59 | 88.0 kb | 5.62:1 | 375,800 | 963,956 | 2.57 | 87.0 kb | 5.68:1 |  

The information for this table was collected using `collect_cohort_stats.sh`, and is visualized below using `plot_cnv_stats_per_cohort.R`.  

![Raw CNV stats](https://storage.googleapis.com/rcnv_project/public/raw_cnv.stats.jpg)  

### Raw data access  

All raw CNV data files and their tabix indexes are stored in a protected Google Cloud bucket, here:  
```
$ gsutil ls gs://rcnv_project/raw_data/cnv/

gs://rcnv_project/raw_data/cnv/BCH.raw.bed.gz
gs://rcnv_project/raw_data/cnv/BCH.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/CHOP.raw.bed.gz
gs://rcnv_project/raw_data/cnv/CHOP.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/Coe.raw.bed.gz
gs://rcnv_project/raw_data/cnv/Coe.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/Cooper.raw.bed.gz
gs://rcnv_project/raw_data/cnv/Cooper.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/Epi25k.raw.bed.gz
gs://rcnv_project/raw_data/cnv/Epi25k.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/GDX.raw.bed.gz
gs://rcnv_project/raw_data/cnv/GDX.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/IU.raw.bed.gz
gs://rcnv_project/raw_data/cnv/IU.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/PGC.raw.bed.gz
gs://rcnv_project/raw_data/cnv/PGC.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/SSC.raw.bed.gz
gs://rcnv_project/raw_data/cnv/SSC.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/SickKids.raw.bed.gz
gs://rcnv_project/raw_data/cnv/SickKids.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/TCGA.raw.bed.gz
gs://rcnv_project/raw_data/cnv/TCGA.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/TSAICG.raw.bed.gz
gs://rcnv_project/raw_data/cnv/TSAICG.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/UKBB.raw.bed.gz
gs://rcnv_project/raw_data/cnv/UKBB.raw.bed.gz.tbi
```

Note that permissions must be granted per user prior to data access.  

## Curation Steps: Rare CNVs  

All raw CNV data was subjected to the same set of global filters:  
 * Restricted to autosomes
 * CNV size ≥ 100kb and ≤ 20Mb
 * Does not have substantial overlap<sup>1</sup> with a common CNV (`AF`>1%) in any population from WGS resolution in gnomAD-SV<sup>2</sup> ([Collins\*, Brand\*, _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/biorxiv/early/2019/03/14/578674))  
 * Does not have substantial overlap<sup>1</sup> with a common CNV (`AF`>1%) in any population from WGS resolution in the 1000 Genomes Project, Phase III ([Sudmant _et al._, _Nature_ (2015)](https://www.nature.com/articles/nature15394))  
 * Does not have substantial overlap<sup>1</sup> with a common CNV (`AF`>1%) from WGS resolution in the NIH Center for Common Disease Genetics SV callset ([Abel _et al._, _bioRxiv_ (2018)](https://www.biorxiv.org/content/10.1101/508515v1))  
 * Does not have substantial overlap<sup>1</sup> with other CNVs in at least 1% of all samples
 * Does not have substantial overlap<sup>1</sup> with other CNVs in at least 1% of all samples within the same dataset or in any of the other array CNV datasets (compared pairwise in serial)
 * Not substantially covered<sup>3</sup> by somatically hypermutable sites (as applied in [Collins\*, Brand\*, _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/biorxiv/early/2019/03/14/578674))  
 * Not substantially covered<sup>3</sup> by segmental duplications and/or simple, low-complexity, or satellite repeats  
 * Not substantially covered<sup>3</sup> by N-masked regions of the hg19 reference genome assembly  

#### Notes on curation  
1. "Substantial" overlap determined based on ≥50% reciprocal overlap using BEDTools ([Quinlan & Hall, _Bioinformatics_ (2010)](https://www.ncbi.nlm.nih.gov/pubmed/20110278)). For overlap-based comparisons, both breakpoints were required to be within ±100kb, and CNV type (DEL vs. DUP) was required to match.  
2. The version of gnomAD-SV used for this analysis (gnomAD-SV v2.1, non-neuro) included 8,342 samples without known neuropsychiatric disorders as [available from the gnomAD website](https://gnomad.broadinstitute.org/downloads/) and described in [Collins\*, Brand\*, _et al._](https://www.biorxiv.org/content/10.1101/578674v1)  
3. "Substantial" coverage determined based on ≥30% coverage per BEDTools coverage ([Quinlan & Hall, _Bioinformatics_ (2010)](https://www.ncbi.nlm.nih.gov/pubmed/20110278)).

### Rare CNV callset properties

The properties of each rare CNV callset are listed below after the above filtering steps.  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| BCH | 3,591 | 2,776 | 0.77 | 333.0 kb | 1:1.40 | 0 | 0 | 0 | NA | NA |  
| CHOP | 153,870 | 81,530 | 0.53 | 195.0 kb | 1:1.05 | 24,161 | 15,205 | 0.63 | 192.0 kb | 1.01:1 |  
| Coe | 29,083 | 17,643 | 0.61 | 264.0 kb | 1:1.17 | 11,256 | 8,141 | 0.72 | 174.0 kb | 1:1.60 |  
| Cooper | 0 | 0 | 0 | NA | NA | 8,329 | 3,841 | 0.46 | 171.0 kb | 1.41:1 |  
| Epi25k | 12,053 | 7,869 | 0.65 | 177.0 kb | 1:1.05 | 8,173 | 4,885 | 0.6 | 186.0 kb | 1:1.11 |  
| GDX | 9,958 | 9,706 | 0.97 | 264.0 kb | 1:1.32 | 0 | 0 | 0 | NA | NA |  
| IU | 1,576 | 1,476 | 0.94 | 208.0 kb | 1:1.38 | 0 | 0 | 0 | NA | NA |  
| PGC | 21,094 | 13,014 | 0.62 | 186.0 kb | 1:1.34 | 20,277 | 11,863 | 0.59 | 180.0 kb | 1:1.31 |  
| SSC | 2,795 | 1,906 | 0.68 | 191.0 kb | 1:1.14 | 0 | 0 | 0 | NA | NA |  
| SickKids | 2,689 | 3,526 | 1.31 | 177.0 kb | 1:1.32 | 0 | 0 | 0 | NA | NA |  
| TCGA | 0 | 0 | 0 | NA | NA | 8,670 | 4,372 | 0.5 | 177.0 kb | 1:2.54 |  
| TSAICG | 2,434 | 1,473 | 0.61 | 192.0 kb | 1:1.32 | 4,093 | 2,441 | 0.6 | 181.0 kb | 1:1.34 |  
| UKBB | 54,071 | 25,545 | 0.47 | 180.0 kb | 2.49:1 | 375,800 | 173,198 | 0.46 | 177.0 kb | 2.56:1 |  

The information for this table was collected using `collect_cohort_stats.sh`, and is visualized below using `plot_cnv_stats_per_cohort.R`.  

![Rare CNV stats](https://storage.googleapis.com/rcnv_project/public/rCNV.stats.jpg)  

## Curation Steps: Very rare and ultra-rare CNVs  

In parallel to the creation of the rare CNV dataset (described above), datasets of "very rare" and "ultra-rare" CNVs were also generated.  

These very rare & ultra-rare CNVs were subjected to the same set of filters as above, but imposed increasingly strict maximum frequencies for the following filters:  
 * Does not have substantial overlap with **a CNV of `AF`>0.1% (for very rare CNVs) or `AF`>0.01% (for ultra-rare CNVs)** in any population from WGS resolution in gnomAD-SV (v2.1, non-neuro; see above) ([Collins\*, Brand\*, _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/biorxiv/early/2019/03/14/578674))  
 * Does not have substantial overlap with **a CNV of `AF`>0.1% (for very rare CNVs) or `AF`>0.01% (for ultra-rare CNVs)** in any population from WGS resolution in the 1000 Genomes Project, Phase III ([Sudmant _et al._, _Nature_ (2015)](https://www.nature.com/articles/nature15394))  
 * Does not have substantial overlap with **a CNV of `AF`>0.1% (for very rare CNVs) or `AF`>0.01% (for ultra-rare CNVs)** from WGS resolution in the NIH Center for Common Disease Genetics SV callset ([Abel _et al._, _bioRxiv_ (2018)](https://www.biorxiv.org/content/10.1101/508515v1))  
 * Does not have substantial overlap with other CNVs in **at least 0.1% (for very rare CNVs) or 0.01% (for ultra-rare CNVs) of all samples** within the same dataset or in any of the other array CNV datasets (compared pairwise in serial)


### Very rare CNV callset properties

The properties of each very rare CNV callset are listed below after the above filtering steps.  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| BCH | 3,591 | 2,247 | 0.63 | 363.0 kb | 1:1.40 | 0 | 0 | 0 | NA | NA |  
| CHOP | 153,870 | 49,574 | 0.32 | 225.0 kb | 1:1.22 | 24,161 | 8,666 | 0.36 | 219.0 kb | 1:1.21 |  
| Coe | 29,083 | 12,008 | 0.41 | 344.0 kb | 1:1.17 | 11,256 | 4,739 | 0.42 | 204.0 kb | 1:1.66 |  
| Cooper | 0 | 0 | 0 | NA | NA | 8,329 | 2,440 | 0.29 | 187.0 kb | 1.38:1 |  
| Epi25k | 12,053 | 4,779 | 0.4 | 202.0 kb | 1:1.20 | 8,173 | 2,856 | 0.35 | 205.0 kb | 1:1.25 |  
| GDX | 9,958 | 7,034 | 0.71 | 326.0 kb | 1:1.29 | 0 | 0 | 0 | NA | NA |  
| IU | 1,576 | 975 | 0.62 | 244.0 kb | 1:1.41 | 0 | 0 | 0 | NA | NA |  
| PGC | 21,094 | 7,849 | 0.37 | 222.0 kb | 1:1.50 | 20,277 | 6,988 | 0.34 | 209.0 kb | 1:1.48 |  
| SSC | 2,795 | 1,099 | 0.39 | 236.0 kb | 1:1.33 | 0 | 0 | 0 | NA | NA |  
| SickKids | 2,689 | 2,225 | 0.83 | 198.0 kb | 1:1.31 | 0 | 0 | 0 | NA | NA |  
| TCGA | 0 | 0 | 0 | NA | NA | 8,670 | 2,393 | 0.28 | 195.0 kb | 1:2.17 |  
| TSAICG | 2,434 | 908 | 0.37 | 225.0 kb | 1:1.52 | 4,093 | 1,462 | 0.36 | 216.0 kb | 1:1.53 |  
| UKBB | 54,071 | 16,718 | 0.31 | 195.0 kb | 2.71:1 | 375,800 | 113,675 | 0.3 | 192.0 kb | 2.75:1 |  

The information for this table was collected using `collect_cohort_stats.sh`, and is visualized below using `plot_cnv_stats_per_cohort.R`.  

![Very rare CNV stats](https://storage.googleapis.com/rcnv_project/public/vCNV.stats.jpg)  


### Ultra-rare CNV callset properties

The properties of each ultra-rare CNV callset are listed below after the above filtering steps.  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| BCH | 3,591 | 1,548 | 0.43 | 454.0 kb | 1:1.35 | 0 | 0 | 0 | NA | NA |  
| CHOP | 153,870 | 27,798 | 0.18 | 267.0 kb | 1:1.34 | 24,161 | 4,884 | 0.2 | 254.0 kb | 1:1.36 |  
| Coe | 29,083 | 7,992 | 0.27 | 432.0 kb | 1:1.09 | 11,256 | 2,577 | 0.23 | 242.0 kb | 1:1.88 |  
| Cooper | 0 | 0 | 0 | NA | NA | 8,329 | 1,340 | 0.16 | 208.0 kb | 1.14:1 |  
| Epi25k | 12,053 | 2,783 | 0.23 | 236.0 kb | 1:1.36 | 8,173 | 1,500 | 0.18 | 232.0 kb | 1:1.32 |  
| GDX | 9,958 | 4,589 | 0.46 | 411.0 kb | 1:1.19 | 0 | 0 | 0 | NA | NA |  
| IU | 1,576 | 568 | 0.36 | 315.0 kb | 1:1.41 | 0 | 0 | 0 | NA | NA |  
| PGC | 21,094 | 4,516 | 0.21 | 253.0 kb | 1:1.71 | 20,277 | 3,932 | 0.19 | 239.0 kb | 1:1.61 |  
| SSC | 2,795 | 638 | 0.23 | 268.0 kb | 1:1.32 | 0 | 0 | 0 | NA | NA |  
| SickKids | 2,689 | 1,472 | 0.55 | 224.0 kb | 1:1.54 | 0 | 0 | 0 | NA | NA |  
| TCGA | 0 | 0 | 0 | NA | NA | 8,670 | 1,203 | 0.14 | 208.0 kb | 1:2.02 |  
| TSAICG | 2,434 | 518 | 0.21 | 251.0 kb | 1:1.85 | 4,093 | 816 | 0.2 | 246.0 kb | 1:1.86 |  
| UKBB | 54,071 | 9,904 | 0.18 | 215.0 kb | 2.47:1 | 375,800 | 66,770 | 0.18 | 212.0 kb | 2.46:1 |  

The information for this table was collected using `collect_cohort_stats.sh`, and is visualized below using `plot_cnv_stats_per_cohort.R`.  

![Ultra-rare CNV stats](https://storage.googleapis.com/rcnv_project/public/uCNV.stats.jpg)  

---  

## Case-control "metacohorts"  

For analyses, we combined CNV data from multiple cohorts into four matched groups, dubbed **metacohorts**, to control for technical differences between individual data sources and cohorts.  

Individual cohorts were assigned to metacohorts on the basis of similarity in CNV counts across the genome, which matched expectations based on platforms and processing pipelines for each callset.  

These metacohorts represent the basic unit on which all burden testing was performed, and are described in the table below.  

For completeness, we also performed identical analyses on a pooled dataset of all samples, dubbed the **megacohort**.  

| Metacohort ID | Case Source(s) | Number of Cases | Control Sources(s) | Number of Controls |  
| :--- | :--- | ---: | :--- | ---: |  
| `meta1` | BCH, Coe, GDX, IU | 44,229 | Coe, Cooper | 19,585 |  
| `meta2` | Epi25k, PGC, SSC, SickKids, TSAICG | 41,065 | Epi25k, PGC, TCGA, TSAICG | 41,213 |  
| `meta3` | CHOP | 153,870 | CHOP | 24,161 |  
| `meta4` | UKBB | 54,071 | UKBB | 375,800 |  
| `mega` | BCH, CHOP, Coe, Epi25k, GDX, IU, PGC, SSC, SickKids, TSAICG, UKBB | 293,235 | CHOP, Coe, Cooper, Epi25k, PGC, TCGA, TSAICG, UKBB | 460,759 |  

### Metacohort rare CNV callset properties  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| meta1 | 44,229 | 31,601 | 0.71 | 264.0 kb | 1:1.24 | 19,585 | 11,982 | 0.61 | 173.0 kb | 1:1.23 |  
| meta2 | 41,065 | 27,788 | 0.68 | 183.0 kb | 1:1.23 | 41,213 | 23,561 | 0.57 | 181.0 kb | 1:1.42 |  
| meta3 | 153,870 | 81,530 | 0.53 | 195.0 kb | 1:1.05 | 24,161 | 15,205 | 0.63 | 192.0 kb | 1.01:1 |  
| meta4 | 54,071 | 25,545 | 0.47 | 180.0 kb | 2.49:1 | 375,800 | 173,198 | 0.46 | 177.0 kb | 2.56:1 |  
| mega | 293,235 | 166,464 | 0.57 | 200.0 kb | 1.03:1 | 460,759 | 223,946 | 0.49 | 178.0 kb | 1.92:1 |  

![Rare CNV stats](https://storage.googleapis.com/rcnv_project/public/rCNV.metacohort.stats.jpg)  

### Metacohort very rare CNV callset properties  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| meta1 | 44,229 | 22,264 | 0.5 | 333.0 kb | 1:1.24 | 19,585 | 7,179 | 0.37 | 197.0 kb | 1:1.24 |  
| meta2 | 41,065 | 16,860 | 0.41 | 214.0 kb | 1:1.37 | 41,213 | 13,699 | 0.33 | 206.0 kb | 1:1.53 |  
| meta3 | 153,870 | 49,574 | 0.32 | 225.0 kb | 1:1.22 | 24,161 | 8,666 | 0.36 | 219.0 kb | 1:1.21 |  
| meta4 | 54,071 | 16,718 | 0.31 | 195.0 kb | 2.71:1 | 375,800 | 113,675 | 0.3 | 192.0 kb | 2.75:1 |  
| mega | 293,235 | 105,416 | 0.36 | 233.0 kb | 1:1.04 | 460,759 | 143,219 | 0.31 | 195.0 kb | 2.03:1 |  

![Very rare CNV stats](https://storage.googleapis.com/rcnv_project/public/vCNV.metacohort.stats.jpg)  

### Metacohort ultra-rare CNV callset properties  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| meta1 | 44,229 | 14,697 | 0.33 | 421.0 kb | 1:1.16 | 19,585 | 3,917 | 0.2 | 230.0 kb | 1:1.43 |  
| meta2 | 41,065 | 9,927 | 0.24 | 245.0 kb | 1:1.56 | 41,213 | 7,451 | 0.18 | 233.0 kb | 1:1.63 |  
| meta3 | 153,870 | 27,798 | 0.18 | 267.0 kb | 1:1.34 | 24,161 | 4,884 | 0.2 | 254.0 kb | 1:1.36 |  
| meta4 | 54,071 | 9,904 | 0.18 | 215.0 kb | 2.47:1 | 375,800 | 66,770 | 0.18 | 212.0 kb | 2.46:1 |  
| mega | 293,235 | 62,326 | 0.21 | 278.0 kb | 1:1.10 | 460,759 | 83,022 | 0.18 | 217.0 kb | 1.86:1 |  

![Ultra-rare CNV stats](https://storage.googleapis.com/rcnv_project/public/uCNV.metacohort.stats.jpg)  

The information for these tables was collected using `collect_cohort_stats.sh` and visualized using `plot_cnv_stats_per_cohort.R`.  
