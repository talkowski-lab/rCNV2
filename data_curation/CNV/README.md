# CNV Data Curation

We aggregated existing microarray-based CNV calls from numerous sources. This process is described below.  

All commands executed to filter the CNV data as descrbed below are contained in `filter_CNV_data.sh`.  

In practice, the commands in `filter_CNV_data.sh` were parallelized in [FireCloud/Terra](https://portal.firecloud.org) with `filter_CNV_data.wdl`.  

### CNV data sources  

We aggregated CNV data from multiple sources, listed below:  

| Dataset | Citation | PMID | Platform(s) | Build | Phenos | N Cases | N Ctrls |
| --- | :--- | :--- | :--- | :--- | :--- | ---: | ---: |
| PGC | [Marshall _et al._, _Nat. Genet._ (2017)](https://www.nature.com/articles/ng.3725) | [27869829](https://www.ncbi.nlm.nih.gov/pubmed/27869829) | Affy 6.0 (37%), Omni Express (31%), Omni Express Plus (12%), Other (20%) | hg18 | Schizophrenia | 21,094 | 20,277 |
| Cooper<sup>1</sup> | [Cooper _et al._, _Nat. Genet._ (2011)](https://www.nature.com/articles/ng.909) | [21841781](https://www.ncbi.nlm.nih.gov/pubmed/21841781) | Ill. 550k-610k (75%), Custom 1.2M (25%) | hg19 | Developmental disorders | 0<sup>1</sup> | 8,329 |
| Coe | [Coe _et al._, _Nat. Genet._ (2014)](https://www.nature.com/articles/ng.3092) | [25217958](https://www.ncbi.nlm.nih.gov/pubmed/25217958) | Cases: SignatureChip OS v2.0 (58%), SignatureChip OS v1.0 (34%), Other (8%); Controls: Affy 6.0 (100%) | hg19 | Developmental disorders | 29,083 | 11,256 |
| SSC<sup>2</sup> | [Sanders _et at._, _Neuron_ (2015)](https://www.sciencedirect.com/science/article/pii/S0896627315007734?) | [26402605](https://www.ncbi.nlm.nih.gov/pubmed/26402605) | Omni 1Mv3 (46%), Omni 2.5 (41%), Omni 1Mv1 (13%) | hg18 | Autism | 2,795 | 0<sup>2</sup> |
| UKBB | [Macé _et al._, _Nat. Comms._ (2017)](https://www.nature.com/articles/s41467-017-00556-x) | [28963451](https://www.ncbi.nlm.nih.gov/pubmed/28963451) | UKBB Affy Axiom (100%) | hg19 | Mixed | 54,071<sup>3</sup> | 375,800<sup>3</sup> |
| CHOP | - | - | Mixed Illumina SNP genotyping platforms | hg19 | Mixed | 153,870<sup>3</sup> | 24,161<sup>3</sup> |
| GDX | - | - | aCGH (?) | hg18 & hg19 | Mixed | 9,959 | 0 |
| TSAICG | [Huang _et al._, _Neuron_ (2017)](https://www.sciencedirect.com/science/article/pii/S0896627317305081) | [28641109](https://www.ncbi.nlm.nih.gov/pubmed/28641109) | OmniExpress (100%) | hg19 | Tourette Syndrome | 2,434 | 4,093 |
| BCH | [Talkowski _et al._, _Cell_ (2012)](https://www.sciencedirect.com/science/article/pii/S0092867412004114) | [22521361](https://www.ncbi.nlm.nih.gov/pubmed/22521361) | aCGH (?) | hg18 | Mixed | 3,591 | 0 |  
| TCGA | [Zack _et al._, _Nat. Genet._ (2013)](https://www.nature.com/articles/ng.2760) | [24071852](https://www.ncbi.nlm.nih.gov/pubmed/24071852) | Affy 6.0 (100%) | hg19 | Cancer | 0<sup>4</sup> | 8,670<sup>4</sup> |  
| Epi25k | [Niestroj _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/10.1101/651299v1) | - | Illumina GSA-MD v1.0 (100%) | hg19 | Epilepsy | 12,758<sup>3</sup> | 8,478<sup>3</sup> |  
| SickKids | [Zarrei _et al._, _NPJ Genomic Medicine_ (2019)](https://www.nature.com/articles/s41525-019-0098-3) | [31602316](https://www.ncbi.nlm.nih.gov/pubmed/31602316) | Affy 6.0 (100%) | hg19 | NDDs | 2,691 | 0<sup>2</sup> |  
| IU<sup>5</sup> | - | - | CMA? | hg19? | Mixed | 1,577 | 0 |  

#### Notes on raw CNV data   
1. Only retained control samples from Cooper _et al._. All cases from Cooper _et al._ also appear in Coe _et al._.  
2. Only retained affected children from Sanders _et al._ and Zarrei _et al._, since all controls were first-degree relatives of affected cases.  
3. Counts represent the number of samples retained after filtering outliers, described below.  
4. Only retained normal samples from tumor:normal pairs, and excluded any donors with known blood cancer.  
5. Excluded samples from Indiana University cohort derived from buccal swab DNA, samples with known aneuploidies or large runs of homozygosity, and samples with no phenotypic indication specified.  

### Raw CNV data processing steps  

All CNV data native to hg18 was lifted over to hg19 coordinate space using UCSC liftOver, requiring at least 50% of the original CNV to map to hg19 successfully in order to be retained.  

Some datasets required manual curation prior to inclusion. Where necessary, these steps are enumerated below:  

 * **SSC**: CNVs were filtered on pCNV ≤ 10<sup>-9</sup>, per recommendation of the authors.  
 * **UKBB**: CNVs were filtered on quality score ≥ 17 and CNV size ≥ 25kb. After CNV filtering, samples with >10 CNV calls were excluded as outliers as well as any samples with known malignant cancers or chromosomal disorders (e.g., Down's Syndrome or sex chromosome aneuploidies).  
 * **CHOP**: CNVs were filtered on quality score ≥ 40 and CNV size ≥ 25kb while requiring at least 10 SNPs per CNV. After CNV filtering, samples with `LRR_SD` < 0.25, >20 CNV calls, or SNP call rate < 98% were excluded as outliers, as well as samples genotyped on arrays with < 175k SNP probes or samples labeled as cancer or Down's Syndrome patients.  
 * **TCGA**: CNVs were filtered on ≥ 10 probes and ≥ 25kb. Deletions were required to have a mean log<sub>2</sub> intensity ≤ -1 and duplications were required to have a mean log<sub>2</sub> intensity of ≥ 0.5849625.  
 * **Epi25k**: CNVs were filtered on ≥ 10 probes and ≥ 20kb. After CNV filtering, samples with >25 CNV calls were excluded as outliers.  

### Raw CNV callset properties  

The properties of each callset are listed below after initial data processing steps but prior to further filtering.  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| PGC | 21,094 | 42,096 | 2.00 | 79.7 kb | 1:1.04 | 20,277 | 40,464 | 2.00 | 78.2 kb | 1:1.03 |
| Cooper | 0 | 0 | - | - | - | 8,329 | 432,478 | 51.92 | 1.9 kb | 8.04:1 |
| Coe | 29,083 | 28,782 | 0.99 | 188.4 kb | 1:1.11 | 11,256 | 273,331 | 24.28 | 53.4 kb | 1.24:1 |
| SSC | 2,795 | 30,867 | 11.04 | 21.0 kb | 3.09:1 | 0 | 0 | - | - | - |
| UKBB | 54,071 | 140,200 | 2.59 | 88.2 kb | 5.62:1 | 375,800 | 964,768 | 2.57 | 87.4 kb | 5.68:1 |
| CHOP | 153,870 | 690,521 | 4.49 | 89.0 kb | 1:1.05 | 24,161 | 114,881 | 4.75 | 86.8 kb | 1.71:1 |
| GDX | 9,959 | 20,789 | 2.09 | 196.3 kb | 1:1.76 | 0 | 0 | - | - | - |
| TSAICG | 2,434 | 3,541 | 1.45 | 91.1 kb | 1.01:1 | 4,093 | 5,834 | 1.43 | 91.3 kb | 1:1.08 |
| BCH | 3,591 | 5,211 | 1.45 | 206.4 kb | 1:1.27 | 0 | 0 | - | - | - |
| TCGA | 0 | 0 | - | - | - | 8,670 | 201,712 | 23.27 | 101.1 kb | 1.19:1 |


The information for this table was collected using `collect_cohort_stats.sh`.  

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
gs://rcnv_project/raw_data/cnv/GDX.raw.bed.gz
gs://rcnv_project/raw_data/cnv/GDX.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/PGC.raw.bed.gz
gs://rcnv_project/raw_data/cnv/PGC.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/SSC.raw.bed.gz
gs://rcnv_project/raw_data/cnv/SSC.raw.bed.gz.tbi
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
 * CNV size ≥ 100kb and ≤ 10Mb
 * Does not have substantial overlap<sup>1</sup> with a common (`AF`>1%) CNV in any population from WGS resolution in gnomAD-SV<sup>2</sup> ([Collins\*, Brand\*, _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/biorxiv/early/2019/03/14/578674))  
 * Does not have substantial overlap<sup>1</sup> with a common (`AF`>1%) CNV in any population from WGS resolution in the 1000 Genomes Project, Phase III ([Sudmant _et al._, _Nature_ (2015)](https://www.nature.com/articles/nature15394))  
 * Does not have substantial overlap<sup>1</sup> with a common (`AF`>1%) CNV from WGS resolution in the NIH Center for Common Disease Genetics SV callset ([Abel _et al._, _bioRxiv_ (2018)](https://www.biorxiv.org/content/10.1101/508515v1))  
 * Does not have substantial overlap<sup>1</sup> with other CNVs in at least 1% of all samples
 * Does not have substantial overlap<sup>1</sup> with other CNVs in at least 1% of all samples within the same dataset or in any of the other array CNV datasets (compared pairwise in serial)
 * Not substantially covered<sup>3</sup> by somatically hypermutable sites (as applied in [Collins\*, Brand\*, _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/biorxiv/early/2019/03/14/578674))  
 * Not substantially covered<sup>3</sup> by segmental duplications and/or simple repeats  
 * Not substantially covered<sup>3</sup> by N-masked regions of the hg19 reference genome assembly  

#### Notes on curation  
1. "Substantial" overlap determined based on ≥50% reciprocal overlap using BEDTools ([Quinlan & Hall, _Bioinformatics_ (2010)](https://www.ncbi.nlm.nih.gov/pubmed/20110278)). For overlap-based comparisons, both breakpoints were required to be within ±50kb, and CNV type (DEL vs. DUP) was required to match.  
2. The version of gnomAD-SV used for this analysis (gnomAD-SV v2.1, non-neuro) included 8,342 samples without known neuropsychiatric disorders as available from the gnomAD website](https://gnomad.broadinstitute.org/downloads/) and described in [Collins*, Brand*, et al.](https://www.biorxiv.org/content/10.1101/578674v1)  
3. "Substantial" coverage determined based on ≥30% coverage per BEDTools coverage ([Quinlan & Hall, _Bioinformatics_ (2010)](https://www.ncbi.nlm.nih.gov/pubmed/20110278)).

### Rare CNV callset properties

The properties of each rare CNV callset are listed below after the above filtering steps.  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| PGC | 21,094 | 13,852 | 0.66 | 190.7 kb | 1:1.37 | 20,277 | 12,710 | 0.63 | 184.9 kb | 1:1.33 |
| Cooper | 0 | 0 | - | - | - | 8,329 | 4,318 | 0.52 | 173.2 kb | 1.32:1 |
| Coe | 29,083 | 18,259 | 0.63 | 263.9 kb | 1:1.22 | 11,256 | 10,154 | 0.90 | 184.9 kb | 1:1.81 |
| SSC | 2,795 | 2,053 | 0.73 | 192.3 kb | 1:1.17 | 0 | 0 | - | - | - |
| UKBB | 54,071 | 26,928 | 0.50 | 179.6 kb | 2.37:1 | 375,800 | 183,149 | 0.49 | 177.5 kb | 2.43:1 |
| CHOP | 153,870 | 112,644 | 0.73 | 199.0 kb | 1:1.50 | 24,161 | 17,795 | 0.74 | 188.2 kb | 1:1.14 |
| GDX | 9,959 | 9,820 | 0.99 | 257.5 kb | 1:1.36 | 0 | 0 | - | - | - |
| TSAICG | 2,434 | 1,505 | 0.62 | 194.8 kb | 1:1.37 | 4,093 | 2,509 | 0.61 | 184.6 kb | 1:1.38 |
| BCH | 3,591 | 2,784 | 0.78 | 329.0 kb | 1:1.45 | 0 | 0 | - | - | - |
| TCGA | 0 | 0 | - | - | - | 8,670 | 5,019 | 0.58 | 175.0 kb | 1:2.22 |


The information for this table was collected using `collect_cohort_stats.sh`.  

## Curation Steps: Ultra-rare CNVs  

In parallel to the creation of the rare CNV dataset (described above), a dataset of ultra-rare CNVs was also generated.

These ultra-rare CNVs were subjected to the same set of filters as above, with the following modifications:  
 * Does not have substantial overlap with **a non-ultra-rare (`AF`>0.01%) CNV** in any population from WGS resolution in gnomAD-SV (v2.1, non-neuro; see above) ([Collins\*, Brand\*, _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/biorxiv/early/2019/03/14/578674))  
 * Does not have substantial overlap with **a non-ultra-rare (`AF`>0.01%) CNV** in any population from WGS resolution in the 1000 Genomes Project, Phase III ([Sudmant _et al._, _Nature_ (2015)](https://www.nature.com/articles/nature15394))  
 * Does not have substantial overlap with **a non-ultra-rare (`AF`>0.01%) CNV** from WGS resolution in the NIH Center for Common Disease Genetics SV callset ([Abel _et al._, _bioRxiv_ (2018)](https://www.biorxiv.org/content/10.1101/508515v1))  
 * Does not have substantial overlap with other CNVs in **at least 0.01% of all samples** within the same dataset or in any of the other array CNV datasets (compared pairwise in serial)


### Ultra-rare CNV callset properties

The properties of each ultra-rare CNV callset are listed below after the above filtering steps.  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| PGC | 21,094 | 5,691 | 0.27 | 245.4 kb | 1:1.63 | 20,277 | 4,963 | 0.24 | 235.9 kb | 1:1.51 |
| Cooper | 0 | 0 | - | - | - | 8,329 | 1,773 | 0.21 | 211.0 kb | 1.22:1 |
| Coe | 29,083 | 9,203 | 0.32 | 361.2 kb | 1:1.16 | 11,256 | 3,288 | 0.29 | 238.1 kb | 1:1.78 |
| SSC | 2,795 | 742 | 0.27 | 260.5 kb | 1:1.23 | 0 | 0 | - | - | - |
| UKBB | 54,071 | 12,733 | 0.24 | 216.2 kb | 2.30:1 | 375,800 | 85,726 | 0.23 | 213.7 kb | 2.31:1 |
| CHOP | 153,870 | 41,951 | 0.27 | 254.1 kb | 1:1.48 | 24,161 | 7,502 | 0.31 | 237.6 kb | 1:1.49 |
| GDX | 9,959 | 5,300 | 0.53 | 355.0 kb | 1:1.24 | 0 | 0 | - | - | - |
| TSAICG | 2,434 | 656 | 0.27 | 240.3 kb | 1:1.55 | 4,093 | 1,017 | 0.25 | 238.9 kb | 1:1.69 |
| BCH | 3,591 | 1,750 | 0.49 | 388.5 kb | 1:1.40 | 0 | 0 | - | - | - |
| TCGA | 0 | 0 | - | - | - | 8,670 | 1,604 | 0.19 | 208.1 kb | 1:2.06 |

The information for this table was collected using `collect_cohort_stats.sh`.  

---  

## Case-control "metacohorts"  

To control for potential technical differences between cohorts, we combined CNV data from multiple cohorts into four matched groups for burden testing, dubbed **metacohorts**.  

Individual cohorts were assigned to metacohorts on the basis of similarity in CNV counts across the genome, which matched expectations based on platforms and processing pipelines for each callset.  

These metacohorts represent the basic unit on which all burden testing was performed, and are described in the table below.  

For completeness, we also performed identical analyses on a pooled dataset of all samples, dubbed the **megacohort**.  

| Metacohort ID | Case Source(s) | Number of Cases | Control Sources(s) | Number of Controls |  
| :--- | :--- | ---: | :--- | ---: |  
| `meta1` | Coe, BCH, GDX | 42,653 | Coe, Cooper | 19,585 |  
| `meta2` | PGC, SSC, TSAICG | 26,323 | TCGA, PGC, TSAICG | 33,040 |  
| `meta3` | CHOP | 153,870 | CHOP | 24,161 |  
| `meta4` | UKBB | 54,071 | UKBB | 375,800 |  
| `mega` | Coe, BCH, GDX, PGC, SSC, TSAICG, CHOP, UKBB | 276,917 | Coe, Cooper, TCGA, PGC, TSAICG, CHOP, UKBB | 452,586 |  


### Metacohort rare CNV callset properties  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| meta1 | 42,653 | 30,863 | 0.72 | 265.4 kb | 1:1.28 | 19,585 | 14,472 | 0.74 | 182.0 kb | 1:1.38 |
| meta2 | 26,323 | 17,410 | 0.66 | 191.4 kb | 1:1.34 | 33,040 | 20,238 | 0.61 | 183.2 kb | 1:1.51 |
| meta3 | 153,870 | 112,644 | 0.73 | 199.0 kb | 1:1.50 | 24,161 | 17,795 | 0.74 | 188.2 kb | 1:1.14 |
| meta4 | 54,071 | 26,928 | 0.50 | 179.6 kb | 2.37:1 | 375,800 | 183,149 | 0.49 | 177.5 kb | 2.43:1 |
| mega | 276,917 | 187,845 | 0.68 | 203.8 kb | 1:1.21 | 452,586 | 235,654 | 0.52 | 179.0 kb | 1.82:1 |

### Metacohort ultra-rare CNV callset properties  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| meta1 | 42,653 | 16,253 | 0.38 | 360.6 kb | 1:1.21 | 19,585 | 5,061 | 0.26 | 228.1 kb | 1:1.35 |
| meta2 | 26,323 | 7,089 | 0.27 | 246.5 kb | 1:1.57 | 33,040 | 7,584 | 0.23 | 229.4 kb | 1:1.63 |
| meta3 | 153,870 | 41,951 | 0.27 | 254.1 kb | 1:1.48 | 24,161 | 7,502 | 0.31 | 237.6 kb | 1:1.49 |
| meta4 | 54,071 | 12,733 | 0.24 | 216.2 kb | 2.30:1 | 375,800 | 85,726 | 0.23 | 213.7 kb | 2.31:1 |
| mega | 276,917 | 78,026 | 0.28 | 262.0 kb | 1:1.17 | 452,586 | 105,873 | 0.23 | 217.0 kb | 1.79:1 |

The information for these tables was collected using `collect_cohort_stats.sh`.  
