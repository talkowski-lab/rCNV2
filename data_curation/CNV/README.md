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

 * **CHOP**: CNVs were filtered on quality score ≥40 and CNV size ≥25kb while requiring at least 10 SNPs per CNV. After CNV filtering, samples with `LRR_SD` <0.25, >20 CNV calls, or SNP call rate <98% were excluded as outliers, as well as samples genotyped on arrays with <175k SNP probes or samples labeled as cancer or Down's Syndrome patients.  
 * **Epi25k**: CNVs were filtered on ≥10 probes and ≥25kb. After CNV filtering, samples with >25 CNV calls were excluded as outliers.  
 * **SSC**: CNVs were filtered on pCNV ≤10<sup>-9</sup>, per recommendation of the authors.  
 * **SickKids**: CNVs were filtered on ≥25kb. After CNV filtering, samples with >80 CNV calls were excluded as outliers.
 * **TCGA**: CNVs were filtered on ≥10 probes and ≥25kb. Deletions were required to have a mean log<sub>2</sub> intensity ≤-1 and duplications were required to have a mean log<sub>2</sub> intensity of ≥0.5849625. After filtering, CNVs per sample were defragmented using `defragment_cnvs.py` with `--max-dist 0.25`.  
 * **UKBB**: CNVs were filtered on quality score ≥17 and CNV size ≥25kb. After CNV filtering, samples with >10 CNV calls were excluded as outliers as well as any samples with known malignant cancers or chromosomal disorders (e.g., Down's Syndrome or sex chromosome aneuploidies).  

### Raw CNV callset properties  

The properties of each callset are listed below after initial data processing steps but prior to further filtering.  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| BCH | 3,591 | 5,211 | 1.45 | 206.0 kb | 1:1.27 | 0 | 0 | 0 | NA | NA |  
| CHOP | 153,870 | 690,521 | 4.49 | 89.0 kb | 1:1.05 | 24,161 | 114,881 | 4.75 | 87.0 kb | 1.71:1 |  
| Coe | 29,083 | 28,782 | 0.99 | 188.0 kb | 1:1.11 | 11,256 | 273,331 | 24.28 | 53.0 kb | 1.24:1 |  
| Cooper | 0 | 0 | 0 | NA | NA | 8,329 | 432,478 | 51.92 | 2.0 kb | 8.04:1 |  
| Epi25k | 12,053 | 94,172 | 7.81 | 97.0 kb | 1:1.56 | 8,173 | 65,329 | 7.99 | 100.0 kb | 1:1.35 |  
| GDX | 9,958 | 20,789 | 2.09 | 196.0 kb | 1:1.76 | 0 | 0 | 0 | NA | NA |  
| IU | 1,576 | 4,940 | 3.13 | 99.0 kb | 1.38:1 | 0 | 0 | 0 | NA | NA |  
| PGC | 21,094 | 42,096 | 2.0 | 80.0 kb | 1:1.04 | 20,277 | 40,464 | 2.0 | 78.0 kb | 1:1.03 |  
| SSC | 2,795 | 30,867 | 11.04 | 21.0 kb | 3.09:1 | 0 | 0 | 0 | NA | NA |  
| SickKids | 2,689 | 98,970 | 36.81 | 52.0 kb | 1:1.02 | 0 | 0 | 0 | NA | NA |  
| TCGA | 0 | 0 | 0 | NA | NA | 8,670 | 201,712 | 23.27 | 101.0 kb | 1.19:1 |  
| TSAICG | 2,434 | 3,541 | 1.45 | 91.0 kb | 1.01:1 | 4,093 | 5,834 | 1.43 | 91.0 kb | 1:1.08 |  
| UKBB | 54,071 | 140,200 | 2.59 | 88.0 kb | 5.62:1 | 375,800 | 964,768 | 2.57 | 87.0 kb | 5.68:1 |  

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
 * CNV size ≥ 100kb and ≤ 10Mb
 * Does not have substantial overlap<sup>1</sup> with a common CNV (`AF`>1%) in any population from WGS resolution in gnomAD-SV<sup>2</sup> ([Collins\*, Brand\*, _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/biorxiv/early/2019/03/14/578674))  
 * Does not have substantial overlap<sup>1</sup> with a common CNV (`AF`>1%) in any population from WGS resolution in the 1000 Genomes Project, Phase III ([Sudmant _et al._, _Nature_ (2015)](https://www.nature.com/articles/nature15394))  
 * Does not have substantial overlap<sup>1</sup> with a common CNV (`AF`>1%) from WGS resolution in the NIH Center for Common Disease Genetics SV callset ([Abel _et al._, _bioRxiv_ (2018)](https://www.biorxiv.org/content/10.1101/508515v1))  
 * Does not have substantial overlap<sup>1</sup> with other CNVs in at least 1% of all samples
 * Does not have substantial overlap<sup>1</sup> with other CNVs in at least 1% of all samples within the same dataset or in any of the other array CNV datasets (compared pairwise in serial)
 * Not substantially covered<sup>3</sup> by somatically hypermutable sites (as applied in [Collins\*, Brand\*, _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/biorxiv/early/2019/03/14/578674))  
 * Not substantially covered<sup>3</sup> by segmental duplications and/or simple, low-complexity, or satellite repeats  
 * Not substantially covered<sup>3</sup> by N-masked regions of the hg19 reference genome assembly  

#### Notes on curation  
1. "Substantial" overlap determined based on ≥50% reciprocal overlap using BEDTools ([Quinlan & Hall, _Bioinformatics_ (2010)](https://www.ncbi.nlm.nih.gov/pubmed/20110278)). For overlap-based comparisons, both breakpoints were required to be within ±50kb, and CNV type (DEL vs. DUP) was required to match.  
2. The version of gnomAD-SV used for this analysis (gnomAD-SV v2.1, non-neuro) included 8,342 samples without known neuropsychiatric disorders as [available from the gnomAD website](https://gnomad.broadinstitute.org/downloads/) and described in [Collins\*, Brand\*, _et al._](https://www.biorxiv.org/content/10.1101/578674v1)  
3. "Substantial" coverage determined based on ≥30% coverage per BEDTools coverage ([Quinlan & Hall, _Bioinformatics_ (2010)](https://www.ncbi.nlm.nih.gov/pubmed/20110278)).

### Rare CNV callset properties

The properties of each rare CNV callset are listed below after the above filtering steps.  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| BCH | 3,591 | 2,745 | 0.76 | 327.0 kb | 1:1.45 | 0 | 0 | 0 | NA | NA |  
| CHOP | 153,870 | 101,279 | 0.66 | 200.0 kb | 1:1.41 | 24,161 | 16,937 | 0.7 | 185.0 kb | 1:1.09 |  
| Coe | 29,083 | 17,601 | 0.61 | 263.0 kb | 1:1.20 | 11,256 | 8,196 | 0.73 | 177.0 kb | 1:1.55 |  
| Cooper | 0 | 0 | 0 | NA | NA | 8,329 | 3,909 | 0.47 | 173.0 kb | 1.44:1 |  
| Epi25k | 12,053 | 8,034 | 0.67 | 177.0 kb | 1:1.08 | 8,173 | 5,054 | 0.62 | 186.0 kb | 1:1.16 |  
| GDX | 9,958 | 9,617 | 0.97 | 256.0 kb | 1:1.35 | 0 | 0 | 0 | NA | NA |  
| IU | 1,576 | 1,483 | 0.94 | 208.0 kb | 1:1.37 | 0 | 0 | 0 | NA | NA |  
| PGC | 21,094 | 13,125 | 0.62 | 188.0 kb | 1:1.35 | 20,277 | 11,943 | 0.59 | 181.0 kb | 1:1.32 |  
| SSC | 2,795 | 1,926 | 0.69 | 194.0 kb | 1:1.12 | 0 | 0 | 0 | NA | NA |  
| SickKids | 2,689 | 4,652 | 1.73 | 185.0 kb | 1:1.34 | 0 | 0 | 0 | NA | NA |  
| TCGA | 0 | 0 | 0 | NA | NA | 8,670 | 4,543 | 0.52 | 179.0 kb | 1:2.52 |  
| TSAICG | 2,434 | 1,473 | 0.61 | 195.0 kb | 1:1.39 | 4,093 | 2,464 | 0.6 | 184.0 kb | 1:1.39 |  
| UKBB | 54,071 | 26,283 | 0.49 | 179.0 kb | 2.48:1 | 375,800 | 178,308 | 0.47 | 177.0 kb | 2.54:1 |  

The information for this table was collected using `collect_cohort_stats.sh`, and is visualized below using `plot_cnv_stats_per_cohort.R`.  

![Rare CNV stats](https://storage.googleapis.com/rcnv_project/public/rare_cnv.stats.jpg)  

## Curation Steps: Ultra-rare CNVs  

In parallel to the creation of the rare CNV dataset (described above), a dataset of ultra-rare CNVs was also generated.

These ultra-rare CNVs were subjected to the same set of filters as above, with the following modifications:  
 * Does not have substantial overlap with **a non-ultra-rare CNV (`AF`>0.01%)** in any population from WGS resolution in gnomAD-SV (v2.1, non-neuro; see above) ([Collins\*, Brand\*, _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/biorxiv/early/2019/03/14/578674))  
 * Does not have substantial overlap with **a non-ultra-rare CNV (`AF`>0.01%)** in any population from WGS resolution in the 1000 Genomes Project, Phase III ([Sudmant _et al._, _Nature_ (2015)](https://www.nature.com/articles/nature15394))  
 * Does not have substantial overlap with **a non-ultra-rare CNV (`AF`>0.01%)** from WGS resolution in the NIH Center for Common Disease Genetics SV callset ([Abel _et al._, _bioRxiv_ (2018)](https://www.biorxiv.org/content/10.1101/508515v1))  
 * Does not have substantial overlap with other CNVs in **at least 0.01% of all samples** within the same dataset or in any of the other array CNV datasets (compared pairwise in serial)


### Ultra-rare CNV callset properties

The properties of each ultra-rare CNV callset are listed below after the above filtering steps.  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| BCH | 3,591 | 1,721 | 0.48 | 393.0 kb | 1:1.40 | 0 | 0 | 0 | NA | NA |  
| CHOP | 153,870 | 40,703 | 0.26 | 253.0 kb | 1:1.51 | 24,161 | 7,289 | 0.3 | 238.0 kb | 1:1.53 |  
| Coe | 29,083 | 9,045 | 0.31 | 362.0 kb | 1:1.17 | 11,256 | 3,116 | 0.28 | 235.0 kb | 1:1.80 |  
| Cooper | 0 | 0 | 0 | NA | NA | 8,329 | 1,698 | 0.2 | 210.0 kb | 1.20:1 |  
| Epi25k | 12,053 | 3,469 | 0.29 | 227.0 kb | 1:1.36 | 8,173 | 1,930 | 0.24 | 227.0 kb | 1:1.34 |  
| GDX | 9,958 | 5,228 | 0.53 | 355.0 kb | 1:1.26 | 0 | 0 | 0 | NA | NA |  
| IU | 1,576 | 674 | 0.43 | 287.0 kb | 1:1.49 | 0 | 0 | 0 | NA | NA |  
| PGC | 21,094 | 5,499 | 0.26 | 244.0 kb | 1:1.68 | 20,277 | 4,788 | 0.24 | 233.0 kb | 1:1.57 |  
| SSC | 2,795 | 769 | 0.28 | 259.0 kb | 1:1.34 | 0 | 0 | 0 | NA | NA |  
| SickKids | 2,689 | 2,114 | 0.79 | 224.0 kb | 1:1.55 | 0 | 0 | 0 | NA | NA |  
| TCGA | 0 | 0 | 0 | NA | NA | 8,670 | 1,568 | 0.18 | 207.0 kb | 1:2.10 |  
| TSAICG | 2,434 | 637 | 0.26 | 242.0 kb | 1:1.67 | 4,093 | 996 | 0.24 | 243.0 kb | 1:1.74 |  
| UKBB | 54,071 | 12,457 | 0.23 | 216.0 kb | 2.27:1 | 375,800 | 83,675 | 0.22 | 213.0 kb | 2.28:1 |  

The information for this table was collected using `collect_cohort_stats.sh`, and is visualized below using `plot_cnv_stats_per_cohort.R`.  

![Ultra-rare CNV stats](https://storage.googleapis.com/rcnv_project/public/ultrarare_cnv.stats.jpg)  

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
| meta1 | 44,229 | 31,446 | 0.71 | 263.0 kb | 1:1.27 | 19,585 | 12,105 | 0.62 | 175.0 kb | 1:1.19 |  
| meta2 | 41,065 | 29,210 | 0.71 | 186.0 kb | 1:1.25 | 41,213 | 24,004 | 0.58 | 182.0 kb | 1:1.45 |  
| meta3 | 153,870 | 101,279 | 0.66 | 200.0 kb | 1:1.41 | 24,161 | 16,937 | 0.7 | 185.0 kb | 1:1.09 |  
| meta4 | 54,071 | 26,283 | 0.49 | 179.0 kb | 2.48:1 | 375,800 | 178,308 | 0.47 | 177.0 kb | 2.54:1 |  
| mega | 293,235 | 188,218 | 0.64 | 202.0 kb | 1:1.15 | 460,759 | 231,354 | 0.5 | 178.0 kb | 1.89:1 |  

![Rare CNV stats](https://storage.googleapis.com/rcnv_project/public/rare_cnv.metacohort.stats.jpg)  

### Metacohort ultra-rare CNV callset properties  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| meta1 | 44,229 | 16,668 | 0.38 | 359.0 kb | 1:1.23 | 19,585 | 4,814 | 0.25 | 225.0 kb | 1:1.36 |  
| meta2 | 41,065 | 12,488 | 0.3 | 236.0 kb | 1:1.54 | 41,213 | 9,282 | 0.23 | 228.0 kb | 1:1.61 |  
| meta3 | 153,870 | 40,703 | 0.26 | 253.0 kb | 1:1.51 | 24,161 | 7,289 | 0.3 | 238.0 kb | 1:1.53 |  
| meta4 | 54,071 | 12,457 | 0.23 | 216.0 kb | 2.27:1 | 375,800 | 83,675 | 0.22 | 213.0 kb | 2.28:1 |  
| mega | 293,235 | 82,316 | 0.28 | 259.0 kb | 1:1.21 | 460,759 | 105,060 | 0.23 | 217.0 kb | 1.74:1 |  

![Ultra-rare CNV stats](https://storage.googleapis.com/rcnv_project/public/ultrarare_cnv.metacohort.stats.jpg)  

The information for these tables was collected using `collect_cohort_stats.sh` and visualized using `plot_cnv_stats_per_cohort.R`.  
