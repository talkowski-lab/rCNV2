# CNV Data Curation

We aggregated existing microarray-based CNV calls from numerous sources. This process is described below.  

All commands executed to filter the CNV data as descrbed below are contained in `filter_CNV_data.sh`.  

In practice, the commands in `filter_CNV_data.sh` were parallelized in [FireCloud/Terra](https://portal.firecloud.org) with `filter_CNV_data.wdl`.  

## CNV data sources  

We aggregated CNV data from multiple sources, listed alphabetically below:  

| Dataset | Citation | PMID | Platform(s) | Build | Phenos | N Cases | N Ctrls |
| --- | :--- | :--- | :--- | :--- | :--- | ---: | ---: |
| BCH | [Talkowski _et al._, _Cell_ (2012)](https://www.sciencedirect.com/science/article/pii/S0092867412004114) | [22521361](https://www.ncbi.nlm.nih.gov/pubmed/22521361) | aCGH | hg18 | Mixed | 3,591 | 0 |  
| BioVU | [Roden _et al._, _Clin. Pharmacol. Ther._ (2008)](https://pubmed.ncbi.nlm.nih.gov/18500243/) | [18500243](https://pubmed.ncbi.nlm.nih.gov/18500243/) | Illumina MEGAEx (100%) | hg19 | Mixed | 32,306<sup>1</sup> | 14,661<sup>1</sup> |  
| CHOP | - | - | Mixed Illumina SNP genotyping platforms | hg19 | Mixed | 153,870<sup>1</sup> | 24,161<sup>1</sup> |
| Coe | [Coe _et al._, _Nat. Genet._ (2014)](https://www.nature.com/articles/ng.3092) | [25217958](https://www.ncbi.nlm.nih.gov/pubmed/25217958) | Cases: SignatureChip OS v2.0 (58%), SignatureChip OS v1.0 (34%), Other (8%); Controls: Affy 6.0 (100%) | hg19 | Developmental disorders | 29,083 | 11,256 |
| Cooper<sup>2</sup> | [Cooper _et al._, _Nat. Genet._ (2011)](https://www.nature.com/articles/ng.909) | [21841781](https://www.ncbi.nlm.nih.gov/pubmed/21841781) | Ill. 550k-610k (75%), Custom 1.2M (25%) | hg19 | Developmental disorders | 0<sup>2</sup> | 8,329 |
| Epi25k | [Niestroj _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/10.1101/651299v1) | - | Illumina GSA-MD v1.0 (100%) | hg19 | Epilepsy | 12,053<sup>1</sup> | 8,173<sup>1</sup> |  
| EstBB | [Leitsalu _et al._, _Int. J. Epidemiol._ (2014)](https://academic.oup.com/ije/article/44/4/1137/666872) | [24518929](https://pubmed.ncbi.nlm.nih.gov/24518929/) | Illumina GSA (100%) | hg19 | Mixed | 63,183<sup>1</sup> | 15,659<sup>1</sup> |  
| GDX | - | - | Affy CytoScan HD & Custom Agilent SNP arrays | hg18 & hg19 | Mixed | 74,028<sup>1</sup> | 0 |
| IU<sup>3</sup> | - | - | CMA (?) | hg19 | Mixed | 1,576<sup>3</sup> | 0 |  
| PGC | [Marshall _et al._, _Nat. Genet._ (2017)](https://www.nature.com/articles/ng.3725) | [27869829](https://www.ncbi.nlm.nih.gov/pubmed/27869829) | Affy 6.0 (37%), Omni Express (31%), Omni Express Plus (12%), Other (20%) | hg18 | Schizophrenia | 21,094 | 20,277 |
| SSC<sup>4</sup> | [Sanders _et al._, _Neuron_ (2015)](https://www.sciencedirect.com/science/article/pii/S0896627315007734?) | [26402605](https://www.ncbi.nlm.nih.gov/pubmed/26402605) | Omni 1Mv3 (46%), Omni 2.5 (41%), Omni 1Mv1 (13%) | hg18 | Autism | 2,795 | 0<sup>4</sup> |
| SickKids | [Zarrei _et al._, _NPJ Genomic Medicine_ (2019)](https://www.nature.com/articles/s41525-019-0098-3) | [31602316](https://www.ncbi.nlm.nih.gov/pubmed/31602316) | Affy 6.0 (100%) | hg19 | NDDs | 2,689<sup>1</sup> | 0<sup>4</sup> |  
| TCGA | [Zack _et al._, _Nat. Genet._ (2013)](https://www.nature.com/articles/ng.2760) | [24071852](https://www.ncbi.nlm.nih.gov/pubmed/24071852) | Affy 6.0 (100%) | hg19 | Cancer | 0<sup>5</sup> | 8,670<sup>5</sup> |  
| TSAICG | [Huang _et al._, _Neuron_ (2017)](https://www.sciencedirect.com/science/article/pii/S0896627317305081) | [28641109](https://www.ncbi.nlm.nih.gov/pubmed/28641109) | OmniExpress (100%) | hg19 | Tourette Syndrome | 2,434 | 4,093 |
| UKBB | [Macé _et al._, _Nat. Commun._ (2017)](https://www.nature.com/articles/s41467-017-00556-x) | [28963451](https://www.ncbi.nlm.nih.gov/pubmed/28963451) | UKBB Affy Axiom (100%) | hg19 | Mixed | 54,071<sup>1</sup> | 375,800<sup>1</sup> |

#### Notes on raw CNV data   
1. Counts represent the number of samples retained after filtering outliers, described below.  
2. Only retained control samples from Cooper _et al._ All cases from Cooper _et al._ also appear in Coe _et al._  
3. Excluded samples from Indiana University (IU) cohort derived from buccal swab DNA, samples with known aneuploidies or large runs of homozygosity, and samples with no phenotypic indication specified.  
4. Only retained affected children from Sanders _et al._ and Zarrei _et al._, since all controls were first-degree relatives of affected cases.  
5. Only retained normal samples from tumor:normal pairs, and excluded any donors with known blood cancer.  

## Raw CNV data processing steps  

All CNV data native to hg18 was lifted over to hg19 using UCSC liftOver, requiring at least 50% of the original CNV to map successfully to hg19 in order to be retained.  

Some datasets required manual curation prior to inclusion. Where necessary, these steps are enumerated below:  

 * **BioVU**: TBD [TODO: ADD TEXT HERE]
 * **CHOP**: CNVs were filtered on quality score ≥40 and CNV size ≥25kb while requiring at least 10 SNPs per CNV. After CNV filtering, samples with `LRR_SD` <0.25, >20 CNV calls, or SNP call rate <98% were excluded as outliers, as well as samples genotyped on arrays with <175k SNP probes or samples labeled as cancer or Down's Syndrome patients. Finally, we identified 19 loci with apparently platform-specific artifactual CNV pileups. CNVs covered ≥10% by any of these artifact regions were removed from the callset. Lists of CHOP-specific blacklisted loci for deletions and duplications are provided as [reference files](https://github.com/talkowski-lab/rCNV2/tree/master/refs).  
 * **Epi25k**: CNVs were filtered on ≥10 probes and ≥25kb. After CNV filtering, samples with >25 CNV calls were excluded as outliers.  
 * **EstBB**: Samples were excluded if were not included in SNP imputation, had genotype calls missing for ≥2% of sites, or belonged to two genotyping batches based on visual inspection of genotyping intensity parameters, followed by further exclusion of genotyping plates (≤24 samples per plate) that contained >3 samples with >200 CNV calls. We retained unrelated samples that had been linked to Estonian health registries and that had ≤50 raw CNV calls. We included CNVs with a quality score ≥15, were covered by ≥10 probes, and were ≥25kb in size. Finally, we pruned related samples and any samples with known malignant cancers or chromosomal disorders (e.g., Down's Syndrome or sex chromosome aneuploidies).
 * **GDX**: All CNVs were required to be ≥20kb and <40Mb in length. Except for the minority of 9,958 samples for which additional CNV call metadata was unavailable, all CNVs were further required to not have been annotated as a suspected false positive or mosaic event, have estimated copy numbers ≤1.5 for deletions or ≥2.5 for duplications, include ≥10 probes and have P(CNV) ≤ 10<sup>-10</sup>. Following CNV call filtering, we excluded all samples that either had >10 calls each, were identified as potential biological replicates, were referred for testing due to being a relative of a known carrier of a medically relevant CNV, or had “advanced maternal age” as their indication for testing.
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

Five studies were unable to be defragmented due to inadequate sample-level information: `Coe` (_controls only_), `Cooper`, `PGC`, `BioVU`, and `TSAICG`.

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
| EstBB | 63,183 | 226,083 | 3.58 | 83.0 kb | 1:1.45 | 15,659 | 56,267 | 3.59 | 83.0 kb | 1:1.46 |  
| BioVU | 32,306 | 58,955 | 1.82 | 133.0 kb | 1:1.32 | 14,661 | 27,411 | 1.87 | 132.0 kb | 1:1.40 |  

The information for this table was collected using `collect_cohort_stats.sh`, and is visualized below using `plot_cnv_stats_per_cohort.R`.  

![Raw CNV stats](https://storage.googleapis.com/rcnv_project/public/raw_cnv.stats.jpg)  

### Raw data access  

All raw CNV data files and their tabix indexes are stored in a protected Google Cloud bucket, here:  
```
$ gsutil ls gs://rcnv_project/raw_data/cnv/

gs://rcnv_project/raw_data/cnv/BCH.raw.bed.gz
gs://rcnv_project/raw_data/cnv/BCH.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/BioVU.raw.bed.gz
gs://rcnv_project/raw_data/cnv/BioVU.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/CHOP.raw.bed.gz
gs://rcnv_project/raw_data/cnv/CHOP.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/Coe.raw.bed.gz
gs://rcnv_project/raw_data/cnv/Coe.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/Cooper.raw.bed.gz
gs://rcnv_project/raw_data/cnv/Cooper.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/Epi25k.raw.bed.gz
gs://rcnv_project/raw_data/cnv/Epi25k.raw.bed.gz.tbi
gs://rcnv_project/raw_data/cnv/EstBB.raw.bed.gz
gs://rcnv_project/raw_data/cnv/EstBB.raw.bed.gz.tbi
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
 * Does not have substantial overlap<sup>1</sup> with a common CNV (lower bound of 95% binomial confidence interval of`AF` > 1%) in any population from WGS resolution in gnomAD-SV<sup>2</sup> v2.1 ([Collins\*, Brand\*, _et al._, _Nature_ (2020)](https://www.nature.com/articles/s41586-020-2287-8))  
 * Does not have substantial overlap<sup>1</sup> with a common CNV (lower bound of 95% binomial confidence interval of`AF` > 1%) in any population from WGS resolution in recent high-coverage resequencing of the 1000 Genomes Project, Phase III ([Byrska-Bishop _et al._, _bioRxiv_ (2021)](https://www.biorxiv.org/content/10.1101/2021.02.06.430068v1))  
 * Does not have substantial overlap<sup>1</sup> with a common CNV (lower bound of 95% binomial confidence interval of`AF` > 1%) in any population from WGS resolution in recent high-coverage resequencing of the Human Genome Diversity Panel ([Almarri _et al._, _Cell_ (2020)](https://pubmed.ncbi.nlm.nih.gov/32531199/))  
 * Does not have substantial overlap<sup>1</sup> with a common CNV (lower bound of 95% binomial confidence interval of`AF` > 1%) from WGS resolution in the NIH Center for Common Disease Genetics SV callset ([Abel _et al._, _Nature_ (2020)](https://www.nature.com/articles/s41586-020-2371-0))  
 * Does not have substantial overlap<sup>1</sup> with other CNVs in at least 1% of all samples within the same dataset or in any of the other array CNV datasets (compared pairwise in serial; 1% cutoff scaled adaptively to each cohort corresponding to the upper bound of the 95% binomial confidence interval of 1% frequency according to that cohort's total sample size)<sup>3</sup>
 * Not substantially covered<sup>4</sup> by somatically hypermutable sites (as applied in [Collins\*, Brand\*, _et al._, _Nature_ (2020)](https://www.nature.com/articles/s41586-020-2287-8))<sup>3</sup>  
 * Not substantially covered<sup>4</sup> by segmental duplications and/or simple, low-complexity, or satellite repeats<sup>3</sup>  
 * Not substantially covered<sup>4</sup> by N-masked regions of the hg19 reference genome assembly<sup>3</sup>  

#### Notes on curation  
1. "Substantial" overlap determined based on ≥50% reciprocal overlap using BEDTools ([Quinlan & Hall, _Bioinformatics_ (2010)](https://www.ncbi.nlm.nih.gov/pubmed/20110278)). For overlap-based comparisons, both breakpoints were required to be within ±100kb, and CNV type (DEL vs. DUP) was required to match.  
2. The version of gnomAD-SV used for this analysis (gnomAD-SV v2.1, non-neuro) included 8,342 samples without known neuropsychiatric disorders as [available from the gnomAD website](https://gnomad.broadinstitute.org/downloads/) and described in [Collins\*, Brand\*, _et al._, _Nature_ (2020)](https://www.nature.com/articles/s41586-020-2287-8)  
3. CNVs with ≥75% reciprocal overlap versus [known genomic disorder CNV loci](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/other#genomic-disorders) reported by at least two of six established resources were exempted from this filter step.  
4. "Substantial" coverage determined based on ≥50% coverage per BEDTools coverage ([Quinlan & Hall, _Bioinformatics_ (2010)](https://www.ncbi.nlm.nih.gov/pubmed/20110278)).

### Rare CNV callset properties

The properties of each rare CNV callset are listed below after the above filtering steps.  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| BCH | 3,591 | 3,028 | 0.84 | 342.0 kb | 1:1.40 | 0 | 0 | 0 | NA | NA |  
| CHOP | 153,870 | 91,922 | 0.6 | 194.0 kb | 1:1.06 | 24,161 | 17,104 | 0.71 | 189.0 kb | 1.02:1 |  
| Coe | 29,104 | 19,650 | 0.68 | 264.0 kb | 1:1.18 | 11,256 | 9,364 | 0.83 | 177.0 kb | 1:1.56 |  
| Cooper | 0 | 0 | 0 | NA | NA | 8,329 | 4,186 | 0.5 | 170.0 kb | 1.39:1 |  
| Epi25k | 12,053 | 9,024 | 0.75 | 176.0 kb | 1:1.06 | 8,173 | 5,638 | 0.69 | 185.0 kb | 1:1.13 |  
| GDX | 74,028 | 45,066 | 0.61 | 283.0 kb | 1:1.47 | 0 | 0 | 0 | NA | NA |  
| IU | 1,577 | 1,652 | 1.05 | 209.0 kb | 1:1.34 | 0 | 0 | 0 | NA | NA |  
| PGC | 21,094 | 14,661 | 0.7 | 190.0 kb | 1:1.34 | 20,277 | 13,273 | 0.65 | 182.0 kb | 1:1.33 |  
| SSC | 2,795 | 2,151 | 0.77 | 195.0 kb | 1:1.15 | 0 | 0 | 0 | NA | NA |  
| SickKids | 2,689 | 3,903 | 1.45 | 182.0 kb | 1:1.33 | 0 | 0 | 0 | NA | NA |  
| TCGA | 0 | 0 | 0 | NA | NA | 8,670 | 4,918 | 0.57 | 182.0 kb | 1:2.56 |  
| TSAICG | 2,434 | 1,607 | 0.66 | 198.0 kb | 1:1.31 | 4,093 | 2,663 | 0.65 | 187.0 kb | 1:1.36 |  
| UKBB | 54,071 | 28,845 | 0.53 | 179.0 kb | 2.52:1 | 375,800 | 195,245 | 0.52 | 177.0 kb | 2.59:1 |  
| EstBB | 63,183 | 32,992 | 0.52 | 185.0 kb | 1:1.13 | 15,659 | 8,228 | 0.53 | 181.0 kb | 1:1.12 |  
| BioVU | 32,306 | 27,353 | 0.85 | 179.0 kb | 1:1.21 | 14,661 | 12,176 | 0.83 | 175.0 kb | 1:1.23 |  


The information for this table was collected using `collect_cohort_stats.sh`, and is visualized below using `plot_cnv_stats_per_cohort.R`.  

![Rare CNV stats](https://storage.googleapis.com/rcnv_project/public/rCNV.stats.jpg)  

---  

## Case-control "metacohorts"  

For analyses, we combined CNV data from multiple cohorts into seven matched groups, dubbed **metacohorts**, to control for technical differences between individual data sources and cohorts.  

Individual cohorts were assigned to metacohorts on the basis of similarity in CNV counts across the genome, which matched expectations based on platforms and processing pipelines for each callset.  

These metacohorts represent the basic unit on which all burden testing was performed, and are described in the table below.  

For completeness, we also performed identical analyses on a pooled dataset of all samples, dubbed the **megacohort**.  

| Metacohort ID | Case Source(s) | Number of Cases | Control Sources(s) | Number of Controls |  
| :--- | :--- | ---: | :--- | ---: |  
| `meta1` | Coe | 29,104 | Coe, Cooper | 19,585 |  
| `meta2` | BCH, GDX, IU | 79,196 | GDX, TCGA | 8,670 |  
| `meta3` | Epi25k, PGC, SSC, SickKids, TSAICG | 41,065 | Epi25k, PGC, TSAICG | 32,543 |  
| `meta4` | CHOP | 153,870 | CHOP | 24,161 |  
| `meta5` | UKBB | 54,071 | UKBB | 375,800 |  
| `meta6` | EstBB | 63,183 | EstBB | 15,659 |  
| `meta7` | BioVU | 32,306 | BioVU | 14,661 |  
| `mega` | BCH, CHOP, Coe, Epi25k, GDX, IU, PGC, SSC, SickKids, TSAICG, UKBB, EstBB, BioVU | 452,795 | CHOP, Coe, Cooper, Epi25k, PGC, TCGA, TSAICG, UKBB, EstBB, BioVU | 491,079 |  

### Metacohort rare CNV callset properties  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| meta1 | 29,104 | 19,650 | 0.68 | 264.0 kb | 1:1.18 | 19,585 | 13,550 | 0.69 | 174.0 kb | 1:1.22 |  
| meta2 | 79,196 | 49,746 | 0.63 | 283.0 kb | 1:1.46 | 8,670 | 4,918 | 0.57 | 182.0 kb | 1:2.56 |  
| meta3 | 41,065 | 31,346 | 0.76 | 186.0 kb | 1:1.24 | 32,543 | 21,574 | 0.66 | 183.0 kb | 1:1.28 |  
| meta4 | 153,870 | 91,922 | 0.6 | 194.0 kb | 1:1.06 | 24,161 | 17,104 | 0.71 | 189.0 kb | 1.02:1 |  
| meta5 | 54,071 | 28,845 | 0.53 | 179.0 kb | 2.52:1 | 375,800 | 195,245 | 0.52 | 177.0 kb | 2.59:1 |  
| meta6 | 63,183 | 32,992 | 0.52 | 185.0 kb | 1:1.13 | 15,659 | 8,228 | 0.53 | 181.0 kb | 1:1.12 |  
| meta7 | 32,306 | 27,353 | 0.85 | 179.0 kb | 1:1.21 | 14,661 | 12,176 | 0.83 | 175.0 kb | 1:1.23 |  
| mega | 452,795 | 281,854 | 0.62 | 202.0 kb | 1:1.07 | 491,079 | 272,795 | 0.56 | 178.0 kb | 1.81:1 |  

![Rare CNV stats](https://storage.googleapis.com/rcnv_project/public/rCNV.metacohort.stats.jpg)  

The information for these tables was collected using `collect_cohort_stats.sh` and visualized using `plot_cnv_stats_per_cohort.R`.  

---  

## Noncoding subsets  

For certain analyses, we restricted this dataset further to rCNVs unlikely to directly disrupt disease-relevant protein-coding genes.  

These subsets are defined as follows:  

1. **Strictly noncoding**: all rCNVs were excluded based on any overlap with any canonical exon from any protein-coding gene ([as described here](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/gene#gene-definitions)).  

2. **Loose noncoding**: to improve power for [noncoding association testing](https://github.com/talkowski-lab/rCNV2/tree/master/analysis/noncoding), we supplemented the `Strictly noncoding` subset (described above) with any rCNVs that overlapped exons from the 4,793 genes meeting both of the following criteria:  
* Known to readily tolerate functional mutations in the general population ([described here](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/gene#gene-set-definitions)); and  
* No known disease association [per OMIM](https://www.omim.org/).  

This second subset ("`loose noncoding`") represents the subset of rare CNVs depleted for strong-effect, disease-relevant coding effects.  

For both sets of noncoding CNVs, we padded all exons from all genes by ±50kb to protect against CNV breakpoint imprecision.  

The code to apply these filters is contained in `extract_noncoding_subsets.sh`.  

### Strictly noncoding rare CNV callset properties  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| meta1 | 29,104 | 4,029 | 0.14 | 187.0 kb | 1.42:1 | 19,585 | 4,329 | 0.22 | 157.0 kb | 2.07:1 |  
| meta2 | 79,196 | 9,331 | 0.12 | 182.0 kb | 1.67:1 | 8,670 | 1,176 | 0.14 | 159.0 kb | 1.70:1 |  
| meta3 | 41,065 | 10,180 | 0.25 | 158.0 kb | 1.80:1 | 32,543 | 7,372 | 0.23 | 163.0 kb | 1.84:1 |  
| meta4 | 153,870 | 31,678 | 0.21 | 167.0 kb | 2.03:1 | 24,161 | 6,535 | 0.27 | 168.0 kb | 2.41:1 |  
| meta5 | 54,071 | 9,340 | 0.17 | 163.0 kb | 6.01:1 | 375,800 | 65,508 | 0.17 | 163.0 kb | 5.86:1 |  
| meta6 | 63,183 | 12,001 | 0.19 | 165.0 kb | 1.99:1 | 15,659 | 3,021 | 0.19 | 164.0 kb | 1.94:1 |  
| meta7 | 32,306 | 9,356 | 0.29 | 167.0 kb | 2.08:1 | 14,661 | 4,198 | 0.29 | 167.0 kb | 1.98:1 |  

![Strictly noncoding rare CNV stats](https://storage.googleapis.com/rcnv_project/public/strict_noncoding.metacohort.stats.jpg)  

### Loose noncoding rare CNV callset properties  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| meta1 | 29,104 | 4,973 | 0.17 | 184.0 kb | 1.38:1 | 19,585 | 5,171 | 0.26 | 158.0 kb | 1.91:1 |  
| meta2 | 79,196 | 10,935 | 0.14 | 182.0 kb | 1.59:1 | 8,670 | 1,364 | 0.16 | 158.0 kb | 1.49:1 |  
| meta3 | 41,065 | 12,053 | 0.29 | 158.0 kb | 1.63:1 | 32,543 | 8,896 | 0.27 | 162.0 kb | 1.69:1 |  
| meta4 | 153,870 | 37,346 | 0.24 | 168.0 kb | 1.92:1 | 24,161 | 7,697 | 0.32 | 168.0 kb | 2.28:1 |  
| meta5 | 54,071 | 11,306 | 0.21 | 162.0 kb | 4.79:1 | 375,800 | 79,015 | 0.21 | 162.0 kb | 4.83:1 |  
| meta6 | 63,183 | 13,797 | 0.22 | 165.0 kb | 1.87:1 | 15,659 | 3,496 | 0.22 | 162.0 kb | 1.82:1 |  
| meta7 | 32,306 | 10,696 | 0.33 | 170.0 kb | 1.89:1 | 14,661 | 4,786 | 0.33 | 170.0 kb | 1.79:1 |  

![Loose noncoding rare CNV stats](https://storage.googleapis.com/rcnv_project/public/loose_noncoding.metacohort.stats.jpg)  

