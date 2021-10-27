# CNV Data Curation

We aggregated existing microarray-based CNV calls from numerous sources. This process is described below.  

All commands executed to filter the CNV data as descrbed below are contained in `filter_CNV_data.sh`.  

In practice, the commands in `filter_CNV_data.sh` were parallelized in [FireCloud/Terra](https://portal.firecloud.org) with `filter_CNV_data.wdl`.  

## CNV data sources  

We aggregated CNV data from multiple sources, listed alphabetically below.  

Note that the sample counts in this table represent the number of samples retained after filtering outliers (described below).  

| Cohort | Citation | PMID | Platform(s) | Build | Phenos | N Cases | N Ctrls |
| --- | :--- | :--- | :--- | :--- | :--- | ---: | ---: |
| Boston Children's Hospital (`BCH`) | [Talkowski _et al._, _Cell_ (2012)](https://pubmed.ncbi.nlm.nih.gov/22521361) | [22521361](https://pubmed.ncbi.nlm.nih.gov/22521361) | aCGH | hg18 | Mixed | 3591 | 0 |  
| BioVU (`BioVU`) | [Roden _et al._, _Clin. Pharmacol. Ther._ (2008)](https://pubmed.ncbi.nlm.nih.gov/18500243) | [18500243](https://pubmed.ncbi.nlm.nih.gov/18500243) | Illumina MegaEx (100%) | hg19 | Mixed | 32306 | 14661 |  
| Children's Hospital of Philadelphia (`CHOP`) | [Li _et al._, _Nat. Commun._ (2020)](https://pubmed.ncbi.nlm.nih.gov/31937769) | [31937769](https://pubmed.ncbi.nlm.nih.gov/31937769) | Mixed Illumina SNP genotyping platforms | hg19 | Mixed | 153870 | 24161 |  
| Eichler Lab (`Coe`) | [Coe _et al._, _Nat. Genet._ (2014)](https://pubmed.ncbi.nlm.nih.gov/25217958) | [25217958](https://pubmed.ncbi.nlm.nih.gov/25217958) | Cases: SignatureChip OS v2.0 (58%), SignatureChip OS v1.0 (34%), Other (8%); Controls: Affy 6.0 (100%) | hg19 | Developmental disorders | 29083 | 11256 |  
| Eichler Lab (`Cooper`) | [Cooper _et al._, _Nat. Genet._ (2011)](https://pubmed.ncbi.nlm.nih.gov/21841781) | [21841781](https://pubmed.ncbi.nlm.nih.gov/21841781) | Illumina 550k-610k (75%), Custom 1.2M (25%) | hg19 | N/A | 0 | 8329 |  
| Epi25 Consortium (`Epi25k`) | [Niestroj _et al._, _Brain_ (2020)](https://pubmed.ncbi.nlm.nih.gov/32568404) | [32568404](https://pubmed.ncbi.nlm.nih.gov/32568404) | Illumina GSA-MD v1.0 (100%) | hg19 | Epilepsy | 12053 | 8173 |  
| Estonian Biobank (`EstBB`) | [Leitsalu _et al._, _Int. J. Epidemiol._ (2014)](https://pubmed.ncbi.nlm.nih.gov/24518929) | [24518929](https://pubmed.ncbi.nlm.nih.gov/24518929) | Illumina GSA (100%) | hg19 | Mixed | 63183 | 15659 |  
| GeneDX (`GDX`) | - | - | Custom Aglient SNP Array (71%), CytoScan HD (29%) | hg18 & hg19 | Mixed | 74208 | 0 |  
| Indiana University (`IU`) | - | - | CytoScan HD (100%) | hg19 | Mixed | 1576 | 0 |  
| Ontario Population Genomics Platform (`Ontario`) | [Uddin _et al._, _Genet. Med._ (2015)](https://pubmed.ncbi.nlm.nih.gov/25503493) | [25503493](https://pubmed.ncbi.nlm.nih.gov/25503493) | CytoScan HD (100%) | hg19 | N/A | 0 | 873 |  
| Psychiatric Genetics Consortium (`PGC`) | [Marshall _et al._, _Nat. Genet._ (2017)](https://pubmed.ncbi.nlm.nih.gov/27869829) | [27869829](https://pubmed.ncbi.nlm.nih.gov/27869829) | Affy 6.0 (37%), Omni Express (31%), Omni Express Plus (12%), Other (20%) | hg18 | Schizophrenia | 21094 | 20277 |  
| Radboud University Medical Center (`RUMC`) | [Vulto-van Silfhout _et al._, _Hum. Mutat._ (2013)](https://pubmed.ncbi.nlm.nih.gov/24038936) | [24038936](https://pubmed.ncbi.nlm.nih.gov/24038936) | Affy 250k (100%) | hg17 | Intellectual disability | 5531 | 0 |  
| SickKids Hospital (`SickKids`) | [Zarrei _et al._, _NPJ Genomic Medicine_ (2019)](https://pubmed.ncbi.nlm.nih.gov/31602316) | [31602316](https://pubmed.ncbi.nlm.nih.gov/31602316) | Affy 6.0 (100%) | hg19 | Developmental disorders | 2689 | 0 |  
| Simons Simplex Collection (`SSC`) | [Sanders _et al._, _Neuron_ (2015)](https://pubmed.ncbi.nlm.nih.gov/26402605) | [26402605](https://pubmed.ncbi.nlm.nih.gov/26402605) | Omni 1Mv3 (46%), Omni 2.5 (41%), Omni 1Mv1 (13%) | hg18 | Autism | 2795 | 0 |  
| The Cancer Genome Atlas (`TCGA`) | [Zack _et al._, _Nat. Genet._ (2013)](https://pubmed.ncbi.nlm.nih.gov/24071852) | [24071852](https://pubmed.ncbi.nlm.nih.gov/24071852) | Affy 6.0 (100%) | hg19 | N/A | 0 | 8670 |  
| The Genetic Etiology of Tourette Syndrome Consortium (`TSAICG`) | [Huang _et al._, _Neuron (2017)](https://pubmed.ncbi.nlm.nih.gov/28641109) | [28641109](https://pubmed.ncbi.nlm.nih.gov/28641109) | OmniExpress (100%) | hg19 | Tourette Syndrome | 2434 | 4093 |  
| UK Biobank (`UKBB`) | [Macé _et al._, _Nat. Commun._ (2017)](https://pubmed.ncbi.nlm.nih.gov/28963451) | [28963451](https://pubmed.ncbi.nlm.nih.gov/28963451) | UKBB Affy Axiom (100%) | hg19 | Mixed | 54071 | 375800 |  

## Raw CNV data processing steps  

All CNV data native to hg17 or hg18 was lifted over to hg19 using UCSC liftOver, requiring at least 50% of the original CNV to map successfully to hg19 in order to be retained.  

Some datasets required manual curation prior to inclusion. Where necessary, these steps are enumerated below:  

 * **BioVU**: TBD [TODO: ADD TEXT HERE]
 * **CHOP**: CNVs were filtered on quality score ≥40 and CNV size ≥25kb while requiring at least 10 SNPs per CNV. After CNV filtering, samples with `LRR_SD` <0.25, >20 CNV calls, or SNP call rate <98% were excluded as outliers, as well as samples genotyped on arrays with <175k SNP probes or samples labeled as cancer or Down's Syndrome patients. Finally, we identified 19 loci with apparently platform-specific artifactual CNV pileups. CNVs covered ≥10% by any of these artifact regions were removed from the callset. Lists of CHOP-specific blacklisted loci for deletions and duplications are provided as [reference files](https://github.com/talkowski-lab/rCNV2/tree/master/refs).  
 * **Cooper**: Only retained control samples from Cooper _et al._ All cases from Cooper _et al._ also appear in Coe _et al._  
 * **Epi25k**: CNVs were filtered on ≥10 probes and ≥25kb. After CNV filtering, samples with >25 CNV calls were excluded as outliers.  
 * **EstBB**: Samples were excluded if were not included in SNP imputation, had genotype calls missing for ≥2% of sites, or belonged to two genotyping batches based on visual inspection of genotyping intensity parameters, followed by further exclusion of genotyping plates (≤24 samples per plate) that contained >3 samples with >200 CNV calls. We retained unrelated samples that had been linked to Estonian health registries and that had ≤50 raw CNV calls. We included CNVs with a quality score ≥15, were covered by ≥10 probes, and were ≥25kb in size. Finally, we pruned related samples and any samples with known malignant cancers or chromosomal disorders (e.g., Down's Syndrome or sex chromosome aneuploidies).
 * **GDX**: All CNVs were required to be ≥20kb and <40Mb in length. Except for the minority of 9,958 samples for which additional CNV call metadata was unavailable, all CNVs were further required to not have been annotated as a suspected false positive or mosaic event, have estimated copy numbers ≤1.5 for deletions or ≥2.5 for duplications, include ≥10 probes and have P(CNV) ≤ 10<sup>-10</sup>. Following CNV call filtering, we excluded all samples that either had >10 calls each, were identified as potential biological replicates, were referred for testing due to being a relative of a known carrier of a medically relevant CNV, or had “advanced maternal age” as their indication for testing.  
 * **IU**: Excluded samples from Indiana University (IU) cohort derived from buccal swab DNA, samples with known aneuploidies or large runs of homozygosity, and samples with no phenotypic indication specified.  
 * **SSC**: CNVs were filtered on pCNV ≤10<sup>-9</sup>, per recommendation of the authors. We only retained affected individuals since all controls were first-degree relatives of affected cases.    
 * **SickKids**: CNVs were filtered on ≥25kb. After CNV filtering, samples with >80 CNV calls were excluded as outliers. Finally, we identified a single locus on chr12 that had CNVs only appearing in ADHD samples at 2.8% frequency; these CNVs were removed from the callset. We only retained affected individuals since all controls were first-degree relatives of affected cases.  
 * **TCGA**: CNVs were filtered on ≥10 probes and ≥25kb. Deletions were required to have a mean intensity ≤ log<sub>2</sub>(0.6) and duplications were required to have a mean intensity ≥  log<sub>2</sub>(1.45). We only retained normal samples from TCGA tumor:normal pairs, and excluded any TCGA donors with known blood cancer.  
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

Seven studies were unable to be defragmented due to inadequate sample-level information: `Coe` (_controls only_), `Cooper`, `PGC`, `BioVU`, `TSAICG`, `Ontario`, and `RUMC`.

### Raw CNV callset properties  

The properties of each callset are listed below after initial data processing steps but prior to further filtering.  

As [described below](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#case-control-metacohorts), we subdivided 13,139 control samples from the UKBB cohort to use as proxy controls for cases from GeneDx and retained the remaining 416,732 UKBB samples as a separate cohort; these cohort subsets are referred to as `UKBB_main` and `UKBB_sub` for the main UKBB cohort and the subset of 13,139 controls, respectively.  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| BCH | 3,591 | 5,159 | 1.44 | 204.0 kb | 1:1.27 | 0 | 0 | 0 | NA | NA |  
| CHOP | 153,870 | 641,938 | 4.17 | 87.0 kb | 1.05:1 | 24,161 | 108,559 | 4.49 | 87.0 kb | 1.86:1 |  
| Coe | 29,104 | 28,667 | 0.98 | 187.0 kb | 1:1.11 | 11,256 | 273,331 | 24.28 | 53.0 kb | 1.24:1 |  
| Cooper | 0 | 0 | 0 | NA | NA | 8,329 | 432,478 | 51.92 | 2.0 kb | 8.04:1 |  
| Epi25k | 12,053 | 92,412 | 7.67 | 97.0 kb | 1:1.54 | 8,173 | 64,036 | 7.84 | 100.0 kb | 1:1.33 |  
| GDX | 74,028 | 145,459 | 1.96 | 202.0 kb | 1:1.59 | 0 | 0 | 0 | NA | NA |  
| IU | 1,577 | 4,890 | 3.1 | 98.0 kb | 1.41:1 | 0 | 0 | 0 | NA | NA |  
| Ontario | 0 | 0 | 0 | NA | NA | 873 | 71,178 | 81.53 | 10.0 kb | 4.11:1 |  
| PGC | 21,094 | 42,096 | 2.0 | 80.0 kb | 1:1.04 | 20,277 | 40,464 | 2.0 | 78.0 kb | 1:1.03 |  
| RUMC | 5,531 | 1,777 | 0.32 | 1090.0 kb | 1:1.00 | 0 | 0 | 0 | NA | NA |  
| SSC | 2,795 | 30,856 | 11.04 | 21.0 kb | 3.09:1 | 0 | 0 | 0 | NA | NA |  
| SickKids | 2,689 | 72,278 | 26.88 | 55.0 kb | 1.04:1 | 0 | 0 | 0 | NA | NA |  
| TCGA | 0 | 0 | 0 | NA | NA | 8,670 | 206,758 | 23.85 | 60.0 kb | 1:1.16 |  
| TSAICG | 2,434 | 3,541 | 1.45 | 91.0 kb | 1.01:1 | 4,093 | 5,834 | 1.43 | 91.0 kb | 1:1.08 |  
| UKBB_main | 54,071 | 140,082 | 2.59 | 88.0 kb | 5.62:1 | 362,661 | 930,043 | 2.56 | 87.0 kb | 5.69:1 |  
| UKBB_sub | 0 | 0 | 0 | NA | NA | 13,139 | 33,913 | 2.58 | 86.0 kb | 5.64:1 |  
| EstBB | 63,183 | 226,083 | 3.58 | 83.0 kb | 1:1.45 | 15,659 | 56,267 | 3.59 | 83.0 kb | 1:1.46 |  
| BioVU | 32,306 | 58,955 | 1.82 | 133.0 kb | 1:1.32 | 14,661 | 27,411 | 1.87 | 132.0 kb | 1:1.40 |  

The information for this table was collected using `collect_cohort_stats.sh`, and is visualized below using `plot_cnv_stats_per_cohort.R`.  

![Raw CNV stats](https://storage.googleapis.com/rcnv_project/public/raw_cnv.stats.jpg)  

### Raw data access  

All raw CNV data files and their tabix indexes are stored in a protected Google Cloud bucket, here:  
```
$ gsutil ls gs://rcnv_project/raw_data/cnv/
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
3. CNVs with ≥75% reciprocal overlap versus [known genomic disorder CNV loci](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/other#genomic-disorders) reported by at least two of six established resources were exempted from this filter step, as we reasoned it was possible that certain genomic disorder CNVs may be represented at frequencies close to 1% in some case-only subsets (e.g., SSC, RUMC, SickKids, GeneDx, _etc._).  
4. "Substantial" coverage determined based on ≥50% coverage per BEDTools coverage ([Quinlan & Hall, _Bioinformatics_ (2010)](https://www.ncbi.nlm.nih.gov/pubmed/20110278)).

### Rare CNV callset properties

The properties of each rare CNV callset are listed below after the above filtering steps.  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| BCH | 3,591 | 3,028 | 0.84 | 342.0 kb | 1:1.40 | 0 | 0 | 0 | NA | NA |  
| CHOP | 153,870 | 91,184 | 0.59 | 193.0 kb | 1:1.08 | 24,161 | 16,819 | 0.7 | 186.0 kb | 1:1.01 |  
| Coe | 29,104 | 19,640 | 0.67 | 264.0 kb | 1:1.18 | 11,256 | 9,279 | 0.82 | 177.0 kb | 1:1.59 |  
| Cooper | 0 | 0 | 0 | NA | NA | 8,329 | 4,183 | 0.5 | 170.0 kb | 1.38:1 |  
| Epi25k | 12,053 | 8,989 | 0.75 | 176.0 kb | 1:1.07 | 8,173 | 5,611 | 0.69 | 185.0 kb | 1:1.14 |  
| GDX | 74,028 | 45,051 | 0.61 | 283.0 kb | 1:1.47 | 0 | 0 | 0 | NA | NA |  
| IU | 1,577 | 1,651 | 1.05 | 209.0 kb | 1:1.34 | 0 | 0 | 0 | NA | NA |  
| Ontario | 0 | 0 | 0 | NA | NA | 873 | 745 | 0.85 | 182.0 kb | 1:1.50 |  
| PGC | 21,094 | 14,660 | 0.69 | 190.0 kb | 1:1.34 | 20,277 | 13,272 | 0.65 | 182.0 kb | 1:1.33 |  
| RUMC | 5,531 | 1,486 | 0.27 | 1059.0 kb | 1.03:1 | 0 | 0 | 0 | NA | NA |  
| SSC | 2,795 | 2,151 | 0.77 | 195.0 kb | 1:1.15 | 0 | 0 | 0 | NA | NA |  
| SickKids | 2,689 | 3,879 | 1.44 | 181.0 kb | 1:1.35 | 0 | 0 | 0 | NA | NA |  
| TCGA | 0 | 0 | 0 | NA | NA | 8,670 | 8,710 | 1.0 | 178.0 kb | 1:1.71 |  
| TSAICG | 2,434 | 1,607 | 0.66 | 198.0 kb | 1:1.31 | 4,093 | 2,663 | 0.65 | 187.0 kb | 1:1.36 |  
| UKBB_main | 54,071 | 28,690 | 0.53 | 180.0 kb | 2.50:1 | 362,661 | 187,626 | 0.52 | 177.0 kb | 2.57:1 |  
| UKBB_sub | 0 | 0 | 0 | NA | NA | 13,139 | 6,758 | 0.51 | 176.0 kb | 2.62:1 |  
| EstBB | 63,183 | 32,868 | 0.52 | 185.0 kb | 1:1.14 | 15,659 | 8,206 | 0.52 | 182.0 kb | 1:1.13 |  
| BioVU | 32,306 | 27,239 | 0.84 | 179.0 kb | 1:1.22 | 14,661 | 12,118 | 0.83 | 176.0 kb | 1:1.25 |  

The information for this table was collected using `collect_cohort_stats.sh`, and is visualized below using `plot_cnv_stats_per_cohort.R`.  

![Rare CNV stats](https://storage.googleapis.com/rcnv_project/public/rCNV.stats.jpg)  

---  

## Case-control "metacohorts"  

For analyses, we combined CNV data from multiple cohorts into seven matched groups, dubbed **metacohorts**, to control for technical differences between individual data sources and cohorts.  

Individual cohorts were assigned to metacohorts on the basis of similarity in sample recruitment protocols, microarray probesets, and other technical factors.  

We did not modify cohorts that had matched cases and controls with one exception: as the UKBB Axiom array design most closely matched the predominant array platform used by the GeneDx cohort, we randomly selected 13,139 controls from the UKBB (dubbed `UKBB_sub`) to group with GeneDx cases in the `meta2` metacohort. We retained all remaining UKBB samples (dubbed `UKBB_main`) in the `meta5` metacohort.  

These metacohorts represent the basic unit on which all burden testing was performed, and are described in the table below.  

| Metacohort ID | Case Source(s) | Number of Cases | Control Sources(s) | Number of Controls |  
| :--- | :--- | ---: | :--- | ---: |  
| `meta1` | BCH, Coe, IU | 34,272 | Coe, Cooper | 19,585 |  
| `meta2` | Epi25k, GDX, TSAICG | 88,515 | Epi25k, Ontario, TSAICG, UKBB_sub | 26,278 |  
| `meta3` | PGC, RUMC, SSC, SickKids | 32,109 | PGC, TCGA | 28,947 |  
| `meta4` | CHOP | 153,870 | CHOP | 24,161 |  
| `meta5` | UKBB_main | 54,071 | UKBB_main | 362,661 |  
| `meta6` | EstBB | 63,183 | EstBB | 15,659 |  
| `meta7` | BioVU | 32,306 | BioVU | 14,661 |  

### Metacohort rare CNV callset properties  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| meta1 | 34,272 | 24,319 | 0.71 | 266.0 kb | 1:1.22 | 19,585 | 13,462 | 0.69 | 174.0 kb | 1:1.24 |  
| meta2 | 88,515 | 55,647 | 0.63 | 256.0 kb | 1:1.39 | 26,278 | 15,777 | 0.6 | 181.0 kb | 1.31:1 |  
| meta3 | 32,109 | 22,176 | 0.69 | 201.0 kb | 1:1.29 | 28,947 | 21,982 | 0.76 | 180.0 kb | 1:1.47 |  
| meta4 | 153,870 | 91,184 | 0.59 | 193.0 kb | 1:1.08 | 24,161 | 16,819 | 0.7 | 186.0 kb | 1:1.01 |  
| meta5 | 54,071 | 28,690 | 0.53 | 180.0 kb | 2.50:1 | 362,661 | 187,626 | 0.52 | 177.0 kb | 2.57:1 |  
| meta6 | 63,183 | 32,868 | 0.52 | 185.0 kb | 1:1.14 | 15,659 | 8,206 | 0.52 | 182.0 kb | 1:1.13 |  
| meta7 | 32,306 | 27,239 | 0.84 | 179.0 kb | 1:1.22 | 14,661 | 12,118 | 0.83 | 176.0 kb | 1:1.25 |  

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
| meta1 | 34,272 | 5,306 | 0.15 | 192.0 kb | 1.32:1 | 19,585 | 4,588 | 0.23 | 158.0 kb | 1.91:1 |  
| meta2 | 88,515 | 12,644 | 0.14 | 174.0 kb | 1.65:1 | 26,278 | 5,896 | 0.22 | 166.0 kb | 2.49:1 |  
| meta3 | 32,109 | 7,070 | 0.22 | 163.0 kb | 1.63:1 | 28,947 | 7,309 | 0.25 | 158.0 kb | 1.89:1 |  
| meta4 | 153,870 | 33,003 | 0.21 | 166.0 kb | 1.85:1 | 24,161 | 6,590 | 0.27 | 165.0 kb | 2.14:1 |  
| meta5 | 54,071 | 9,876 | 0.18 | 164.0 kb | 5.07:1 | 362,661 | 66,811 | 0.18 | 164.0 kb | 4.98:1 |  
| meta6 | 63,183 | 12,339 | 0.2 | 165.0 kb | 1.89:1 | 15,659 | 3,094 | 0.2 | 165.0 kb | 1.86:1 |  
| meta7 | 32,306 | 9,381 | 0.29 | 167.0 kb | 2.03:1 | 14,661 | 4,211 | 0.29 | 167.0 kb | 1.93:1 |  

![Strictly noncoding rare CNV stats](https://storage.googleapis.com/rcnv_project/public/strict_noncoding.metacohort.stats.jpg)  

### Loose noncoding rare CNV callset properties  

| Dataset | N Cases | Case CNVs | CNVs /Case | Case Median Size | Case DEL:DUP | N Ctrls | Ctrl CNVs | CNVs /Ctrl | Ctrl Median Size | Ctrl DEL:DUP |  
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| meta1 | 34,272 | 6,449 | 0.19 | 189.0 kb | 1.30:1 | 19,585 | 5,485 | 0.28 | 159.0 kb | 1.79:1 |  
| meta2 | 88,515 | 14,786 | 0.17 | 174.0 kb | 1.57:1 | 26,278 | 7,049 | 0.27 | 165.0 kb | 2.32:1 |  
| meta3 | 32,109 | 8,344 | 0.26 | 162.0 kb | 1.50:1 | 28,947 | 8,709 | 0.3 | 159.0 kb | 1.74:1 |  
| meta4 | 153,870 | 38,817 | 0.25 | 168.0 kb | 1.78:1 | 24,161 | 7,782 | 0.32 | 166.0 kb | 2.04:1 |  
| meta5 | 54,071 | 11,876 | 0.22 | 162.0 kb | 4.26:1 | 362,661 | 80,208 | 0.22 | 162.0 kb | 4.33:1 |  
| meta6 | 63,183 | 14,209 | 0.22 | 165.0 kb | 1.78:1 | 15,659 | 3,588 | 0.23 | 165.0 kb | 1.75:1 |  
| meta7 | 32,306 | 10,941 | 0.34 | 168.0 kb | 1.91:1 | 14,661 | 4,888 | 0.33 | 168.0 kb | 1.80:1 |  

![Loose noncoding rare CNV stats](https://storage.googleapis.com/rcnv_project/public/loose_noncoding.metacohort.stats.jpg)  

