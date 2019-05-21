# Phenotype Data Curation

We grouped samples into phenotype-matched subsets using a standardized hierarchical ontology.  

This process is described below.  

### Phenotype curation  

Given the variability of phenotypic detail available between different cohorts, we applied a standardized hierarchical phenotype consolidation scheme uniformly across all samples.  

We used the [Human Phenotype Ontology (HPO)](https://hpo.jax.org/app/) to harmonize phenotype data across all [cohorts included in the overall CNV dataset](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV/).  

For each sample, we recursively matched all phenotype information available against HPO keywords, and assigned the corresponding HPO terms to each sample when a match was found.  

The code to perform this HPO conversion is contained in `convert_phenotypes_to_hpo.sh`.  

### Phenotype hierarchy  

After assigning HPO terms to each sample, we next subgrouped samples based on concordant HPO terms into a reduced hierarchy of "phenotype groups" considered in this study.  

After tabulating the number of samples matching each HPO code, we retained all HPO terms with at least 1,000 samples. 

We next compared shared sample memberships between all pairs of related HPO terms:
 * If fewer than 1,000 samples differed between a pair of HPO terms, we retained the more general (i.e., higher-level) term, and excluded the more specific (i.e., lower-level) term.
 * If both terms were equally high-level, we retained the term with the larger sample size.

This process yielded a hierarchical phenotype classification system with 30 distinct HPO-based sample groups, each which has ≥ 1,000 samples and differs from another phenotype group by at least 1,000 distinct samples.  

The code to perform the phenotype group consolidation process described above is contained in `convert_phenotypes_to_hpo.sh`    

The final set of phenotypes considered is detailed below:  

| HPO Term | Description | Samples | HPO Tier | Parent Terms | Child Terms |  
| :--- | :--- | ---: | ---: | :--- | :--- |  
| HP:0000118 | Phenotypic abnormality | 223,116 | 1 | NA | HP:0012638, HP:0004305, HP:0003011, HP:0012758, HP:0025031, HP:0001510, HP:0000752, HP:0001250, HP:0100753, HP:0002715, HP:0012759, HP:0000736, HP:0012639, HP:0000707, HP:0000152, HP:0000708, HP:0001507, HP:0001627, HP:0001249, HP:0001626, HP:0007367, HP:0002597, HP:0012443, HP:0002197, HP:0000924, HP:0100022, HP:0031466, HP:0000729 |  
| HP:0000707 | Abnormality of the nervous system | 116,783 | 2 | HP:0000118 | HP:0100022, HP:0012758, HP:0007367, HP:0100753, HP:0000752, HP:0000708, HP:0012443, HP:0012759, HP:0000736, HP:0002197, HP:0012638, HP:0031466, HP:0004305, HP:0012639, HP:0001249, HP:0001250, HP:0000729 |  
| HP:0012638 | Abnormality of nervous system physiology | 91,067 | 3 | HP:0000707, HP:0000118 | HP:0100022, HP:0000708, HP:0012759, HP:0000736, HP:0002197, HP:0031466, HP:0004305, HP:0012758, HP:0100753, HP:0000752, HP:0001249, HP:0001250, HP:0000729 |  
| HP:0000708 | Behavioral abnormality | 54,041 | 4 | HP:0000707, HP:0000118, HP:0012638 | HP:0100753, HP:0031466, HP:0000736, HP:0000729 |  
| UNKNOWN | NA | 49,137 | 2 | HP:0000118 | NA |  
| HP:0012759 | Neurodevelopmental abnormality | 36,378 | 4 | HP:0000707, HP:0000118, HP:0012638 | HP:0012758, HP:0001249 |  
| HP:0012758 | Neurodevelopmental delay | 34,633 | 5 | HP:0000707, HP:0012759, HP:0012638, HP:0000118 | NA |  
| HP:0002715 | Abnormality of the immune system | 32,932 | 2 | HP:0000118 | NA |  
| HP:0012639 | Abnormality of nervous system morphology | 26,277 | 3 | HP:0000707, HP:0000118 | HP:0012443, HP:0007367 |  
| HP:0007367 | Atrophy/Degeneration affecting the central nervous system | 24,257 | 5 | HP:0000707, HP:0012639, HP:0000118 | NA |  
| HP:0001626 | Abnormality of the cardiovascular system | 22,372 | 2 | HP:0000118 | HP:0002597, HP:0001627 |  
| HP:0100753 | Schizophrenia | 22,301 | 5 | HP:0000707, HP:0000708, HP:0012638, HP:0000118 | NA |  
| HP:0000729 | Autistic behavior | 16,926 | 5 | HP:0000707, HP:0000708, HP:0012638, HP:0000118 | NA |  
| HP:0100022 | Abnormality of movement | 15,553 | 4 | HP:0000707, HP:0000118, HP:0012638 | HP:0000752, HP:0004305 |  
| HP:0002597 | Abnormality of the vasculature | 15,073 | 3 | HP:0000118, HP:0001626 | NA |  
| HP:0000752 | Hyperactivity | 12,941 | 5 | HP:0000707, HP:0100022, HP:0012638, HP:0000118 | NA |  
| HP:0000736 | Short attention span | 12,908 | 5 | HP:0000707, HP:0000708, HP:0012638, HP:0000118 | NA |  
| HP:0000152 | Abnormality of head or neck | 7,045 | 2 | HP:0000118 | NA |  
| HP:0001627 | Abnormal heart morphology | 6,751 | 4 | HP:0000118, HP:0001626 | NA |  
| HP:0001250 | Seizures | 4,965 | 4 | HP:0000707, HP:0000118, HP:0012638 | HP:0002197 |  
| HP:0001507 | Growth abnormality | 3,262 | 2 | HP:0000118 | HP:0001510 |  
| HP:0004305 | Involuntary movements | 2,538 | 5 | HP:0000707, HP:0100022, HP:0012638, HP:0000118 | NA |  
| HP:0031466 | Impairment in personality functioning | 2,171 | 5 | HP:0000707, HP:0000708, HP:0012638, HP:0000118 | NA |  
| HP:0001249 | Intellectual disability | 2,072 | 5 | HP:0000707, HP:0012759, HP:0012638, HP:0000118 | NA |  
| HP:0000924 | Abnormality of the skeletal system | 1,986 | 2 | HP:0000118 | NA |  
| HP:0012443 | Abnormality of brain morphology | 1,950 | 5 | HP:0000707, HP:0012639, HP:0000118 | NA |  
| HP:0002197 | Generalized-onset seizure | 1,781 | 5 | HP:0000707, HP:0001250, HP:0012638, HP:0000118 | NA |  
| HP:0025031 | Abnormality of the digestive system | 1,588 | 2 | HP:0000118 | NA |  
| HP:0003011 | Abnormality of the musculature | 1,376 | 2 | HP:0000118 | NA |  
| HP:0001510 | Growth delay | 1,035 | 3 | HP:0000118, HP:0001507 | NA |  

### HPO terms per cohort

The number of samples per phenotype group per cohort is outlined in the table, below.  

The code to generate this table is provided in `gather_hpo_per_cohort_table.py`.  

| HPO | description | Total | PGC | Cooper | Coe | SSC | UKBB | CHOP | GDX | TSAICG | BCH | TCGA |  
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| HP:0000118 | Phenotypic abnormality | 223,116 | 21,094 | 0 | 29,104 | 2,795 | 0 | 154,140 | 9,958 | 2,434 | 3,591 | 0 |  
| HP:0000707 | Abnormality of the nervous system | 116,783 | 21,094 | 0 | 29,104 | 2,795 | 0 | 54,406 | 4,036 | 2,434 | 2,914 | 0 |  
| HP:0012638 | Abnormality of nervous system physiology | 91,067 | 21,094 | 0 | 29,104 | 2,795 | 0 | 28,948 | 3,846 | 2,434 | 2,846 | 0 |  
| HP:0000708 | Behavioral abnormality | 54,041 | 21,094 | 0 | 1,834 | 2,795 | 0 | 26,364 | 730 | 0 | 1,224 | 0 |  
| UNKNOWN | NA | 49,137 | 0 | 0 | 0 | 0 | 0 | 43,873 | 4,971 | 0 | 293 | 0 |  
| HP:0012759 | Neurodevelopmental abnormality | 36,378 | 0 | 0 | 29,104 | 814 | 0 | 671 | 3,423 | 0 | 2,366 | 0 |  
| HP:0012758 | Neurodevelopmental delay | 34,633 | 0 | 0 | 29,104 | 0 | 0 | 0 | 3,320 | 0 | 2,209 | 0 |  
| HP:0002715 | Abnormality of the immune system | 32,932 | 0 | 0 | 19 | 0 | 0 | 32,843 | 36 | 0 | 34 | 0 |  
| HP:0012639 | Abnormality of nervous system morphology | 26,277 | 0 | 0 | 352 | 0 | 0 | 24,936 | 712 | 0 | 277 | 0 |  
| HP:0007367 | Atrophy/Degeneration affecting the central nervous system | 24,257 | 0 | 0 | 4 | 0 | 0 | 24,247 | 4 | 0 | 2 | 0 |  
| HEALTHY_CONTROL | Unaffected control sample | 561,196 | 20,277 | 8,329 | 11,256 | 0 | 480,501 | 28,070 | 0 | 4,093 | 0 | 8,670 |  
| HP:0001626 | Abnormality of the cardiovascular system | 22,372 | 0 | 0 | 748 | 0 | 0 | 21,003 | 409 | 0 | 212 | 0 |  
| HP:0100753 | Schizophrenia | 22,301 | 21,094 | 0 | 2 | 0 | 0 | 1,200 | 5 | 0 | 0 | 0 |  
| HP:0000729 | Autistic behavior | 16,926 | 0 | 0 | 1,569 | 2,795 | 0 | 11,071 | 638 | 0 | 853 | 0 |  
| HP:0100022 | Abnormality of movement | 15,553 | 0 | 0 | 317 | 0 | 0 | 12,394 | 83 | 2,434 | 325 | 0 |  
| HP:0002597 | Abnormality of the vasculature | 15,073 | 0 | 0 | 50 | 0 | 0 | 14,948 | 59 | 0 | 16 | 0 |  
| HP:0000752 | Hyperactivity | 12,941 | 0 | 0 | 281 | 0 | 0 | 12,359 | 63 | 0 | 238 | 0 |  
| HP:0000736 | Short attention span | 12,908 | 0 | 0 | 264 | 0 | 0 | 12,359 | 58 | 0 | 227 | 0 |  
| HP:0000152 | Abnormality of head or neck | 7,045 | 0 | 0 | 4,535 | 0 | 0 | 0 | 1,662 | 0 | 848 | 0 |  
| HP:0001627 | Abnormal heart morphology | 6,751 | 0 | 0 | 166 | 0 | 0 | 6,055 | 334 | 0 | 196 | 0 |  
| HP:0001250 | Seizures | 4,965 | 0 | 0 | 1,807 | 318 | 0 | 1,762 | 539 | 0 | 539 | 0 |  
| HP:0001507 | Growth abnormality | 3,262 | 0 | 0 | 1,221 | 0 | 0 | 757 | 1,066 | 0 | 218 | 0 |  
| HP:0004305 | Involuntary movements | 2,538 | 0 | 0 | 10 | 0 | 0 | 35 | 5 | 2,434 | 54 | 0 |  
| HP:0031466 | Impairment in personality functioning | 2,171 | 0 | 0 | 18 | 0 | 0 | 1,730 | 14 | 0 | 409 | 0 |  
| HP:0001249 | Intellectual disability | 2,072 | 0 | 0 | 129 | 814 | 0 | 671 | 158 | 0 | 300 | 0 |  
| HP:0000924 | Abnormality of the skeletal system | 1,986 | 0 | 0 | 338 | 0 | 0 | 0 | 1,292 | 0 | 356 | 0 |  
| HP:0012443 | Abnormality of brain morphology | 1,950 | 0 | 0 | 319 | 0 | 0 | 689 | 686 | 0 | 256 | 0 |  
| HP:0002197 | Generalized-onset seizure | 1,781 | 0 | 0 | 1,776 | 0 | 0 | 0 | 2 | 0 | 3 | 0 |  
| HP:0025031 | Abnormality of the digestive system | 1,588 | 0 | 0 | 91 | 0 | 0 | 1,258 | 146 | 0 | 93 | 0 |  
| HP:0003011 | Abnormality of the musculature | 1,376 | 0 | 0 | 242 | 0 | 0 | 0 | 705 | 0 | 429 | 0 |  
| HP:0001510 | Growth delay | 1,035 | 0 | 0 | 223 | 0 | 0 | 0 | 725 | 0 | 87 | 0 |  


--- 

## Case-control "metacohorts"  

To control for potential technical differences between cohorts, we combined CNV data from multiple cohorts into four matched groups for burden testing, dubbed **metacohorts**.  

These metacohorts represent the basic unit on which all burden testing was performed, and are described in the table below.  

For completeness, we also performed identical analyses on a pooled dataset of all samples, dubbed the **megacohort**.  

| Metacohort ID | Case Source(s) | Number of Cases | Control Sources(s) | Number of Controls |  
| :--- | :--- | ---: | :--- | ---: |  
| `meta1` | Coe, BCH, GDX | 42,643 | Coe, Cooper | 19,585 |  
| `meta2` | PGC, SSC, TSAICG | 26,323 | TCGA, PGC, TSAICG | 33,040 |  
| `meta3` | CHOP | 154,140 | CHOP | 28,070 |  
| `meta4` | UKBB | 0 (TBD) | UKBB | 480,501 |  
| `mega` | Coe, BCH, GDX, PGC, SSC, TSAICG, CHOP, UKBB | 223,116 | Coe, Cooper, TCGA, PGC, TSAICG, CHOP, UKBB | 561,196 |  

The count of samples for each phenotype group per metacohort is listed below:  

| HPO | description | Total | meta1 | meta2 | meta3 | meta4 | mega |  
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: |  
| HP:0000118 | Phenotypic abnormality | 223,116 | 42,653 | 26,323 | 154,140 | 0 | 223,116 |  
| HP:0000707 | Abnormality of the nervous system | 116,783 | 36,054 | 26,323 | 54,406 | 0 | 116,783 |  
| HP:0012638 | Abnormality of nervous system physiology | 91,067 | 35,796 | 26,323 | 28,948 | 0 | 91,067 |  
| HP:0000708 | Behavioral abnormality | 54,041 | 3,788 | 23,889 | 26,364 | 0 | 54,041 |  
| UNKNOWN | NA | 49,137 | 5,264 | 0 | 43,873 | 0 | 49,137 |  
| HP:0012759 | Neurodevelopmental abnormality | 36,378 | 34,893 | 814 | 671 | 0 | 36,378 |  
| HP:0012758 | Neurodevelopmental delay | 34,633 | 34,633 | 0 | 0 | 0 | 34,633 |  
| HP:0002715 | Abnormality of the immune system | 32,932 | 89 | 0 | 32,843 | 0 | 32,932 |  
| HP:0012639 | Abnormality of nervous system morphology | 26,277 | 1,341 | 0 | 24,936 | 0 | 26,277 |  
| HP:0007367 | Atrophy/Degeneration affecting the central nervous system | 24,257 | 10 | 0 | 24,247 | 0 | 24,257 |  
| HEALTHY_CONTROL | Unaffected control sample | 1,122,392 | 19,585 | 33,040 | 28,070 | 480,501 | 561,196 |  
| HP:0001626 | Abnormality of the cardiovascular system | 22,372 | 1,369 | 0 | 21,003 | 0 | 22,372 |  
| HP:0100753 | Schizophrenia | 22,301 | 7 | 21,094 | 1,200 | 0 | 22,301 |  
| HP:0000729 | Autistic behavior | 16,926 | 3,060 | 2,795 | 11,071 | 0 | 16,926 |  
| HP:0100022 | Abnormality of movement | 15,553 | 725 | 2,434 | 12,394 | 0 | 15,553 |  
| HP:0002597 | Abnormality of the vasculature | 15,073 | 125 | 0 | 14,948 | 0 | 15,073 |  
| HP:0000752 | Hyperactivity | 12,941 | 582 | 0 | 12,359 | 0 | 12,941 |  
| HP:0000736 | Short attention span | 12,908 | 549 | 0 | 12,359 | 0 | 12,908 |  
| HP:0000152 | Abnormality of head or neck | 7,045 | 7,045 | 0 | 0 | 0 | 7,045 |  
| HP:0001627 | Abnormal heart morphology | 6,751 | 696 | 0 | 6,055 | 0 | 6,751 |  
| HP:0001250 | Seizures | 4,965 | 2,885 | 318 | 1,762 | 0 | 4,965 |  
| HP:0001507 | Growth abnormality | 3,262 | 2,505 | 0 | 757 | 0 | 3,262 |  
| HP:0004305 | Involuntary movements | 2,538 | 69 | 2,434 | 35 | 0 | 2,538 |  
| HP:0031466 | Impairment in personality functioning | 2,171 | 441 | 0 | 1,730 | 0 | 2,171 |  
| HP:0001249 | Intellectual disability | 2,072 | 587 | 814 | 671 | 0 | 2,072 |  
| HP:0000924 | Abnormality of the skeletal system | 1,986 | 1,986 | 0 | 0 | 0 | 1,986 |  
| HP:0012443 | Abnormality of brain morphology | 1,950 | 1,261 | 0 | 689 | 0 | 1,950 |  
| HP:0002197 | Generalized-onset seizure | 1,781 | 1,781 | 0 | 0 | 0 | 1,781 |  
| HP:0025031 | Abnormality of the digestive system | 1,588 | 330 | 0 | 1,258 | 0 | 1,588 |  
| HP:0003011 | Abnormality of the musculature | 1,376 | 1,376 | 0 | 0 | 0 | 1,376 |  
| HP:0001510 | Growth delay | 1,035 | 1,035 | 0 | 0 | 0 | 1,035 |  

