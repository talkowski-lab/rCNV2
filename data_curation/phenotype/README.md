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
| HP:0000118 | Phenotypic abnormality | 223,116 | 1 | NA | HP:0007367, HP:0031466, HP:0001627, HP:0012759, HP:0000736, HP:0000752, HP:0012638, HP:0001626, HP:0001249, HP:0000152, HP:0004305, HP:0001510, HP:0012639, HP:0012443, HP:0000924, HP:0001507, HP:0002715, HP:0012758, HP:0002597, HP:0002197, HP:0100022, HP:0100753, HP:0000729, HP:0001250, HP:0025031, HP:0000708, HP:0003011, HP:0000707 |  
| HP:0000707 | Abnormality of the nervous system | 116,756 | 2 | HP:0000118 | HP:0012639, HP:0012443, HP:0001250, HP:0007367, HP:0012758, HP:0001249, HP:0000708, HP:0012759, HP:0000736, HP:0004305, HP:0000752, HP:0012638, HP:0002197, HP:0100022, HP:0100753, HP:0000729, HP:0031466 |  
| HP:0012638 | Abnormality of nervous system physiology | 91,042 | 3 | HP:0000707, HP:0000118 | HP:0000708, HP:0012759, HP:0000736, HP:0004305, HP:0000752, HP:0001250, HP:0012758, HP:0002197, HP:0100022, HP:0001249, HP:0100753, HP:0000729, HP:0031466 |  
| HP:0000708 | Behavioral abnormality | 54,033 | 4 | HP:0000707, HP:0012638, HP:0000118 | HP:0100753, HP:0000729, HP:0000736, HP:0031466 |  
| UNKNOWN | NA | 49,175 | 2 | HP:0000118 | NA |  
| HP:0012759 | Neurodevelopmental abnormality | 36,359 | 4 | HP:0000707, HP:0012638, HP:0000118 | HP:0001249, HP:0012758 |  
| HP:0012758 | Neurodevelopmental delay | 34,615 | 5 | HP:0000707, HP:0012638, HP:0012759, HP:0000118 | NA |  
| HP:0002715 | Abnormality of the immune system | 32,930 | 2 | HP:0000118 | NA |  
| HP:0012639 | Abnormality of nervous system morphology | 26,271 | 3 | HP:0000707, HP:0000118 | HP:0012443, HP:0007367 |  
| HP:0007367 | Atrophy/Degeneration affecting the central nervous system | 24,257 | 5 | HP:0000707, HP:0012639, HP:0000118 | NA |  
| HP:0001626 | Abnormality of the cardiovascular system | 22,367 | 2 | HP:0000118 | HP:0001627, HP:0002597 |  
| HP:0100753 | Schizophrenia | 22,301 | 5 | HP:0000707, HP:0012638, HP:0000708, HP:0000118 | NA |  
| HP:0000729 | Autistic behavior | 16,919 | 5 | HP:0000707, HP:0012638, HP:0000708, HP:0000118 | NA |  
| HP:0100022 | Abnormality of movement | 15,552 | 4 | HP:0000707, HP:0012638, HP:0000118 | HP:0004305, HP:0000752 |  
| HP:0002597 | Abnormality of the vasculature | 15,070 | 3 | HP:0000118, HP:0001626 | NA |  
| HP:0000752 | Hyperactivity | 12,940 | 5 | HP:0000707, HP:0012638, HP:0100022, HP:0000118 | NA |  
| HP:0000736 | Short attention span | 12,907 | 5 | HP:0000707, HP:0012638, HP:0000708, HP:0000118 | NA |  
| HP:0000152 | Abnormality of head or neck | 7,037 | 2 | HP:0000118 | NA |  
| HP:0001627 | Abnormal heart morphology | 6,749 | 4 | HP:0000118, HP:0001626 | NA |  
| HP:0001250 | Seizures | 4,961 | 4 | HP:0000707, HP:0012638, HP:0000118 | HP:0002197 |  
| HP:0001507 | Growth abnormality | 3,254 | 2 | HP:0000118 | HP:0001510 |  
| HP:0004305 | Involuntary movements | 2,538 | 5 | HP:0000707, HP:0012638, HP:0100022, HP:0000118 | NA |  
| HP:0031466 | Impairment in personality functioning | 2,171 | 5 | HP:0000707, HP:0012638, HP:0000708, HP:0000118 | NA |  
| HP:0001249 | Intellectual disability | 2,069 | 5 | HP:0000707, HP:0012638, HP:0012759, HP:0000118 | NA |  
| HP:0000924 | Abnormality of the skeletal system | 1,976 | 2 | HP:0000118 | NA |  
| HP:0012443 | Abnormality of brain morphology | 1,944 | 5 | HP:0000707, HP:0012639, HP:0000118 | NA |  
| HP:0002197 | Generalized-onset seizure | 1,781 | 5 | HP:0000707, HP:0012638, HP:0001250, HP:0000118 | NA |  
| HP:0025031 | Abnormality of the digestive system | 1,586 | 2 | HP:0000118 | NA |  
| HP:0003011 | Abnormality of the musculature | 1,373 | 2 | HP:0000118 | NA |  
| HP:0001510 | Growth delay | 1,030 | 3 | HP:0000118, HP:0001507 | NA |  

### HPO terms per cohort

The number of samples per phenotype group per cohort is outlined in the table, below.  

The code to generate this table is provided in `gather_hpo_per_cohort_table.py`.  

| HPO | description | Total | PGC | Cooper | Coe | SSC | UKBB | CHOP | GDX | TSAICG | BCH | TCGA |  
| :--- | :--- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |  
| HP:0000118 | Phenotypic abnormality | 223,116 | 21,094 | 0 | 29,104 | 2,795 | 0 | 154,140 | 9,958 | 2,434 | 3,591 | 0 |  
| HP:0000707 | Abnormality of the nervous system | 116,756 | 21,094 | 0 | 29,104 | 2,795 | 0 | 54,406 | 4,009 | 2,434 | 2,914 | 0 |  
| HP:0012638 | Abnormality of nervous system physiology | 91,042 | 21,094 | 0 | 29,104 | 2,795 | 0 | 28,948 | 3,821 | 2,434 | 2,846 | 0 |  
| HP:0000708 | Behavioral abnormality | 54,033 | 21,094 | 0 | 1,834 | 2,795 | 0 | 26,364 | 722 | 0 | 1,224 | 0 |  
| UNKNOWN | NA | 49,175 | 0 | 0 | 0 | 0 | 0 | 43,873 | 5,009 | 0 | 293 | 0 |  
| HP:0012759 | Neurodevelopmental abnormality | 36,359 | 0 | 0 | 29,104 | 814 | 0 | 671 | 3,404 | 0 | 2,366 | 0 |  
| HP:0012758 | Neurodevelopmental delay | 34,615 | 0 | 0 | 29,104 | 0 | 0 | 0 | 3,302 | 0 | 2,209 | 0 |  
| HP:0002715 | Abnormality of the immune system | 32,930 | 0 | 0 | 19 | 0 | 0 | 32,843 | 34 | 0 | 34 | 0 |  
| HP:0012639 | Abnormality of nervous system morphology | 26,271 | 0 | 0 | 352 | 0 | 0 | 24,936 | 706 | 0 | 277 | 0 |  
| HP:0007367 | Atrophy/Degeneration affecting the central nervous system | 24,257 | 0 | 0 | 4 | 0 | 0 | 24,247 | 4 | 0 | 2 | 0 |  
| HEALTHY_CONTROL | Unaffected control sample | 561,196 | 20,277 | 8,329 | 11,256 | 0 | 480,501 | 28,070 | 0 | 4,093 | 0 | 8,670 |  
| HP:0001626 | Abnormality of the cardiovascular system | 22,367 | 0 | 0 | 748 | 0 | 0 | 21,003 | 404 | 0 | 212 | 0 |  
| HP:0100753 | Schizophrenia | 22,301 | 21,094 | 0 | 2 | 0 | 0 | 1,200 | 5 | 0 | 0 | 0 |  
| HP:0000729 | Autistic behavior | 16,919 | 0 | 0 | 1,569 | 2,795 | 0 | 11,071 | 631 | 0 | 853 | 0 |  
| HP:0100022 | Abnormality of movement | 15,552 | 0 | 0 | 317 | 0 | 0 | 12,394 | 82 | 2,434 | 325 | 0 |  
| HP:0002597 | Abnormality of the vasculature | 15,070 | 0 | 0 | 50 | 0 | 0 | 14,948 | 56 | 0 | 16 | 0 |  
| HP:0000752 | Hyperactivity | 12,940 | 0 | 0 | 281 | 0 | 0 | 12,359 | 62 | 0 | 238 | 0 |  
| HP:0000736 | Short attention span | 12,907 | 0 | 0 | 264 | 0 | 0 | 12,359 | 57 | 0 | 227 | 0 |  
| HP:0000152 | Abnormality of head or neck | 7,037 | 0 | 0 | 4,535 | 0 | 0 | 0 | 1,654 | 0 | 848 | 0 |  
| HP:0001627 | Abnormal heart morphology | 6,749 | 0 | 0 | 166 | 0 | 0 | 6,055 | 332 | 0 | 196 | 0 |  
| HP:0001250 | Seizures | 4,961 | 0 | 0 | 1,807 | 318 | 0 | 1,762 | 535 | 0 | 539 | 0 |  
| HP:0001507 | Growth abnormality | 3,254 | 0 | 0 | 1,221 | 0 | 0 | 757 | 1,058 | 0 | 218 | 0 |  
| HP:0004305 | Involuntary movements | 2,538 | 0 | 0 | 10 | 0 | 0 | 35 | 5 | 2,434 | 54 | 0 |  
| HP:0031466 | Impairment in personality functioning | 2,171 | 0 | 0 | 18 | 0 | 0 | 1,730 | 14 | 0 | 409 | 0 |  
| HP:0001249 | Intellectual disability | 2,069 | 0 | 0 | 129 | 814 | 0 | 671 | 155 | 0 | 300 | 0 |  
| HP:0000924 | Abnormality of the skeletal system | 1,976 | 0 | 0 | 338 | 0 | 0 | 0 | 1,282 | 0 | 356 | 0 |  
| HP:0012443 | Abnormality of brain morphology | 1,944 | 0 | 0 | 319 | 0 | 0 | 689 | 680 | 0 | 256 | 0 |  
| HP:0002197 | Generalized-onset seizure | 1,781 | 0 | 0 | 1,776 | 0 | 0 | 0 | 2 | 0 | 3 | 0 |  
| HP:0025031 | Abnormality of the digestive system | 1,586 | 0 | 0 | 91 | 0 | 0 | 1,258 | 144 | 0 | 93 | 0 |  
| HP:0003011 | Abnormality of the musculature | 1,373 | 0 | 0 | 242 | 0 | 0 | 0 | 702 | 0 | 429 | 0 |  
| HP:0001510 | Growth delay | 1,030 | 0 | 0 | 223 | 0 | 0 | 0 | 720 | 0 | 87 | 0 |  


--- 

## Case-control "metacohorts"  

To further control for potential technical differences between cohorts, we broke down available case and control datasets into four matched groups for burden testing.  

These metacohorts represent the base unit on which all burden testing was performed.  

These case-control pairings (dubbed **metacohorts**) are described below.  

| Metacohort ID | Case Source(s) | Number of Cases | Control Sources(s) | Number of Controls |  
| :--- | :--- | ---: | :--- | ---: |  
| `meta1` | Coe, BCH, GDX | 42,643 | Coe, Cooper | 19,585 |  
| `meta2` | PGC, SSC, TSAICG | 26,323 | TCGA, PGC, TSAICG | 33,040 |  
| `meta3` | CHOP | 154,140 | CHOP | 28,070 |  
| `meta4` | UKBB | 0 (TBD) | UKBB | 480,501 |  

The count of samples for each phenotype group per metacohort is listed below:  

| HPO | description | Total | meta1 | meta2 | meta3 | meta4 |  
| :--- | :--- | ---: | ---: | ---: | ---: | ---: |  
| HP:0000118 | Phenotypic abnormality | 223,116 | 42,653 | 26,323 | 154,140 | 0 |  
| HP:0000707 | Abnormality of the nervous system | 116,756 | 36,027 | 26,323 | 54,406 | 0 |  
| HP:0012638 | Abnormality of nervous system physiology | 91,042 | 35,771 | 26,323 | 28,948 | 0 |  
| HP:0000708 | Behavioral abnormality | 54,033 | 3,780 | 23,889 | 26,364 | 0 |  
| UNKNOWN | NA | 49,175 | 5,302 | 0 | 43,873 | 0 |  
| HP:0012759 | Neurodevelopmental abnormality | 36,359 | 34,874 | 814 | 671 | 0 |  
| HP:0012758 | Neurodevelopmental delay | 34,615 | 34,615 | 0 | 0 | 0 |  
| HP:0002715 | Abnormality of the immune system | 32,930 | 87 | 0 | 32,843 | 0 |  
| HP:0012639 | Abnormality of nervous system morphology | 26,271 | 1,335 | 0 | 24,936 | 0 |  
| HP:0007367 | Atrophy/Degeneration affecting the central nervous system | 24,257 | 10 | 0 | 24,247 | 0 |  
| HEALTHY_CONTROL | Unaffected control sample | 561,196 | 19,585 | 33,040 | 28,070 | 480,501 |  
| HP:0001626 | Abnormality of the cardiovascular system | 22,367 | 1,364 | 0 | 21,003 | 0 |  
| HP:0100753 | Schizophrenia | 22,301 | 7 | 21,094 | 1,200 | 0 |  
| HP:0000729 | Autistic behavior | 16,919 | 3,053 | 2,795 | 11,071 | 0 |  
| HP:0100022 | Abnormality of movement | 15,552 | 724 | 2,434 | 12,394 | 0 |  
| HP:0002597 | Abnormality of the vasculature | 15,070 | 122 | 0 | 14,948 | 0 |  
| HP:0000752 | Hyperactivity | 12,940 | 581 | 0 | 12,359 | 0 |  
| HP:0000736 | Short attention span | 12,907 | 548 | 0 | 12,359 | 0 |  
| HP:0000152 | Abnormality of head or neck | 7,037 | 7,037 | 0 | 0 | 0 |  
| HP:0001627 | Abnormal heart morphology | 6,749 | 694 | 0 | 6,055 | 0 |  
| HP:0001250 | Seizures | 4,961 | 2,881 | 318 | 1,762 | 0 |  
| HP:0001507 | Growth abnormality | 3,254 | 2,497 | 0 | 757 | 0 |  
| HP:0004305 | Involuntary movements | 2,538 | 69 | 2,434 | 35 | 0 |  
| HP:0031466 | Impairment in personality functioning | 2,171 | 441 | 0 | 1,730 | 0 |  
| HP:0001249 | Intellectual disability | 2,069 | 584 | 814 | 671 | 0 |  
| HP:0000924 | Abnormality of the skeletal system | 1,976 | 1,976 | 0 | 0 | 0 |  
| HP:0012443 | Abnormality of brain morphology | 1,944 | 1,255 | 0 | 689 | 0 |  
| HP:0002197 | Generalized-onset seizure | 1,781 | 1,781 | 0 | 0 | 0 |  
| HP:0025031 | Abnormality of the digestive system | 1,586 | 328 | 0 | 1,258 | 0 |  
| HP:0003011 | Abnormality of the musculature | 1,373 | 1,373 | 0 | 0 | 0 |  
| HP:0001510 | Growth delay | 1,030 | 1,030 | 0 | 0 | 0 |  
