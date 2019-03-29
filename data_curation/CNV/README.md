# CNV Data Curation

## CNV Data Curation  

We aggregated existing microarray-based CNV calls from numerous sources. This process is described below.  

All commands executed to filter the CNV data as descrbed below are contained in `filter_CNV_data.sh`.  

### Glossary & abbreviations  

Several abbreviations and acronyms are used below, as follows:  

| Abbreviation | Definition |
| --- | :--- |
| CNV | Copy Number Variant |
| CTRL | Control (Unaffected Sample) |
| SCZ | Schizophrenia |
| DD | Developmental Disorders |
| ASD | Autism Spectrum Disorders |
| TS | Tourette Syndrome |

### CNV data sources  

We aggregated CNV data from multiple sources, listed below:  

| Dataset | Citation | PMID | Platform(s) | Build | Phenos | N Cases | N Ctrls |
| --- | :--- | :--- | :--- | :--- | :--- | ---: | ---: |
| PGC | [Marshall _et al._, _Nat. Genet._ (2017)](https://www.nature.com/articles/ng.3725) | [27869829](https://www.ncbi.nlm.nih.gov/pubmed/27869829) | Affy 6.0 (37%), Omni Express (31%), Omni Express Plus (12%), Other (20%) | hg18 | SCZ | 21,094 | 20,277 |
| Coe | [Coe _et al._, _Nat. Genet._ (2014)](https://www.nature.com/articles/ng.3092) | [25217958](https://www.ncbi.nlm.nih.gov/pubmed/25217958) | Cases: SignatureChip OS v2.0 (58%), SignatureChip OS v1.0 (34%), Other (8%); Controls: various low-density SNP genotyping arrays | hg19 | DD | 29,085 | 19,584 |
| Cooper<sup>1</sup> | [Cooper _et al._, _Nat. Genet._ (2011)](https://www.nature.com/articles/ng.909) | [21841781](https://www.ncbi.nlm.nih.gov/pubmed/21841781) | Various low-density SNP arrays | hg19 | DD | 0<sup>1</sup> | 8,329 |
| SSC<sup>2</sup> | [Sanders _et at._, _Neuron_ (2015)](https://www.sciencedirect.com/science/article/pii/S0896627315007734?) | [26402605](https://www.ncbi.nlm.nih.gov/pubmed/26402605) | Omni 1Mv3 (46%), Omni 2.5 (41%), Omni 1Mv1 (13%) | hg18 | ASD | 2,795 | 0<sup>2</sup> |
| UKBB | [Macé _et al._, _Nat. Comms._ (2017)](https://www.nature.com/articles/s41467-017-00556-x) | [28963451](https://www.ncbi.nlm.nih.gov/pubmed/28963451) | TBD | TBD | Mixed | TBD | TBD |
| CHOP | - | - | TBD | TBD | TBD | TBD | TBD |
| GeneDX | - | - | TBD | hg18 & hg19 | Mixed | TBD | TBD |
| TSAICG | [Huang _et al._, _Neuron_ (2017)](https://www.sciencedirect.com/science/article/pii/S0896627317305081) | [28641109](https://www.ncbi.nlm.nih.gov/pubmed/28641109) | OmniExpress | TBD | TS | 2,434 | 4,093 |
| BCH | [Talkowski _et al._, _Cell_ (2012)](https://www.sciencedirect.com/science/article/pii/S0092867412004114) | [22521361](https://www.ncbi.nlm.nih.gov/pubmed/22521361) | TBD | TBD | TBD | TBD | TBD |  

#### Notes on raw CNV data   
1. Only retained control samples from Cooper _et al._. All cases from Cooper _et al._ also appear in Coe _et al._.  
2. Only retained affected children from Sanders _et al._, since all controls were first-degree relatives of affected cases.  

### Raw CNV data processing steps  

All CNV data native to hg18 was lifted over to hg19 coordinate space using UCSC liftOver, requiring at least 50% of the original CNV to map to hg19 successfully in order to be retained.  

Some datasets required manual curation prior to inclusion. Where necessary, these steps are enumerated below:  

 * **SSC**: CNVs were filtered on pCNV ≤ 10<sup>-9</sup>, per recommendation of the authors.  

### Raw CNV callset properties  

Prior to filtering, the datasets were as follows:  

| Dataset | N Cases | Case CNVs | Case Median Size | Case DEL:DUP | CNVs /Case | N Ctrls | Ctrl CNVs | Ctrl Median Size | Ctrl DEL:DUP | CNVs /Ctrl | 
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| PGC | 21,094 | 42,096 | 79.7 kb | 1:1.04 | 2.0 | 20,277 | 40,464 | 78.2 kb | 1:1.03 | 2.0 |
| Coe | 29,085 | TBD | TBD | TBD | TBD | 19,584 | TBD | TBD | TBD | TBD |
| Cooper | 0 | 0 | - | - |  - | 8,329 | TBD | TBD | TBD | TBD |
| SSC | 2,795 | 30,867 | 21.0 kb | 3.09:1 | 11.0 | 0 | 0 | - | - | - |
| UKBB | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| CHOP | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| GeneDX | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| TSAICG | 2,434 | TBD | TBD | TBD | TBD | 4,093 | TBD | TBD | TBD | TBD |
| BCH | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |

### Raw data access  

All raw CNV data is stored in a protected Google Cloud bucket, here:  
```
$ gsutil ls gs://rcnv_project/raw_data/cnv/

gs://rcnv_project/raw_data/cnv/PGC.raw.bed.gz
gs://rcnv_project/raw_data/cnv/SSC.raw.bed.gz
```

Note that permissions must be granted per user prior to data access.  

## Curation Steps: Rare CNVs  

All raw CNV data was subjected to the same set of global filters:  
1. Restricted to autosomes
2. Does not have substantial overlap<sup>1</sup> a common (AF≥1%) CNV from WGS resolution in gnomAD-SV<sup>2</sup> ([Collins\*, Brand\*, _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/biorxiv/early/2019/03/14/578674))  
3. Does not have substantial overlap<sup>1</sup> with other CNVs in at least 1% of all samples within the same dataset
4. Does not have substantial overlap<sup>1</sup> with other CNVs in at least 1% of all samples across all datasets
5. CNV size ≥ 50kb and ≤ 10Mb
6. Not substantially covered<sup>3</sup> by somatically hypermutable sites (described in [Collins\*, Brand\*, _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/biorxiv/early/2019/03/14/578674))  
7. Not substantially covered<sup>3</sup> by segmental duplications and/or simple repeats  
8. Not substantially covered<sup>3</sup> by N-masked regions of the hg19 reference genome assembly  

#### Notes on curation  
1. "Substantial" overlap determined based on ≥50% reciprocal overlap using BEDTools ([Quinlan & Hall, _Bioinformatics_ (2010)](https://www.ncbi.nlm.nih.gov/pubmed/20110278)). For overlap-based comparisons, both breakpoints were required to be within ±50kb, and CNV type (DEL vs. DUP) was required to match.  
2. The version of gnomAD-SV used for this analysis included 12,549 samples as reported in [the gnomAD-SV preprint](https://www.biorxiv.org/content/biorxiv/early/2019/03/14/578674), and not the smaller release subset of 10,738 genomes [available from the gnomAD website](https://gnomad.broadinstitute.org/downloads/)
3. "Substantial" coverage determined based on ≥30% coverage per BEDTools coverage ([Quinlan & Hall, _Bioinformatics_ (2010)](https://www.ncbi.nlm.nih.gov/pubmed/20110278)).

### Rare CNV callset properties

After the filtering steps described above, the datasets were as follows:

| Dataset | N Cases | Case CNVs | Case Median Size | Case DEL:DUP | CNVs /Case | N Ctrls | Ctrl CNVs | Ctrl Median Size | Ctrl DEL:DUP | CNVs /Ctrl | 
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| PGC | 21,094 | TBD | TBD | TBD | TBD | 20,277 | TBD | TBD | TBD | TBD |
| Coe | 29,085 | TBD | TBD | TBD | TBD | 19,584 | TBD | TBD | TBD | TBD |
| Cooper | 0 | 0 | - | - |  - | 8,329 | TBD | TBD | TBD | TBD |
| SSC | TBD | TBD | TBD | TBD | TBD | 0 | 0 | - | - | - |
| UKBB | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| CHOP | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| GeneDX | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| TSAICG | 2,434 | TBD | TBD | TBD | TBD | 4,093 | TBD | TBD | TBD | TBD |
| BCH | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |


## Curation Steps: Ultra-rare CNVs  

In parallel to the creation of the rare CNV dataset (described above), a dataset of ultra-rare CNVs was also generated.

These ultra-rare CNVs were subjected to the same set of filters as above, with the following modifications:  
2. Does not have substantial overlap **an ultra-rare (AF≥0.01%) CNV** from WGS resolution in gnomAD-SV ([Collins\*, Brand\*, _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/biorxiv/early/2019/03/14/578674))  
3. Does not have substantial overlap with other CNVs in **at least 0.01%** of all samples within the same dataset (to a minimum of one sample)  
4. Does not have substantial overlap with other CNVs in **at least 0.01%** of all samples across all datasets (to a minimum of one sample)  

### Ultra-rare CNV callset properties

After the filtering steps described above, the datasets were as follows:

| Dataset | N Cases | Case CNVs | Case Median Size | Case DEL:DUP | CNVs /Case | N Ctrls | Ctrl CNVs | Ctrl Median Size | Ctrl DEL:DUP | CNVs /Ctrl | 
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| PGC | 21,094 | TBD | TBD | TBD | TBD | 20,277 | TBD | TBD | TBD | TBD |
| Coe | 29,085 | TBD | TBD | TBD | TBD | 19,584 | TBD | TBD | TBD | TBD |
| Cooper | 0 | 0 | - | - |  - | 8,329 | TBD | TBD | TBD | TBD |
| SSC | TBD | TBD | TBD | TBD | TBD | 0 | 0 | - | - | - |
| UKBB | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| CHOP | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| GeneDX | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
| TSAICG | 2,434 | TBD | TBD | TBD | TBD | 4,093 | TBD | TBD | TBD | TBD |
| BCH | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD | TBD |
