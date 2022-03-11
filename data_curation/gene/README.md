# Gene Curation  

We curated and annotated all protein-coding genes in the genome for association testing. This process is described below.  

All commands executed to filter the CNV data as descrbed below are contained in `format_genes.sh`.    

## Gene definitions  

We restricted all gene-based analyses to canonical transcripts from autosomal protein-coding genes as defined by [Gencode v19](https://www.gencodegenes.org/human/release_19.html).  

The Gencode definition of "canonical transcript" is [provided here](http://www.ensembl.org/Help/Glossary?id=346).  

We extracted all canonical transcripts from protein-coding genes using `get_canonical_transcripts.py`.   

## Exon-level expression filtering  

We also restricted all gene-based analyses to exons expressed in at least 20% of transcripts in at least one human tissue catalogued by the [Genotype-Tissue Expression (GTEx) Project](https://gtexportal.org/).  

For this filter, we computed the per-exon max `pext` score as defined in [Cummings _et al._, _bioRxiv_ (2019)](https://www.biorxiv.org/content/10.1101/554444v1) across all available tissues.  

Bases from exons with missing `pext` scores were ignored when calculating mean `pext` scores, and exons entirely missing `pext` scores were considered as passing.  

Genes with no exons remaining after exon-level expression filtering were removed outright from all subsequent analyses.  

## Gene features  

For certain analyses, we used per-gene features across a variety of categories.  

These "gene features" are described below.  

For each feature set, we controlled for inter-feature correlation structure by decomposing all features with principal components analysis and retaining the top "eigenfeatures" that explained at least 99% of inter-gene variance.  

In practice, the feature collection and decomposition process was parallelized in [FireCloud/Terra](https://portal.firecloud.org) with `get_gene_metadata.wdl`.  

#### Gene features: genomic  

We collected the following genome-based gene features:  

| Feature name | Abbreviation | Description | Source |  
| :--- | :--- | :--- | :--- |  
| Gene length | `gene_length` | Length of gene body corresponding to canonical transcript | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) |  
| Dist. to nearest gene | `nearest_gene` | Distance to nearest other gene. Defaults to zero for overlapping genes. | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | 
| Genes within 1Mb | `genes_within_1mb` | Count of all other genes within ±1Mb of gene body | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) |  
| Length of coding sequence | `cds_length` | Length of coding sequence | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | 
| Number of exons | `n_exons` | Number of non-overlapping exons | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | 
| Min. exon size | `min_exon_size` | Size of smallest exon | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | 
| Median exon size | `med_exon_size` | Median exon size | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | 
| Max. exon size | `max_exon_size` | Size of largest exon | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | 
| Min. intron size | `min_intron_size` | Size of smallest intron. Defauts to zero for single-exon genes. | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | 
| Median intron size | `med_intron_size` | Median intron size. Defauts to zero for single-exon genes. | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | 
| Max. intron size | `max_intron_size` | Size of largest intron. Defauts to zero for single-exon genes. | Gencode v19 [(Harrow _et al._, _Genome Res._, 2012)](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | 
| Dist. to telomere | `telomere_dist` | Absolute distance to nearest telomere | [UCSC Genome Browser](http://genome.ucsc.edu) |    
| Dist. to centromere | `centromere_dist` | Absolute distance to nearest centromere | [UCSC Genome Browser](http://genome.ucsc.edu) |  
| Norm. dist. to telomere | `telomere_dist_norm` | Distance to nearest telomere, normalized to chromosome arm length | [UCSC Genome Browser](http://genome.ucsc.edu) |    
| Norm. dist. to centromere | `centromere_dist_norm` | Distance to nearest centromere, normalized to chromosome arm length | [UCSC Genome Browser](http://genome.ucsc.edu) |  
| GC content | `gc_pct` | Fraction of G/C nucleotides within gene body | [UCSC Genome Browser](http://genome.ucsc.edu) |  
| Repetitive sequence content | `repeat_cov` | Fraction of gene body covered by segmental duplications and/or simple/low-complexity/satellite repeats | [UCSC Genome Browser](http://genome.ucsc.edu) |  
| Sequence alignability | `align_100mer` | Mean 100mer alignability score over the gene body | [UCSC Genome Browser](http://genome.ucsc.edu) |  
| Synonymous mutation rate | `gnomad_mu_syn` | Expected rate of synonymous mutations | gnomAD v2.1 [(Karczewski _et al._, _Nature_, 2020)](https://pubmed.ncbi.nlm.nih.gov/32461654) |  
| Missense mutation rate | `gnomad_mu_mis` | Expected rate of missense mutations | gnomAD v2.1 [(Karczewski _et al._, _Nature_, 2020)](https://pubmed.ncbi.nlm.nih.gov/32461654) |  
| Loss-of-function mutation rate | `gnomad_mu_lof` | Expected rate of loss-of-function mutations | gnomAD v2.1 [(Karczewski _et al._, _Nature_, 2020)](https://pubmed.ncbi.nlm.nih.gov/32461654) |  

#### Gene features: expression  

We collected the following expression-based gene features:  

| Feature name | Abbreviation | Description | Source |  
| :--- | :--- | :--- | :--- |  
| Number of expressing tissues | `n_tissues_expressed` | Number of tissues with median expression (TPM) > 1 in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Min. average expression | `median_expression_min` | Minimum of median expression values across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| First quartile of average expression | `median_expression_q1` | First quartile of median expression values across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Mean average expression | `median_expression_mean` | Mean of median expression values across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Third quartile of average expression | `median_expression_q3` | Third quartile of median expression values across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Max. average expression | `median_expression_max` | Maximum of median expression values across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Std. dev of average expression | `median_expression_sd` | Standard deviation of median expression values across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Min. expression MAD | `expression_mad_min` | Minimum of median absolute deviation (MAD) of expression values  across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| First quartile of expression MAD | `expression_mad_q1` | First quartile of median absolute deviation (MAD) of expression values across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Mean expression MAD | `expression_mad_mean` | Mean of median absolute deviation (MAD) of expression values  across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Third quartile of expression MAD | `expression_mad_q3` | Third quartile of median absolute deviation (MAD) of expression values  across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Max. expression MAD | `expression_mad_max` | Maximum of median absolute deviation (MAD) of expression values  across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Expression MAD std. dev. | `expression_mad_sd` | Standard deviation of median absolute deviation (MAD) of expression values across all tissues in GTEx v7 | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  
| Expression components 1-20 | `expression_component_1` ... `expression_component_20` | Top 20 principal components of gene X tissue matrix from GTEx v7, using mean expression across individuals per gene per tissue | GTEx v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) |  

All expression values were derived from GTEx v7 (rather than the final v8 release) because GTEx v7 was the last version native to hg19/Gencode v19, which was a direct match for the Gencode version used in this study.  

Note that all expression values from GTEx were transformed as log<sub>10</sub>(x + 1).  

#### Gene features: constraint

We collected the following gene features related to evolutionary constraint:  

| Feature name | Abbreviation | Description | Source |  
| :--- | :--- | :--- | :--- |  
| Probability of loss-of-function intolerance | `gnomad_pLI` | Probability of loss-of-function intolerance calculated in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _Nature_, 2020)](https://pubmed.ncbi.nlm.nih.gov/32461654) |  
| Probability of complete haplosufficiency | `gnomad_pNull` | Probability of haplosufficiency calculated in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _Nature_, 2020)](https://pubmed.ncbi.nlm.nih.gov/32461654) |  
| Probability of recessive lethality | `gnomad_pRec` | Probability of recessive lethality calculated in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _Nature_, 2020)](https://pubmed.ncbi.nlm.nih.gov/32461654) |  
| gnomAD missense obs:exp | `gnomad_oe_mis` | Observed : expected ratio for missense SNVs in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _Nature_, 2020)](https://pubmed.ncbi.nlm.nih.gov/32461654) |  
| gnomAD loss-of-function obs:exp | `gnomad_oe_lof` | Observed : expected ratio for loss-of-function SNVs in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _Nature_, 2020)](https://pubmed.ncbi.nlm.nih.gov/32461654) |  
| gnomAD missense obs:exp (upper 90% CI) | `gnomad_oe_mis_upper` | Upper 90% confidence interval of obs:exp ratio for missense SNVs in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _Nature_, 2020)](https://pubmed.ncbi.nlm.nih.gov/32461654) |  
| gnomAD loss-of-function obs:exp (upper 90% CI) | `gnomad_oe_lof_upper` | Upper 90% confidence interval of obs:exp ratio for loss-of-function SNVs in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _Nature_, 2020)](https://pubmed.ncbi.nlm.nih.gov/32461654) |  
| gnomAD missense Z-score | `gnomad_mis_z` | Z-score for missense obs:exp ratio in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _Nature_, 2020)](https://pubmed.ncbi.nlm.nih.gov/32461654) |  
| gnomAD loss-of-function Z-score | `gnomad_lof_z` | Z-score for loss-of-function obs:exp ratio in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _Nature_, 2020)](https://pubmed.ncbi.nlm.nih.gov/32461654) |  
| ExAC deletion intolerance Z-score | `exac_del_z` | Z-score for obs:exp ratio of coding deletions in ExAC | ExAC v1.0 [(Ruderfer _et al._, _Nat. Genet._, 2016)](https://www.ncbi.nlm.nih.gov/pubmed/27533299) |  
| ExAC duplication intolerance Z-score | `exac_dup_z` | Z-score for obs:exp ratio of coding duplications in ExAC | ExAC v1.0 [(Ruderfer _et al._, _Nat. Genet._, 2016)](https://www.ncbi.nlm.nih.gov/pubmed/27533299) |  
| ExAC CNV intolerance Z-score | `exac_cnv_z` | Z-score for obs:exp ratio of coding CNVs in ExAC | ExAC v1.0 [(Ruderfer _et al._, _Nat. Genet._, 2016)](https://www.ncbi.nlm.nih.gov/pubmed/27533299) |  
| Haploinsufficiency score | `hurles_hi` | Haploinsufficiency scores predicted based on haplosufficient genes | [Huang _et al._, _PLOS Genetics_, 2010](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1001154) |  
| RVIS | `rvis` | Residual variation intolerance score (RVIS) computed on gnomAD v2.0 | [Petrovski _et al._, _PLOS Genetics_, 2013](https://www.ncbi.nlm.nih.gov/pubmed/23990802) |  
| RVIS percentile | `rvis_pct` | RVIS percentile computed on gnomAD v2.0 | [Petrovski _et al._, _PLOS Genetics_, 2013](https://www.ncbi.nlm.nih.gov/pubmed/23990802) |  
| sHet | `sHet` | Estimated selection coefficient on heterozygous loss-of-function mutations in ExAC | [Cassa _et al._, _Nat. Genet., 2017](https://pubmed.ncbi.nlm.nih.gov/28369035) |  
| Enhancer Domain Score | `EDS` | Predicted probability of being a disease gene based on enhancer domain features | [Wang and Goldstein, _Am. J. Hum. Genet., 2020](https://pubmed.ncbi.nlm.nih.gov/32032514/) |  
| CCDG deletion sensitivity score | `ccdg_del` | Deletion sensitivity score estimated from population CNV data in 17,795 individuals | [Abel _et al._, _Nature_, 2020](https://pubmed.ncbi.nlm.nih.gov/32460305) |  
| CCDG duplication sensitivity score | `ccdg_dup` | Duplication sensitivity score estimated from population CNV data in 17,795 individuals | [Abel _et al._, _Nature_, 2020](https://pubmed.ncbi.nlm.nih.gov/32460305) |  
| Promoter GC content | `promoter_gc_pct` | Percentage of C/G nucleotides in gene promoter | [UCSC Genome Browser](http://genome.ucsc.edu) |  
| Promoter CpG count | `promoter_cpg_count` | Count of CpG dinucleotides in gene promoter | [UCSC Genome Browser](http://genome.ucsc.edu) |  
| Promoter conservation | `promoter_phastcons` | Mean 100-way phastCons score in gene promoter | [UCSC Genome Browser](http://genome.ucsc.edu) |  
| Exon conservation | `exon_phastcons` | Mean 100-way hastCons score in exons (weighted by exon size) | [UCSC Genome Browser](http://genome.ucsc.edu) |  
| Gene body conservation | `gene_body_phastcons` | Mean 100-way phastCons score across entire gene body | [UCSC Genome Browser](http://genome.ucsc.edu) |  

#### Gene features: chromatin  

We collected the following chromatin-based gene features:  

| Feature name | Abbreviation | Description | Source |  
| :--- | :--- | :--- | :--- |  
| Mean active TSS coverage | chromhmm_1_TssA_mean | Mean gene coverage by active TSS from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Standard deviation of active TSS coverage | chromhmm_1_TssA_sd | Standard deviation of gene coverage by active TSS from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Mean flanking TSS coverage | chromhmm_2_TssFlnk_mean | Mean gene coverage by flanking TSS from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Standard deviation of flanking TSS coverage | chromhmm_2_TssFlnk_sd | Standard deviation of gene coverage by flanking TSS from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Mean flanking TSS upstream coverage | chromhmm_3_TssFlnkU_mean | Mean gene coverage by flanking TSS upstream from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Standard deviation of flanking TSS upstream coverage | chromhmm_3_TssFlnkU_sd | Standard deviation of gene coverage by flanking TSS upstream from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Mean flanking TSS downstream coverage | chromhmm_4_TssFlnkD_mean | Mean gene coverage by flanking TSS downstream from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Standard deviation of flanking TSS downstream coverage | chromhmm_4_TssFlnkD_sd | Standard deviation of gene coverage by flanking TSS downstream from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Mean strong transcription coverage | chromhmm_5_Tx_mean | Mean gene coverage by strong transcription from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Standard deviation of strong transcription coverage | chromhmm_5_Tx_sd | Standard deviation of gene coverage by strong transcription from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Mean weak transcription coverage | chromhmm_6_TxWk_mean | Mean gene coverage by weak transcription from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Standard deviation of weak transcription coverage | chromhmm_6_TxWk_sd | Standard deviation of gene coverage by weak transcription from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Mean genic enhancer1 coverage | chromhmm_7_EnhG1_mean | Mean gene coverage by genic enhancer1 from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Standard deviation of genic enhancer1 coverage | chromhmm_7_EnhG1_sd | Standard deviation of gene coverage by genic enhancer1 from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Mean genic enhancer2 coverage | chromhmm_8_EnhG2_mean | Mean gene coverage by genic enhancer2 from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Standard deviation of genic enhancer2 coverage | chromhmm_8_EnhG2_sd | Standard deviation of gene coverage by genic enhancer2 from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Mean active enhancer 1 coverage | chromhmm_9_EnhA1_mean | Mean gene coverage by active enhancer 1 from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Standard deviation of active enhancer 1 coverage | chromhmm_9_EnhA1_sd | Standard deviation of gene coverage by active enhancer 1 from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Mean active enhancer 2 coverage | chromhmm_10_EnhA2_mean | Mean gene coverage by active enhancer 2 from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Standard deviation of active enhancer 2 coverage | chromhmm_10_EnhA2_sd | Standard deviation of gene coverage by active enhancer 2 from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Mean weak enhancer coverage | chromhmm_11_EnhWk_mean | Mean gene coverage by weak enhancer from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Standard deviation of weak enhancer coverage | chromhmm_11_EnhWk_sd | Standard deviation of gene coverage by weak enhancer from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Mean znf genes & repeats coverage | chromhmm_12_ZNF/Rpts_mean | Mean gene coverage by znf genes & repeats from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Standard deviation of znf genes & repeats coverage | chromhmm_12_ZNF/Rpts_sd | Standard deviation of gene coverage by znf genes & repeats from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Mean heterochromatin coverage | chromhmm_13_Het_mean | Mean gene coverage by heterochromatin from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Standard deviation of heterochromatin coverage | chromhmm_13_Het_sd | Standard deviation of gene coverage by heterochromatin from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Mean bivalent/poised TSS coverage | chromhmm_14_TssBiv_mean | Mean gene coverage by bivalent/poised TSS from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Standard deviation of bivalent/poised TSS coverage | chromhmm_14_TssBiv_sd | Standard deviation of gene coverage by bivalent/poised TSS from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Mean bivalent enhancer coverage | chromhmm_15_EnhBiv_mean | Mean gene coverage by bivalent enhancer from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Standard deviation of bivalent enhancer coverage | chromhmm_15_EnhBiv_sd | Standard deviation of gene coverage by bivalent enhancer from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Mean repressed polycomb coverage | chromhmm_16_ReprPC_mean | Mean gene coverage by repressed polycomb from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Standard deviation of repressed polycomb coverage | chromhmm_16_ReprPC_sd | Standard deviation of gene coverage by repressed polycomb from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Mean weak repressed polycomb coverage | chromhmm_17_ReprPCWk_mean | Mean gene coverage by weak repressed polycomb from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Standard deviation of weak repressed polycomb coverage | chromhmm_17_ReprPCWk_sd | Standard deviation of gene coverage by weak repressed polycomb from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Mean quiescent/low coverage | chromhmm_18_Quies_mean | Mean gene coverage by quiescent/low from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Standard deviation of quiescent/low coverage | chromhmm_18_Quies_sd | Standard deviation of gene coverage by quiescent/low from ChromHMM across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Chromatin components 1-20 | `chromatin_component_1` ... `chromatin_component_20` | Top 20 principal components of gene X ChromHMM state matrix across 98 tissues from the Roadmap Epigenomics Project | Roadmap Epigenomics Project [(Kundaje _et al._, _Nature_, 2015)](https://www.nature.com/articles/nature14248) |  
| Episcore | `episcore` | A predictive score for haploinsufficient genes based on their epigenetic signatures across tissues | [Han _et al._, _Nat. Commun._, 2018](https://pubmed.ncbi.nlm.nih.gov/29849042) |  

All chromatin data was based on the [Roadmap Epigenomics dataset](https://www.nature.com/articles/nature14248) using the expanded 18-state ChromHMM model on 98 tissues [as described here](https://egg2.wustl.edu/roadmap/web_portal/chr_state_learning.html).  

#### Gene features: protein

We collected the following gene features related to the protein encoded by each gene:  

| Feature name | Abbreviation | Description | Source |  
| :--- | :--- | :--- | :--- |  
| Amino acid length | `aa_length` | Length of canonical peptide sequence | Swissprot [Bairoch & Apweiler, _Nucleic Acids Research_, 2000)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102476) |  
| Molecular mass | `protein_mass` | Mass of protein in Daltons | Swissprot [Bairoch & Apweiler, _Nucleic Acids Research_, 2000)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102476) |  
| Number of protein interactors | `ppi_degree` | Number of proteins known to interact with this protein as per the Intact database | Swissprot [Bairoch & Apweiler, _Nucleic Acids Research_, 2000)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102476) |  
| Intramembrane domains | `intramembrane_domains` | Number of intramembrane domains | Swissprot [Bairoch & Apweiler, _Nucleic Acids Research_, 2000)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102476) |  
| Transmembrane domains | `transmembrane_domains` | Number of transmembrane domains | Swissprot [Bairoch & Apweiler, _Nucleic Acids Research_, 2000)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102476) |  
| Alpha helixes | `alpha_helixes` | Number of alpha helixes | Swissprot [Bairoch & Apweiler, _Nucleic Acids Research_, 2000)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102476) |  
| Beta sheets | `beta_sheets` | Number of beta sheets | Swissprot [Bairoch & Apweiler, _Nucleic Acids Research_, 2000)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102476) |  
| Metal binding sites | `metal_binding_sites` | Number of known binding sites for metal ions | Swissprot [Bairoch & Apweiler, _Nucleic Acids Research_, 2000)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102476) |  
| Nucleotide binding sites | `nucleotide_binding_sites` | Number of known nucleotide binding sites | Swissprot [Bairoch & Apweiler, _Nucleic Acids Research_, 2000)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102476) |  
| Active sites | `active_sites` | Number of known active sites | Swissprot [Bairoch & Apweiler, _Nucleic Acids Research_, 2000)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC102476) |  

#### Gene features: human genetic variation

We collected the following gene features related to reported genetic variation in humans:  

| Feature name | Abbreviation | Description | Source |  
| :--- | :--- | :--- | :--- |  
| Loss-of-function short variants in gnomAD | `gnomad_obs_syn` | Number of synonymous short variants observed in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _Nature_, 2020)](https://pubmed.ncbi.nlm.nih.gov/32461654) |  
| Loss-of-function short variants in gnomAD | `gnomad_obs_mis` | Number of missense short variants observed in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _Nature_, 2020)](https://pubmed.ncbi.nlm.nih.gov/32461654) |  
| Loss-of-function short variants in gnomAD | `gnomad_obs_lof` | Number of loss-of-function short variants observed in gnomAD | gnomAD v2.1 [(Karczewski _et al._, _Nature_, 2020)](https://pubmed.ncbi.nlm.nih.gov/32461654) |  
| All loss-of-function structural variants in gnomAD-SV | `gnomad_sv_lof_any` | Total number of loss-of-function structural variants observed in the non-neuro subset of gnomAD-SV | gnomAD-SV v2.1, non-neuro [(Collins _et al._, _Nature_, 2020)](https://www.nature.com/articles/s41586-020-2287-8) |  
| Loss-of-function deletions in gnomAD-SV | `gnomad_sv_lof_del` | Number of loss-of-function deletions observed in the non-neuro subset of gnomAD-SV | gnomAD-SV v2.1, non-neuro [(Collins _et al._, _Nature_, 2020)](https://www.nature.com/articles/s41586-020-2287-8) |  
| Non-deletion loss-of-function structural variants in gnomAD-SV | `gnomad_sv_lof_del` | Number of non-deletion loss-of-function structural variants observed in the non-neuro subset of gnomAD-SV | gnomAD-SV v2.1, non-neuro [(Collins _et al._, _Nature_, 2020)](https://www.nature.com/articles/s41586-020-2287-8) |  
| Whole-gene copy gains in gnomAD-SV | `gnomad_sv_cg` | Total number of whole-gene copy gains observed in the non-neuro subset of gnomAD-SV | gnomAD-SV v2.1, non-neuro [(Collins _et al._, _Nature_, 2020)](https://www.nature.com/articles/s41586-020-2287-8) |  
| Intragenic exonic duplications in gnomAD-SV | `gnomad_sv_ied` | Total number of intragenic exonic duplications observed in the non-neuro subset of gnomAD-SV | gnomAD-SV v2.1, non-neuro [(Collins _et al._, _Nature_, 2020)](https://www.nature.com/articles/s41586-020-2287-8) |  
| _De novo_ loss-of-function mutations in DDD | `ddd_dn_lof` | Number of _de novo_ loss-of-function mutations reported in the DDD | Deciphering Developmental Disorders [(Kaplanis _et al._, _bioRxiv_, 2020)](https://www.biorxiv.org/content/10.1101/797787v3) |  
| _De novo_ missense mutations in DDD | `ddd_dn_mis` | Number of _de novo_ loss-of-function mutations reported in the DDD | Deciphering Developmental Disorders [(Kaplanis _et al._, _bioRxiv_, 2020)](https://www.biorxiv.org/content/10.1101/797787v3) |  
| _De novo_ synonymous mutations in DDD | `ddd_dn_syn` | Number of _de novo_ loss-of-function mutations reported in the DDD | Deciphering Developmental Disorders [(Kaplanis _et al._, _bioRxiv_, 2020)](https://www.biorxiv.org/content/10.1101/797787v3) |  
| _De novo_ loss-of-function mutations in ASC probands | `asc_dn_lof` | Number of _de novo_ loss-of-function mutations reported in affected probands in the ASC | ASC [(Satterstrom _et al._, _Cell_, 2020)](https://pubmed.ncbi.nlm.nih.gov/31981491/) |  
| _De novo_ missense mutations in ASC probands | `asc_dn_mis` | Number of _de novo_ missense mutations reported in affected probands in the ASC | ASC [(Satterstrom _et al._, _Cell_, 2020)](https://pubmed.ncbi.nlm.nih.gov/31981491/) |  
| _De novo_ synonymous mutations in ASC probands | `asc_dn_syn` | Number of _de novo_ synonymous mutations reported in affected probands in the ASC | ASC [(Satterstrom _et al._, _Cell_, 2020)](https://pubmed.ncbi.nlm.nih.gov/31981491/) |  
| _De novo_ loss-of-function mutations in ASC siblings | `asc_unaffected_dn_lof` | Number of _de novo_ loss-of-function mutations reported in unaffected siblings in the ASC | ASC [(Satterstrom _et al._, _Cell_, 2020)](https://pubmed.ncbi.nlm.nih.gov/31981491/) |  
| _De novo_ missense mutations in ASC siblings | `asc_unaffected_dn_mis` | Number of _de novo_ missense mutations reported in unaffected siblings in the ASC | ASC [(Satterstrom _et al._, _Cell_, 2020)](https://pubmed.ncbi.nlm.nih.gov/31981491/) |  
| _De novo_ synonymous mutations in ASC siblings | `asc_unaffected_dn_syn` | Number of _de novo_ synonymous mutations reported in unaffected siblings in the ASC | ASC [(Satterstrom _et al._, _Cell_, 2020)](https://pubmed.ncbi.nlm.nih.gov/31981491/) |  
| Translocations in congenital anomalies | `redin_tloc` | Number of _de novo_ or disease-segregating gene-disruptive simple translocations in congenital anomalies | DGAP [(Redin _et al._, _Nat. Genet._, 2017)](https://pubmed.ncbi.nlm.nih.gov/27841880/) |  
| Inversions in congenital anomalies | `redin_inv` | Number of _de novo_ or disease-segregating gene-disruptive simple inversions in congenital anomalies | DGAP [(Redin _et al._, _Nat. Genet._, 2017)](https://pubmed.ncbi.nlm.nih.gov/27841880/) |  
| Complex rearrangements in congenital anomalies | `redin_cpx` | Number of _de novo_ or disease-segregating gene-disruptive complex rearrangements in congenital anomalies | DGAP [(Redin _et al._, _Nat. Genet._, 2017)](https://pubmed.ncbi.nlm.nih.gov/27841880/) |  
| Balanced chromosomal abnormalities in congenital anomalies | `redin_any_bca` | Number of _de novo_ or disease-segregating gene-disruptive balanced chromosomal abnormalities in congenital anomalies | DGAP [(Redin _et al._, _Nat. Genet._, 2017)](https://pubmed.ncbi.nlm.nih.gov/27841880/) |  

### Eigenfeature calculation  

Following annotation of all features described above, we collapsed correlated annotations to retain the principal components, or "eigenfeatures," that captured at least 99% of inter-gene variance.  

Prior to principal components analysis, we normalized the data in two steps, as follows:  
1. All variables were transformed using Box-Cox power transformations (unless otherwise specified in the table below), and
2. Following transformation, all variables were centered (mean = 0) and scaled (standard deviation = 1).  

Exceptions to Box-Cox power transformations are listed below:  

| Feature | Transformation | Reason |  
| :--- | :--- | :--- |  
| All `expression_component` features | None | By definition, features defined as principal components of other data need no additional transformation |  
| All `chromatin_component` features | None | By definition, features defined as principal components of other data need no additional transformation |  
| All de novo mutation data from ASC & DDD | None | Mutation counts are too sparse for transformation |  
| All Z-score features | None | By definition, Z-score features need no additional transformation |  
| All gnomAD-SV mutation count data | None | Mutation counts are too sparse for transformation |  
| All BCA counts from Redin _et al._ | None | Mutation counts are too sparse for transformation |  
| gnomAD constraint probabilities (pLI, pRec., pNull) | None | Beta-distributed probabilities are not well suited for power transformation |  
| `promoter_phastcons` | None | Data not well suited for power transformation |  


## Gene set definitions  

Throughout our analyses, we reference various subsets of genes. These are defined below.  

Note that all gene sets were required to have a unique match to a gene name from one of the canonical, protein-coding genes curated for this study (as defined above).  

All gene lists in the table below are available from `gs://rcnv_project/cleaned_data/genes/gene_lists/`, which is accessible to users with appropriate GCP permissions.  

| Gene set | Genes | Filename prefix | Source | Description |
| :--- | ---: | :--- | :--- | :---- | 
| All genes | 18,641 | `gencode.v19.canonical.pext_filtered` | Gencode v19 [Harrow _et al._, _Genome Res._, 2012](https://www.ncbi.nlm.nih.gov/pubmed/22955987) | Canonical transcripts from autosomal, protein-coding genes |
| LoF-constrained genes | 3,036 | `gnomad.v2.1.1.lof_constrained` | gnomAD v2.1.1 [Karczewski _et al._, _Nature_, 2020](https://www.nature.com/articles/s41586-020-2308-7) | pLI ≥ 0.9 or in the first LOEUF sextile |
| LoF-constrained (strict) genes | 2,710 | `gnomad.v2.1.1.lof_constrained_strict` | gnomAD v2.1.1 [Karczewski _et al._, _Nature_, 2020](https://www.nature.com/articles/s41586-020-2308-7) | pLI ≥ 0.9 and in the first LOEUF sextile |
| Missense-constrained genes | 3,019 | `gnomad.v2.1.1.mis_constrained` | gnomAD v2.1.1 [Karczewski _et al._, _Nature_, 2020](https://www.nature.com/articles/s41586-020-2308-7) | Missense Z ≥ 3 or in the first MOEUF sextile |
| Likely unconstrained genes | 5,578 | `gnomad.v2.1.1.likely_unconstrained` | gnomAD v2.1.1 [Karczewski _et al._, _Nature_, 2020](https://www.nature.com/articles/s41586-020-2308-7) | LOEUF ≥ 1, mis. OEUF ≥ 1, synonymous Z-score ~ [-3, 3], pLI ≤ 0.1, LoF O/E in upper 50% of all genes, mis. OE in upper 50% of all genes, observed LoF > 0, observed mis. > 0 |
| Mutation-tolerant genes | 2,207 | `gnomad.v2.1.1.mutation_tolerant` | gnomAD v2.1.1 [Karczewski _et al._, _Nature_, 2020](https://www.nature.com/articles/s41586-020-2308-7) | pLI ≥ 0.01, the last third of LOEUF, missense Z-score ≤ 0, missense OEUF ≥ 1, and synonymous Z-score ~ (-3, 3) |
| Cell essential genes | 664 | `cell_essential` | [Hart _et al._, _G3_, 2017](https://www.g3journal.org/content/7/8/2719) | Genes essential in human cell lines as determined by CRISPR/Cas9 screens. Curated in [Karczewski _et al._, _Nature_, 2020](https://www.nature.com/articles/s41586-020-2308-7) |
| Cell nonessential genes | 711 | `cell_nonessential` | [Hart _et al._, _G3_, 2017](https://www.g3journal.org/content/7/8/2719) | Genes nonessential in human cell lines as determined by CRISPR/Cas9 screens. Curated in [Karczewski _et al._, _Nature_, 2020](https://www.nature.com/articles/s41586-020-2308-7) |
| Mouse heterozygous LoF lethal | 370 | `mouse_het_lethal` | [Motenko _et al._, _Mamm. Genome_, 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4534495/) | Genes lethal in mouse models when heterozygously inactivated in mice. Curated in [Karczewski _et al._, _Nature_, 2020](https://www.nature.com/articles/s41586-020-2308-7) |
| Mouse dispensable | 216 | `mouse_dispensable` | [Motenko _et al._, _Mamm. Genome_, 2015](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4534495/) | Genes with no phenotype reported when heterozygously inactivated in mice. Curated in [Karczewski _et al._, _Nature_, 2020](https://www.nature.com/articles/s41586-020-2308-7) |
| Olfactory receptors | 229 | `olfactory_receptors` | [Braschi _et al._, _Nucleic Acids Res_, 2019](https://pubmed.ncbi.nlm.nih.gov/30304474/) | Genes from any HUGO-recognized family of olfactory receptor genes. |
| Low expression-variable genes | 2,682 | `gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.low_expression_variable` | GTex v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) | ≥ 5 mean TPM and 25<sup>th</sup> percentile ≥ 0 in ≥ 3 tissues; mean ≥ 0.5% of samples are low expression outliers; mean TPM CV ≥ 0.5 |
| Low expression-invariant genes | 617 | `gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.low_expression_invariant` | GTex v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) | ≥ 5 mean TPM and 25<sup>th</sup> percentile ≥ 0 in ≥ 3 tissues; mean ≤ 0.1% of samples are low expression outliers; mean TPM CV ≤ 0.3; no samples with TPM < 1 in any tissue |
| High expression-variable genes | 1,924 | `gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.high_expression_variable` | GTex v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) | ≥ 5 mean TPM and 25<sup>th</sup> percentile ≥ 0 in ≥ 3 tissues; mean ≥ 5% of samples are high expression outliers; mean TPM CV ≥ 0.5 |
| High expression-invariant genes | 279 | `gencode.v19.canonical.pext_filtered.GTEx_v7_variable_expressors.high_expression_invariant` | GTex v7 [(GTEx Consortium, _Nature_, 2017)](https://www.ncbi.nlm.nih.gov/pubmed/29022597) | ≥ 5 mean TPM and 25<sup>th</sup> percentile ≥ 0 in ≥ 3 tissues; mean ≤ 1% of samples are low expression outliers; mean TPM CV ≤ 0.3 |
| Phenotypic abnormality (HP:0000118)-associated genes | 4,054 | `HP0000118.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000118 |
| Abnormality of the nervous system (HP:0000707)-associated genes | 3,084 | `HP0000707.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000707 |
| Behavioral abnormality (HP:0000708)-associated genes | 1,481 | `HP0000708.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000708 |
| NA (UNKNOWN)-associated genes | 4,054 | `UNKNOWN.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to UNKNOWN |
| Abnormal nervous system morphology (HP:0012639)-associated genes | 2,347 | `HP0012639.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0012639 |
| Impairment in personality functioning (HP:0031466)-associated genes | 552 | `HP0031466.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0031466 |
| Neurodevelopmental abnormality (HP:0012759)-associated genes | 2,179 | `HP0012759.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0012759 |
| Abnormality of the immune system (HP:0002715)-associated genes | 1,548 | `HP0002715.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0002715 |
| Abnormality of metabolism/homeostasis (HP:0001939)-associated genes | 2,020 | `HP0001939.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001939 |
| Autoimmunity (HP:0002960)-associated genes | 158 | `HP0002960.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0002960 |
| Abnormality of the musculoskeletal system (HP:0033127)-associated genes | 3,025 | `HP0033127.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0033127 |
| Abnormality of the cardiovascular system (HP:0001626)-associated genes | 2,132 | `HP0001626.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001626 |
| Morphological central nervous system abnormality (HP:0002011)-associated genes | 2,212 | `HP0002011.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0002011 |
| Bipolar affective disorder (HP:0007302)-associated genes | 23 | `HP0007302.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0007302 |
| Abnormal peripheral nervous system morphology (HP:0000759)-associated genes | 443 | `HP0000759.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000759 |
| Abnormality of the vasculature (HP:0002597)-associated genes | 1,480 | `HP0002597.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0002597 |
| Abnormal joint morphology (HP:0001367)-associated genes | 892 | `HP0001367.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001367 |
| Autistic behavior (HP:0000729)-associated genes | 465 | `HP0000729.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000729 |
| Atrophy/Degeneration affecting the central nervous system (HP:0007367)-associated genes | 694 | `HP0007367.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0007367 |
| Abnormal fear/anxiety-related behavior (HP:0100852)-associated genes | 309 | `HP0100852.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0100852 |
| Abnormal circulating cholesterol concentration (HP:0003107)-associated genes | 91 | `HP0003107.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0003107 |
| Schizophrenia (HP:0100753)-associated genes | 54 | `HP0100753.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0100753 |
| Abnormality of prenatal development or birth (HP:0001197)-associated genes | 779 | `HP0001197.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001197 |
| Peripheral neuropathy (HP:0009830)-associated genes | 374 | `HP0009830.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0009830 |
| Seizure (HP:0001250)-associated genes | 1,501 | `HP0001250.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001250 |
| Atherosclerosis (HP:0002621)-associated genes | 67 | `HP0002621.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0002621 |
| Abnormality of head or neck (HP:0000152)-associated genes | 2,591 | `HP0000152.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000152 |
| Abnormality of higher mental function (HP:0011446)-associated genes | 2,143 | `HP0011446.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0011446 |
| Hyperactivity (HP:0000752)-associated genes | 385 | `HP0000752.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000752 |
| Hyperlipidemia (HP:0003077)-associated genes | 96 | `HP0003077.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0003077 |
| Abnormality of movement (HP:0100022)-associated genes | 1,635 | `HP0100022.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0100022 |
| Abnormality of the face (HP:0000271)-associated genes | 2,338 | `HP0000271.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000271 |
| Abnormal heart morphology (HP:0001627)-associated genes | 1,128 | `HP0001627.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001627 |
| Abnormal inflammatory response (HP:0012647)-associated genes | 864 | `HP0012647.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0012647 |
| Abnormality of the musculature (HP:0003011)-associated genes | 2,280 | `HP0003011.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0003011 |
| Abnormality of the integument (HP:0001574)-associated genes | 2,084 | `HP0001574.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001574 |
| Intellectual disability (HP:0001249)-associated genes | 1,515 | `HP0001249.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001249 |
| Abnormality of the digestive system (HP:0025031)-associated genes | 2,292 | `HP0025031.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0025031 |
| Rheumatoid arthritis (HP:0001370)-associated genes | 19 | `HP0001370.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001370 |
| Abnormality of the respiratory system (HP:0002086)-associated genes | 1,587 | `HP0002086.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0002086 |
| Abnormality of brain morphology (HP:0012443)-associated genes | 2,029 | `HP0012443.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0012443 |
| Growth abnormality (HP:0001507)-associated genes | 2,075 | `HP0001507.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001507 |
| Abnormality of the genitourinary system (HP:0000119)-associated genes | 2,034 | `HP0000119.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0000119 |
| Gout (HP:0001997)-associated genes | 16 | `HP0001997.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0001997 |
| Abnormality of digestive system morphology (HP:0025033)-associated genes | 820 | `HP0025033.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0025033 |
| Abnormal axial skeleton morphology (HP:0009121)-associated genes | 2,064 | `HP0009121.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0009121 |
| Sleep disturbance (HP:0002360)-associated genes | 452 | `HP0002360.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0002360 |
| Abnormality of the peripheral nervous system (HP:0410008)-associated genes | 328 | `HP0410008.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0410008 |
| Abnormal muscle physiology (HP:0011804)-associated genes | 1,972 | `HP0011804.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0011804 |
| Abnormal central motor function (HP:0011442)-associated genes | 1,476 | `HP0011442.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0011442 |
| Involuntary movements (HP:0004305)-associated genes | 871 | `HP0004305.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0004305 |
| Abnormality of body weight (HP:0004323)-associated genes | 1,389 | `HP0004323.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0004323 |
| Abnormality of skin morphology (HP:0011121)-associated genes | 1,588 | `HP0011121.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0011121 |
| Abnormal myelination (HP:0012447)-associated genes | 449 | `HP0012447.HPOdb` | HPO database (accessed Mar 2022) [Köhler _et al._, _Nucleic Acids Res._, 2018](https://academic.oup.com/nar/article/47/D1/D1018/5198478) | Genes linked to HP:0012447 |
| ClinGen dominant haploinsufficient genes (high confidence) | 208 | `ClinGen.hc_haploinsufficient` | ClinGen gene curation map (accessed 10 UTC) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant haploinsufficient genes scored at high confidence by ClinGen |
| ClinGen dominant haploinsufficient genes (medium confidence) | 43 | `ClinGen.mc_haploinsufficient` | ClinGen gene curation map (accessed 10 UTC) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant haploinsufficient genes scored at medium confidence by ClinGen |
| ClinGen dominant haploinsufficient genes (high or medium confidence) | 251 | `ClinGen.hmc_haploinsufficient` | ClinGen gene curation map (accessed 10 UTC) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant haploinsufficient genes scored at high or medium confidence by ClinGen |
| ClinGen dominant haploinsufficient genes (low confidence) | 78 | `ClinGen.lc_haploinsufficient` | ClinGen gene curation map (accessed 10 UTC) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant haploinsufficient genes scored at low confidence by ClinGen |
| ClinGen dominant haploinsufficient genes (any confidence) | 329 | `ClinGen.all_haploinsufficient` | ClinGen gene curation map (accessed 10 UTC) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant haploinsufficient genes scored at any confidence by ClinGen |
| ClinGen dominant triplosensitive genes (high confidence) | 1 | `ClinGen.hc_triplosensitive` | ClinGen gene curation map (accessed 10 UTC) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant triplosensitive genes scored at high confidence by ClinGen |
| ClinGen dominant triplosensitive genes (medium confidence) | 1 | `ClinGen.mc_triplosensitive` | ClinGen gene curation map (accessed 10 UTC) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant triplosensitive genes scored at medium confidence by ClinGen |
| ClinGen dominant triplosensitive genes (high or medium confidence) | 2 | `ClinGen.hmc_triplosensitive` | ClinGen gene curation map (accessed 10 UTC) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant triplosensitive genes scored at high or medium confidence by ClinGen |
| ClinGen dominant triplosensitive genes (low confidence) | 13 | `ClinGen.lc_triplosensitive` | ClinGen gene curation map (accessed 10 UTC) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant triplosensitive genes scored at low confidence by ClinGen |
| ClinGen dominant triplosensitive genes (any confidence) | 15 | `ClinGen.all_triplosensitive` | ClinGen gene curation map (accessed 10 UTC) [Strande _et al._, _Am. J. Hum. Genet._, 2017](https://www.ncbi.nlm.nih.gov/pubmed/28552198) | Dominant triplosensitive genes scored at any confidence by ClinGen |
| DECIPHER dominant dev. disorder genes (loss-of-function mechanism; high confidence) | 300 | `DDG2P.hc_lof` | DECIPHER/DDG2P (accessed 10 UTC) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with a loss-of-function mechanism scored as "confirmed" in DECIPHER |
| DECIPHER dominant dev. disorder genes (loss-of-function mechanism; medium confidence) | 135 | `DDG2P.mc_lof` | DECIPHER/DDG2P (accessed 10 UTC) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with a loss-of-function mechanism scored as "probable" in DECIPHER |
| DECIPHER dominant dev. disorder genes (loss-of-function mechanism; high or medium confidence) | 435 | `DDG2P.hmc_lof` | DECIPHER/DDG2P (accessed 10 UTC) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with a loss-of-function mechanism scored as "confirmed" or "probable" in DECIPHER |
| DECIPHER dominant dev. disorder genes (loss-of-function mechanism; low confidence) | 63 | `DDG2P.lc_lof` | DECIPHER/DDG2P (accessed 10 UTC) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with a loss-of-function mechanism scored as "possible" in DECIPHER |
| DECIPHER dominant dev. disorder genes (loss-of-function mechanism; any confidence) | 498 | `DDG2P.all_lof` | DECIPHER/DDG2P (accessed 10 UTC) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with a loss-of-function mechanism scored as "confirmed", "probable", or "possible" in DECIPHER |
| DECIPHER dominant dev. disorder genes (activating/gain-of-function mechanism; high confidence) | 82 | `DDG2P.hc_gof` | DECIPHER/DDG2P (accessed 10 UTC) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an activating/gain-of-function mechanism scored as "confirmed" in DECIPHER |
| DECIPHER dominant dev. disorder genes (activating/gain-of-function mechanism; medium confidence) | 25 | `DDG2P.mc_gof` | DECIPHER/DDG2P (accessed 10 UTC) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an activating/gain-of-function mechanism scored as "probable" in DECIPHER |
| DECIPHER dominant dev. disorder genes (activating/gain-of-function mechanism; high or medium confidence) | 107 | `DDG2P.hmc_gof` | DECIPHER/DDG2P (accessed 10 UTC) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an activating/gain-of-function mechanism scored as "confirmed" or "probable" in DECIPHER |
| DECIPHER dominant dev. disorder genes (activating/gain-of-function mechanism; low confidence) | 9 | `DDG2P.lc_gof` | DECIPHER/DDG2P (accessed 10 UTC) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an activating/gain-of-function mechanism scored as "possible" in DECIPHER |
| DECIPHER dominant dev. disorder genes (activating/gain-of-function mechanism; any confidence) | 116 | `DDG2P.all_gof` | DECIPHER/DDG2P (accessed 10 UTC) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an activating/gain-of-function mechanism scored as "confirmed", "probable", or "possible" in DECIPHER |
| DECIPHER dominant dev. disorder genes (other/unknown coding mechanism; high confidence) | 145 | `DDG2P.hc_other` | DECIPHER/DDG2P (accessed 10 UTC) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an other/unknown coding mechanism scored as "confirmed" in DECIPHER |
| DECIPHER dominant dev. disorder genes (other/unknown coding mechanism; medium confidence) | 96 | `DDG2P.mc_other` | DECIPHER/DDG2P (accessed 10 UTC) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an other/unknown coding mechanism scored as "probable" in DECIPHER |
| DECIPHER dominant dev. disorder genes (other/unknown coding mechanism; high or medium confidence) | 241 | `DDG2P.hmc_other` | DECIPHER/DDG2P (accessed 10 UTC) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an other/unknown coding mechanism scored as "confirmed" or "probable" in DECIPHER |
| DECIPHER dominant dev. disorder genes (other/unknown coding mechanism; low confidence) | 36 | `DDG2P.lc_other` | DECIPHER/DDG2P (accessed 10 UTC) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an other/unknown coding mechanism scored as "possible" in DECIPHER |
| DECIPHER dominant dev. disorder genes (other/unknown coding mechanism; any confidence) | 277 | `DDG2P.all_other` | DECIPHER/DDG2P (accessed 10 UTC) [Wright _et al._, _Lancet_, 2015](https://www.ncbi.nlm.nih.gov/pubmed/25529582) | Dominant dev. disorder genes with an other/unknown coding mechanism scored as "confirmed", "probable", or "possible" in DECIPHER |
| High-confidence oncogenes | 122 | `COSMIC.hc_oncogenes` | COSMIC v91 [Sondka _et al._, _Nat. Rev. Cancer_, 2018](https://www.nature.com/articles/s41568-018-0060-1) | Tier 1 dominant oncogenes due to amplification, missense, or 'other' mutations |  
