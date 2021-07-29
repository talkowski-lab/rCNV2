# Other Data Curation  

Below, we describe curation steps for other datasets not already listed in a different subdirectory, including:  

* Various microarray probesets  
* Structural variant calls from the Human Genome Diversity Panel
* Previously reported genomic disorders  
* Predicted NAHR-mediated CNV loci  
* _De novo_ point mutations from exome sequencing studies of neurodevelopmental disorders
* Balanced chromosomal abnormality breakpoints from sequencing studies of congenital anomalies


---   

## Microarray probesets used for control samples  

We collected the microarray SNP probesets used for all control samples from each cohort included in our analyses.  

The code to curate these data are provided in `curate_control_probesets.sh`.  

These datasets were accessed as follows:  

| Array | Abbreviation | Source |  
| :--- | :--- | :--- |  
| Affymetrix Human Mapping 500K | `affy_500k` | UCSC Genome Browser (table `snpArrayAffy250Nsp` combined with table `snpArrayAffy250Sty`) |  
| Affymetrix Genome-Wide Human SNP Array 5.0 | `affy_5.0` | UCSC Genome Browser (table `snpArrayAffy5`) |  
| Affymetrix Genome-Wide Human SNP Array 5.0 | `affy_6.0` | UCSC Genome Browser (table `snpArrayAffy6`) |  
| Affymetrix CytoScan HD Array | `affy_cyto_hd` | Downloaded from ThermoFisher website ([link](https://www.thermofisher.com/order/catalog/product/901835#/901835)) |  
| Agilent Custom 180k SNP Array | `agilent_gdx_180k` | Provided by GeneDx, Inc. |  
| Illumina Human Hap 300 | `illumina_300k` | UCSC Genome Browser (table `snpArrayIllumina300`) |  
| Illumina Human Hap 550 | `illumina_550k` | UCSC Genome Browser (table `snpArrayIllumina550`) |  
| Illumina Human610-Quad | `illumina_610k_quad` | Downloaded from Illumina website ([link](https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/human610/human610-quadv1_h.zip)) |  
| Illumina Human Hap 650 | `illumina_650k` | UCSC Genome Browser (table `snpArrayIllumina650`) |  
| Illumina Human1M-Duo | `illumina_1m_duo` | UCSC Genome Browser (table `snpArrayIllumina1M`) |  
| Illumina Infinium Global Screening Array | `illumina_gsa` | Provided by the Estonian Biobank |  
| Illumina Infinium OmniExpress | `omniexpress` | Downloaded from Illumina website ([link](https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/humanomniexpress-24/v1-3/infinium-omniexpress-24-v1-3-a1-bed.zip)) |  
| Illumina Infinium Omni2.5 | `omni_2.5` | Downloaded from Illumina website ([link](https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/humanomni25/v1-5/infinium-omni2-5-8v1-5-a1-bed.zip)) |  
| Illumina Multi-Ethnic Global Array | `illumina_mega` | Downloaded from Illumina website ([link](https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/multiethnic-global/multi-ethnic-global-8-d1-bed.zip)) |  
| ThermoFisher Applied Biosystems UK Biobank Axiom | `ukbbaxiom` | Downloaded from ThermoFisher website ([link](https://www.thermofisher.com/order/catalog/product/902502?us&en#/902502?us&en)) |  

---  

## Structural variant callset from whole-genome sequencing of the Human Genome Diversity Panel  

We curated structural variants from the Human Genome Diversity Panel as described in [Almarri _et al._, _Cell_ (2020)](https://pubmed.ncbi.nlm.nih.gov/32531199/).  

The code to curate the HGDP SV callset is provided in `curate_hgdp_svs.sh`.  

---  

## Genomic disorders  

We curated lists of previously reported genomic disorders (GDs), defined as genomic intervals where rare CNVs have been associated with one or more diseases.  

The code for this process is provided in `curate_known_gds.sh`.  

For this purpose, we integrated lists of GDs from the following publications and public resources:  

| Resource | Deletion GDs | Duplication GDs | Citation |  
| :--- | ---: | ---: | --- |  
| DECIPHER CNV Syndromes | 40 | 14 | [Firth _et al._, _Am. H. Hum. Genet._ (2009)](http://dx.doi.org/10.1016/j.ajhg.2009.03.010) |  
| ClinGen Pathogenic CNV Regions\* | 43\* | 26\*  | [Riggs _et al._, _Clin. Genet._ (2012)](https://www.ncbi.nlm.nih.gov/pubmed/22097934) |  
| UK BioBank | 24 | 30 | [Owen _et al._, _BMC Genomics_ (2018)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6278042/) |  
| Girirajan _et al._ | 39 | 33 | [Girirajan _et al._, _New Engl. J. Med._ (2012)](https://www.nejm.org/doi/full/10.1056/NEJMoa1200395) |  
| Dittwald _et al._ | 48 | 18 | [Dittwald _et al._, _Genome Res._ (2013)](https://www.ncbi.nlm.nih.gov/pubmed/23657883) |  
| Stefansson _et al._ | 17 | 12 | [Stefansson _et al._, _Nature_ (2014)](https://pubmed.ncbi.nlm.nih.gov/24352232/) |  

_\* Note: ClinGen counts represent regions scored at high and medium confidence only._  

From these existing datasets, we curated three distinct GD lists for the analyses in this study:  

1. **High-confidence GDs**: regions covered by at least four sources.  
2. **Medium-confidence GDs**: regions covered by two or three sources.  
3. **Low-confidence GDs**: regions covered by only one source  

After overlapping reported GDs from the sources above, we trimmed overlapping segmental duplications overlapping the boundaries of the GD interval.  

---  

## Predicted NAHR-mediated CNVs  

We defined a set of loci where NAHR-mediated CNVs might be predicted to occur based on the genomic properties of the flanking regions.  

The code for this process is provided in `predict_nahr_cnvs.sh`.  

To build this set of loci, we first defined pairs of [segmental duplications](https://genome.ucsc.edu/cgi-bin/hgTables) meeting the following criteria:  
1. Same chromosome
2. Distance ≥ 100kb & ≤ 10Mb
3. Both segmental duplications ≥ 1kb
4. Strict homology (no indels) ≥ 90%
5. Direct orentation of repeats (_i.e._, same strand)
6. Total intervening sequence has ≤30% coverage by the blacklist used [during CNV curation](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#curation-steps-rare-cnvs)  
7. At least 100kb of intervening sequence remaining after subtracting blacklist regions

After defining candidate pairs of segmental duplications (above), we collapsed overlapping pairs into predicted NAHR-mediated CNV loci while requiring:
1. Both ends of their respective intervals to be within 1Mb of each other  
2. \>50% reciprocal overlap of intervening sequence  

For each cluster of segmental duplication pairs, we retained the pair with the smallest intervening (_i.e._, spanning) distance, and used the innermost coordinates for analysis purposes.  


---  

## _De novo_ mutations in neurodevelopmental disorders  

TBD

---  

## Chromosomal rearrangement breakpoints in congenital anomalies  

We curated breakpoints from balanced chromosomal abnormalities (BCAs) defined at sequence resolution in subjects with congenital anomalies from [Redin _et al._, _Nat. Genet._ (2016)](https://www.nature.com/articles/ng.3720).  

We extracted all breakpoints as reported in the supplement of Redin _et al._, and further restricted the dataset as follows:  
1. Breakpoint must be resolved at sequence resolution; 
2. Rearrangement must either be confirmed as _de novo_ in affected subject, or confirmed to segregate with disease in the family (if multiple members were affected).  
