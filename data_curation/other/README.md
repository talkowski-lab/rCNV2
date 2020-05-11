# Other Data Curation  

Below, we describe curation steps for other datasets not already listed in a different subdirectory, including:  

* Previously reported genomic disorders  
* Predicted NAHR-mediated CNV loci  


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

_\* Note: ClinGen counts represent regions scored at high and medium confidence only._  

From these existing datasets, we curated two distinct GD lists for the analyses in this study:  

1. **Established GDs**: regions covered by at least four sources.  
2. **Candidate GDs**: regions covered by two or three sources.  

After overlapping reported GDs from the sources above, we further:
1. excluded all consensus GDs with ≥50% coverage by common CNVs from three WGS-based resources as [described elsewhere](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV/#curation-steps-rare-cnvs) or by common CNVs in [control samples from any metacohort in this study](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV/#case-control-metacohorts) (prior to [CNV frequency filtering](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV/#curation-steps-rare-cnvs)), and
2. trimmed overlapping segmental duplications overlapping the boundaries of the GD interval.  

---  

## Predicted NAHR-mediated CNVs  

We defined a set of loci where NAHR-mediated CNVs might be predicted to occur based on the genomic properties of the flanking regions.  

The code for this process is provided in `predict_nahr_cnvs.sh`.  

To build this set of loci, we first defined pairs of [segmental duplications](https://genome.ucsc.edu/cgi-bin/hgTables) meeting the following criteria:  
1. Same chromosome
2. Distance ≥ 100kb & ≤ 10Mb
3. Both segmental duplications ≥ 1kb
4. Strict homology (no indels) ≥ 95%
5. Direct orentation of repeats (_i.e._, same strand)
6. Total intervening sequence has ≤30% coverage by the blacklist used [during CNV curation](https://github.com/talkowski-lab/rCNV2/tree/master/data_curation/CNV#curation-steps-rare-cnvs)  
7. At least 100kb of intervening sequence remaining after subtracting blacklist regions

After defining candidate pairs of segmental duplications (above), we collapsed overlapping pairs into predicted NAHR-mediated CNV loci while requiring:
1. Both ends of their respective intervals to be within 1Mb of each other  
2. \>50% reciprocal overlap of intervening sequence  

For each cluster of segmental duplication pairs, we retained the pair with the smallest intervening (_i.e._, spanning) distance, and used the innermost coordinates for analysis purposes.  
