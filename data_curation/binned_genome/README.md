# Genome-Wide Bin Curation  

## Genome-Wide Bin Curation  

We segmented the genome into a set of regular, sequential bins for association testing. This process is described below.  

All commands executed to filter the CNV data as descrbed below are contained in `create_genome_bins.sh`.  

### Dependencies  

#### Athena  
Several commands used in the creation and processing of binned genome files require Athena, a toolkit currently in development by the Talkowski Lab.  

Athena is not currently publicly available. For now, the current Athena build is installed into the rCNV Docker, and can be accessed as follows:  
```
$ docker pull talkowski/rcnv

$ docker run --rm -it talkowski/rcnv

(rcnv) root@hash:/# athena --help

Usage: athena [OPTIONS] COMMAND [ARGS]...

Options:
  --help  Show this message and exit.

Commands:
  annotate    Annotate bins
  count-sv    Intersect SV and bins
  decomp      Decompose bin annotations
  make-bins   Create sequential bins
  query       Mutation rate lookup
  vcf-filter  Filter an input VCF
  vcf-stats   Get SV size & spacing 
```
