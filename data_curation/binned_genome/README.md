# Genome-Wide Bin Curation  

## Genome-Wide Bin Curation  

We segmented the genome into a set of regularized sliding windows for association testing. This process is described in the methods of Collins _et al._, 2022.    

All commands executed to curate genome-wide sliding windows are contained in `create_genome_bins.sh`.  

In practice, the commands in `create_genome_bins.sh` were parallelized in [FireCloud/Terra](https://portal.firecloud.org) with `create_genome_bins.wdl`.  

