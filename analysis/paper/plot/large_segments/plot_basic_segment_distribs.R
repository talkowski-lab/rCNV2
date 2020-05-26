#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot basic distributions for large segments for rCNV2 manuscript


options(stringsAsFactors=F, scipen=1000)


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(funr, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog loci.bed segs.tsv out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop(paste("Three positional arguments required: loci.bed, segs.tsv, output_prefix\n", sep=" "))
}

# Writes args & opts to vars
loci.in <- args$args[1]
segs.in <- args$args[2]
out.prefix <- args$args[3]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d1.master_segments.bed.gz"
# out.prefix <- "~/scratch/final_segs"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/large_segments/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load loci & segment table
loci <- load.loci(loci.in)
segs <- load.segment.table(segs.in)

# Merge loci & segment data for genome-wide significant sites only
gw <- merge.loci.segs(loci, segs)

# Swarmplot of segment size
pdf(paste(out.prefix, "gw.segment_size_distribs.pdf", sep="."),
    height=2, width=2.3)
gw.simple.vioswarm(gw, y=log10(gw$size), add.y.axis=F, parmar=c(1.2, 4, 0.3, 0.3))
axis(2, at=log10(logscale.minor), tck=-0.03, col=blueblack, labels=NA)
axis(2, at=log10(c(300000, 600000, 1000000, 3000000)), tick=F, 
     las=2, line=-0.65, labels=c("300kb", "600kb", "1Mb", "3Mb"))
mtext(2, line=2.75, text=bquote("log"[10] * "(Segment Size)"))
dev.off()

# Swarmplot of genes per segment
pdf(paste(out.prefix, "gw.n_genes_per_segment.pdf", sep="."),
    height=2, width=2.3)
gw.simple.vioswarm(gw, y=gw$n_genes, 
                   add.y.axis=T, ytitle="Genes in Segment", 
                   parmar=c(1.2, 3, 0.3, 0.3))
dev.off()

# Swarmplot of gene density
pdf(paste(out.prefix, "gw.basic_gene_density.pdf", sep="."),
    height=2, width=2.3)
gw.simple.vioswarm(gw, y=100000*(gw$n_genes/gw$size), 
                   add.y.axis=T, ytitle="Genes per 100kb", 
                   parmar=c(1.2, 3, 0.3, 0.3))
dev.off()

# Scatterplot of size vs genes
pdf(paste(out.prefix, "gw.genes_vs_size.pdf", sep="."),
    height=2, width=2)
gw.scatter(gw, x=log10(gw$size), y=gw$n_genes, 
           x.at=log10(logscale.minor), x.labs=rep(NA, length(logscale.minor)),
           ytitle="Genes in Segment", parmar=c(2.75, 2.75, 0.3, 0.3))
axis(1, at=log10(c(300000,1000000, 3000000)), tick=F, 
     line=-0.7, labels=c("300kb", "1Mb", "3Mb"))
mtext(1, line=1.35, text=bquote("log"[10] * "(Segment Size)"))
dev.off()
 
