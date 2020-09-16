#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot summary distributions of outcome of gene-based fine-mapping for rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load table of PIPs from all genes considered in fine-mapping
load.allgenes <- function(allgenes.in){
  ag <- read.table(allgenes.in, header=T, sep="\t", comment.char="")
  colnames(ag)[1] <- gsub("X.", "", colnames(ag)[1], fixed=T)
  return(ag)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot histogram of PIPs for all genes in all credible sets
plot.pip.hist <- function(allgenes, conf.cutoff=0.15, vconf.cutoff=0.85,
                          breaks=seq(0, 1, 0.025), parmar=c(2.15, 2.85, 1, 0.6)){
  # Get plotting data
  pips <- lapply(c("DEL", "DUP"), function(cnv){
    allgenes$PIP[which(allgenes$cnv==cnv & !is.na(allgenes$credible_set))]
  })
  counts.del <- hist(pips[[1]], breaks=breaks, plot=F)$counts
  counts.dup <- hist(pips[[2]], breaks=breaks, plot=F)$counts
  counts <- counts.del + counts.dup
  bar.colors <- lapply(c("DEL", "DUP"), function(cnv){
    sapply(breaks[-1], function(x){
      x <- round(x, 2)
      if(x <= conf.cutoff){
        cnv.whites[cnv]
      }else if(x <= vconf.cutoff){
        control.cnv.colors[cnv]
      }else{
        cnv.colors[cnv]
      }
    })
  })
  
  # Prepare plot area
  par(bty="n", mar=parmar)
  plot(NA, xlim=c(0, 1), ylim=c(0, 1.025*max(counts)),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  # rect(xleft=c(conf.cutoff, vconf.cutoff), xright=c(vconf.cutoff, 1),
  #      ybottom=0, ytop=par("usr")[4], col=c(cnv.whites[cnv], control.cnv.colors[cnv]),
  #      border=NA, bty="n")
  abline(v=c(conf.cutoff, vconf.cutoff), lty=2, col=blueblack)
  text(x=c(conf.cutoff, vconf.cutoff)+0.05, y=sum(par("usr")[3:4])/2, srt=90,
       labels=c("Confident", "Highly Confident"), cex=5/6)
  
  # Add stacked histogram
  rect(xleft=breaks[-length(breaks)], xright=breaks[-1],
       ybottom=c(rep(0, length(counts.del)), counts.del), 
       ytop=c(counts.del, counts.del+counts.dup),
       col=unlist(bar.colors), border=NA, bty="n", xpd=T)
  rect(xleft=breaks[-length(breaks)], xright=breaks[-1],
       ybottom=rep(0, length(counts.del)), 
       ytop=counts.del+counts.dup,
       col=NA, border=blueblack, xpd=T)
  
  # Add axes
  x.ax.at <- axTicks(1)
  y.ax.at <- axTicks(2)
  axis(1, x.ax.at, tck=-0.025, labels=NA, col=blueblack)
  sapply(x.ax.at, function(x){axis(1, x, tick=F, line=-0.7)})
  mtext(1, line=1.25, text="Posterior Inclusion Prob. (PIP)")
  axis(2, c(0, 10e10), tck=0, labels=NA, col=blueblack)
  axis(2, y.ax.at, tck=-0.025, labels=NA, col=blueblack)
  axis(2, y.ax.at, tick=F, las=2, line=-0.6)
  mtext(2, line=1.9, text="Genes in Credible Sets")
}


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
args <- parse_args(OptionParser(usage=paste("%prog credsets.bed assocs.bed allgenes.tsv out.prefix"),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop(paste("Four positional arguments required: credsets.bed, assocs.bed,",
             "allgenes.tsv, and out.prefix\n"))
}

# Writes args & opts to vars
credsets.in <- args$args[1]
assocs.in <- args$args[2]
allgenes.in <- args$args[3]
out.prefix <- args$args[4]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# credsets.in <- "~/scratch/rCNV.final_genes.credible_sets.bed.gz"
# assocs.in <- "~/scratch/rCNV.final_genes.associations.bed.gz"
# allgenes.in <- "~/scratch/all_genes_from_fine_mapping.tsv"
# out.prefix <- "~/scratch/finemap_distribs_test"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/gene_association/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load credible sets and associations
credsets <- load.credsets(credsets.in)
assocs <- load.associations(assocs.in)

# Load PIPs for all genes considered in fine-mapping
allgenes <- load.allgenes(allgenes.in)

# Plot histograms of PIPs for all genes in credible sets
pdf(paste(out.prefix, "finemapped_distribs.credset_pip.pdf", sep="."),
    height=2.2, width=2.7)
plot.pip.hist(allgenes)
dev.off()
