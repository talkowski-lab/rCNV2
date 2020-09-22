#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot gene set enrichments vs. dosage sensitivity scores for rCNV2 manuscript


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load multiple gene sets from .tsv input
load.genesets <- function(genesets.in, elig.genes=NULL){
  glist <- read.table(genesets.in, header=F, sep="\t")
  gsets <- lapply(glist[, 2], function(path){
    genes <- sort(unique(as.character(read.table(path, header=F)[, 1])))
    if(!is.null(elig.genes)){
      genes <- genes[which(genes %in% elig.genes)]
    }
    return(genes)
  })
  names(gsets) <- glist[, 1]
  return(gsets)
}

# Split genes into groups based on score bin
bin.genes <- function(scores, score, n.bins=10){
  breaks <- seq(0, 1, length.out=n.bins+1)
  lapply(2:length(breaks), function(i){
    scores$gene[which(scores[, score] >= breaks[i-1]
                      & scores[, score] <= breaks[i])]
  })
}

# Compute the enrichment of a single set of genes vs. score bins
calc.enrichment.single <- function(genes, score.groups){
  n.genes <- length(genes)
  n.univ <- length(unlist(score.groups))
  baseline <- n.genes / n.univ
  enrich <- as.data.frame(t(sapply(score.groups, function(denom){
    n.denom <- length(denom)
    n.ovr <- length(which(genes %in% denom))
    frac <- n.ovr / n.genes
    binom.p <- binom.test(n.ovr, n.denom, p=baseline)$p.value
    return(c(n.ovr, n.denom, frac, binom.p))
  })))
  colnames(enrich) <- c("hits", "denom", "frac", "binom.p")
  enrich$bonf <- (enrich$binom.p <= 0.05 / length(score.groups))
  return(enrich)
}

# Compute enrichments of all gene sets vs. score bins
calc.enrichments <- function(gsets, scores, score, n.bins=10){
  score.groups <- bin.genes(scores, score, n.bins)
  enrich.df <- lapply(gsets, calc.enrichment.single, score.groups)
  names(enrich.df) <- names(gsets)
  return(enrich.df)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot enrichments as heatmap
plot.enrichments <- function(enrich.df, score, max.pct=0.5,
                             parmar=c(0.25, 8, 2.5, 0.25)){
  # Set score-dependent plotting parameters
  if(score == "pHI"){
    score.pal <- colorRampPalette(c("white", cnv.colors[1]))(101)
    bonf.col <- cnv.blacks[1]
  }else if(score == "pTS"){
    score.pal <- colorRampPalette(c("white", cnv.colors[2]))(101)
    bonf.col <- cnv.blacks[2]
  }
  
  # Collect plot parameters
  n.rows <- length(enrich.df)
  n.cols <- nrow(enrich.df[[1]])
  box.lefts <- 0:(n.cols-1)
  box.rights <- 1:n.cols
  x.ax.at <- seq(0, n.cols, length.out=6)
  x.ax.labels <- seq(0, 1, 0.2)

  # Prep plot area
  par(bty="n", mar=parmar)
  plot(NA, xlim=c(0, n.cols), ylim=c(n.rows, 0),
       xaxt="n", yaxt="n", xlab="", ylab="")
  
  # Add rectangles
  sapply(1:length(enrich.df), function(i){
    col.idxs <- (floor(100*enrich.df[[i]]$frac)+1)/max.pct
    col.idxs[which(col.idxs>=101)] <- 101
    rect(xleft=box.lefts, xright=box.rights,
         ybottom=i-1, ytop=i,
         border=bluewhite, col=score.pal[col.idxs])
    # rect(xleft=box.lefts[enrich.df[[i]]$bonf], 
    #      xright=box.rights[enrich.df[[i]]$bonf],
    #      ybottom=i-1, ytop=i,
    #      border=bonf.col, col=NA)
  })
  
  # Add labels
  axis(2, at=(1:n.rows)-0.5, tick=F, line=-1, las=2, cex=5.5/6, labels=names(enrich.df))
  
  # Add x-axis label
  axis(3, at=x.ax.at, labels=NA, tck=-0.025, col=blueblack)
  axis(3, at=x.ax.at, tick=F, line=-0.65, labels=x.ax.labels)
  mtext(3, line=1.5, text=paste("Genes per", score, "bin (%)"))
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
args <- parse_args(OptionParser(usage=paste("%prog scores.tsv genesets.tsv score", 
                                            "out_prefix"),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop(paste("Four positional arguments required: scores.tsv, genesets.tsv,", 
             "score, and out.pdf\n"))
}

# Writes args & opts to vars
scores.in <- args$args[1]
genesets.in <- args$args[2]
score <- args$args[3]
out.prefix <- args$args[4]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# scores.in <- "~/scratch/rCNV.gene_scores.tsv.gz"
# genesets.in <- "~/scratch/phi_vs_gene_sets.input.tsv"
# score <- "pHI"
# out.prefix <- "~/scratch/geneset_enrichments_test"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/gene_scores/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load scores
scores <- load.scores(scores.in)

# Load gene sets
gsets <- load.genesets(genesets.in, elig.genes=scores$gene)

# Stratify gene sets per score bin
enrich.df <- calc.enrichments(gsets, scores, score, n.bins=10)

# Plot enrichments
pdf(paste(out.prefix, score, "geneset_enrichments.pdf", sep="."),
    height=2, width=3)
plot.enrichments(enrich.df, score)
dev.off()
