#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot mini version of pHI & pTS scores for graphical abstract of rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
load.scores <- function(scores.in){
  scores <- read.table(scores.in, header=T, sep="\t", comment.char="", check.names=F)
  colnames(scores)[1] <- gsub("#", "", colnames(scores)[1], fixed=T)
  scores[, -1] <- apply(scores[, -1], 2, as.numeric)
  scores$color <- apply(scores[, -1], 1, function(vals){
    if(vals[1]>=0.9 & vals[2]<0.9){
      cnv.colors[which(names(cnv.colors) == "DEL")]
    }else if(vals[1]<0.9 & vals[2]>=0.9){
      cnv.colors[which(names(cnv.colors) == "DUP")]
    }else if(vals[1]>=0.9 & vals[2]>=0.9){
      cnv.colors[which(names(cnv.colors) == "CNV")]
    }else{
      "gray70"
    }
  })
  return(scores)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
plot.mini.scatter <- function(scores, n.points=1000, seed=2020){
  # Downsample points
  set.seed(seed)
  scores <- scores[sample(1:nrow(scores), n.points, replace=F), ]
  # Prep plot area
  par(mar=c(1.1, 1.1, 0.2, 0.2), bty="n", bg=NA)
  plot(NA, xlim=c(0, 1), ylim=c(0, 1),
       xaxs="i", xlab="", xaxt="n",
       yaxs="i", ylab="", yaxt="n")
  # Add shading
  rect(xleft=par("usr")[1], xright=0.9,
       ybottom=par("usr")[3], ytop=0.9,
       bty="n", border=NA, col=bluewhite)
  rect(xleft=par("usr")[1], xright=0.9,
       ybottom=0.9, ytop=par("usr")[4],
       bty="n", border=NA, col=control.cnv.colors[2])
  rect(xleft=0.9, xright=par("usr")[2],
       ybottom=0.9, ytop=par("usr")[4],
       bty="n", border=NA, col=control.cnv.colors[3])
  rect(xleft=0.9, xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=0.9,
       bty="n", border=NA, col=control.cnv.colors[1])
  # abline(h=seq(0, 1, 0.2), v=seq(0, 1, 0.2), col="white")
  # Add points
  points(x=scores$pHI, y=scores$pTS, pch=19, cex=0.25, 
         col=scores$color)
  # Axes
  # axis(2, labels=NA, tck=-0.03, col=blueblack)
  # axis(3, labels=NA, tck=-0.03, col=blueblack)
  mtext(2, line=0.1, text="Duplication", col=blueblack)
  mtext(1, line=0.1, text="Deletion", col=blueblack)
  # Cleanup
  abline(h=0.9, col=blueblack, lty=2)
  abline(v=0.9, col=redblack, lty=2)
  box(bty="o", xpd=T, col=blueblack)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog scores.tsv out.png", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop(paste("Two positional arguments required: scores.tsv and output.png\n", sep=" "))
}

# Writes args & opts to vars
scores.in <- args$args[1]
out.png <- args$args[2]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# scores.in <- "~/scratch/rCNV.gene_scores.tsv.gz"
# out.png <- "~/scratch/test_mini_scatter.png"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Load scores
scores <- load.scores(scores.in)

# Plot scores
png(out.png, height=1.15*300, width=1.15*300, res=300, family="sans")
plot.mini.scatter(scores, n.points=1250)
dev.off()
