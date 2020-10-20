#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Volcano plots of track-level rCNV association tests for rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


##########################
### PLOTTING FUNCTIONS ###
##########################
# Volcano plot for a single CNV type
plot.volcano <- function(stats, cnv, pt.cex=0.2, parmar=c(2.25, 2.25, 0.25, 0.25)){
  # Get plot data
  plot.df <- stats[which(stats$cnv==cnv), c("meta.lnOR", "meta.phred_p")]
  xlims <- range(stats$meta.lnOR, na.rm=T)
  ylims <- range(stats$meta.phred_p, na.rm=T)
  x.label <- paste(c("DEL"="Deletion", "DUP"="Duplication")[cnv], "Odds Ratio")
  
  # Assign point colors
  plot.df$color <- apply(plot.df, 1, function(vals){
    sig <- (vals[2] > -log10(0.05))
    if(vals[1] > 0){
      if(sig==TRUE){
        return(cnv.colors[cnv])
      }else{
        return(control.cnv.colors[cnv])
      }
    }else{
      if(sig==TRUE){
        return(ns.color.dark)
      }else{
        return(ns.color)
      }
    }
  })
  
  # Prep plot area
  par(bty="n", mar=parmar)
  plot(NA, xlim=xlims, ylim=ylims, xaxt="n", yaxt="n", xlab="", ylab="")
  x.ax.at <- log(c(1/(2^(4:0)), 2^(1:4)))
  y.ax.at <- -10:10
  
  # Add background shading
  rect(xleft=par("usr")[1], xright=par("usr")[2], 
       ybottom=par("usr")[3], ytop=par("usr")[4],
       bty="n", border=NA, col=bluewhite)
  abline(h=y.ax.at, v=x.ax.at, col="white")
  
  # Add points & nomsig marker
  points(plot.df[, 1:2], pch=19, cex=pt.cex, col=plot.df$color)
  abline(h=-log10(0.05), lty=5, col=cnv.blacks[cnv])
  
  # Add text to top corners noting number of nomsig tracks
  n.nomsig.ctrl <- length(which(plot.df$color==ns.color.dark))
  n.nomsig.case <- length(which(plot.df$color==cnv.colors[cnv]))
  y.buf <- 0.075*diff(par("usr")[3:4])
  text(x=par("usr")[1], y=par("usr")[4]-y.buf, col=ns.color.dark, pos=4,
       labels=prettyNum(n.nomsig.ctrl))
  text(x=par("usr")[2], y=par("usr")[4]-y.buf, col=cnv.colors[cnv], pos=2,
       labels=prettyNum(n.nomsig.case))
  
  # Add axes
  axis(1, at=c(-10e10, 10e10), col=blueblack, tck=0, labels=NA)
  axis(1, at=x.ax.at, tck=-0.025, col=blueblack, labels=NA)
  sapply(x.ax.at, function(x){
    axis(1, at=x, tick=F, line=-0.8, cex.axis=5.5/6, labels=exp(x))
  })
  mtext(1, line=1.2, text=x.label)
  axis(2, at=c(-10e10, 10e10), col=blueblack, tck=0, labels=NA)
  axis(2, at=y.ax.at, tck=-0.025, col=blueblack, labels=NA)
  axis(2, at=y.ax.at, tick=F, line=-0.65, las=2, cex.axis=5.5/6)
  mtext(2, line=1, text=bquote(-log[10](italic(P)[.(cnv)])))
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
args <- parse_args(OptionParser(usage=paste("%prog stats.tsv out.prefix"),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop(paste("Two positional arguments required: stats.tsv and out.prefix\n"))
}

# Writes args & opts to vars
stats.in <- args$args[1]
out.prefix <- args$args[2]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# stats.in <- "~/scratch/rCNV.burden_stats.tsv.gz"
# out.prefix <- "~/scratch/volcano_test"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/noncoding_association/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load track stats
stats <- load.track.stats(stats.in)

# Generate one volcano plot each for DEL & DUP
sapply(c("DEL", "DUP"), function(cnv){
  pdf(paste(out.prefix, "track_volcano", cnv, "pdf", sep="."),
      height=2.25, width=2.25)
  plot.volcano(stats, cnv)
  dev.off()
})

