#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot distributions of clustered CRBs for rCNV2 manuscript


options(stringsAsFactors=F, scipen=1000)


##########################
### PLOTTING FUNCTIONS ###
##########################
# Scatterplot of CRB size vs. number of elements
crb.scatter <- function(crbs, pt.cex=0.2, parmar=c(2.25, 3.25, 0.25, 0.25)){
  # Get plot data
  plot.df <- color.points.by.density(log10(crbs$size), log10(crbs$n_elements))
  xlims <- range(plot.df[, 1], na.rm=T)
  ylims <- range(plot.df[, 2], na.rm=T)
  
  # Prep plot area
  par(bty="n", mar=parmar)
  plot(NA, xlim=xlims, ylim=ylims, xaxt="n", yaxt="n", xlab="", ylab="")
  ax.at <- -10:10
  
  # Add background shading
  rect(xleft=par("usr")[1], xright=par("usr")[2], 
       ybottom=par("usr")[3], ytop=par("usr")[4],
       bty="n", border=NA, col=bluewhite)
  abline(h=log10(logscale.demi), v=ax.at, col="white")
  
  # Add points
  points(plot.df[, 1:2], pch=19, cex=pt.cex, col=plot.df$col)
  
  # Add axes
  axis(1, at=c(-10e10, 10e10), col=blueblack, tck=0, labels=NA)
  axis(1, at=ax.at, tck=-0.025, col=blueblack, labels=NA)
  sapply(1:length(logscale.major.bp), function(i){
    axis(1, at=log10(logscale.major.bp[i]), tick=F, line=-0.8, cex.axis=5.5/6, 
         labels=logscale.major.bp.labels[i])
  })
  mtext(1, line=1.2, text="CRB Size")
  axis(2, at=c(-10e10, 10e10), col=blueblack, tck=0, labels=NA)
  axis(2, at=log10(logscale.demi), tck=-0.025, col=blueblack, labels=NA)
  sapply(logscale.demi, function(y){
    axis(2, at=log10(y), tick=F, line=-0.65, las=2, cex.axis=5.5/6,
         labels=prettyNum(y, big.mark=","))
  })
  mtext(2, line=2.35, text="Elements in CRB")
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
args <- parse_args(OptionParser(usage=paste("%prog crbs.bed out.prefix"),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop(paste("Two positional arguments required: crbs.bed and out.prefix\n"))
}

# Writes args & opts to vars
crbs.in <- args$args[1]
out.prefix <- args$args[2]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# crbs.in <- "~/scratch/rCNV.crbs.bed.gz"
# out.prefix <- "~/scratch/crb_stats_test"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/noncoding_association/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load CRBs
crbs <- load.crbs(crbs.in)

# Print CRB stats for manuscript summary
print(paste("Mean CRB size:", round(mean(crbs$size)), "bp"))
print(paste("Median CRB size:", round(median(crbs$size)), "bp"))
print(paste("Mean elements per CRB:", round(mean(crbs$n_elements))))
print(paste("Median elements per CRB:", round(median(crbs$n_elements))))

# Scatterplot of CRB size vs number of elements in CRB
pdf(paste(out.prefix, "crb_stats.scatterplot.pdf", sep="."),
    height=2.75, width=2.9)
crb.scatter(crbs)
dev.off()

