#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot grid summarizing large segment association across phenotypes


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load flat table of effect sizes and P-values per region per phenotype per CNV
load.sumstats <- function(sumstats.in){
  sumstats <- read.table(sumstats.in, header=T, sep="\t")
  sumstats[, 4:8] <- apply(sumstats[, 4:8], 2, as.numeric)
  return(sumstats)
}

# Load list of clusters to plot on a single line
load.clusters <- function(clusters.in, loci){
  if(!is.null(clusters.in)){
    clusters <- strsplit(read.table(clusters.in, header=F, sep="\t")[, 1], split=",")
  }else{
    clusters <- lapply(loci$region_id, function(x){x})
  }
  names(clusters) <- unlist(lapply(clusters, function(x){sort(sapply(strsplit(x, split="_"), function(l){l[4]}))[1]}))
  return(clusters)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot a single semicircle
plot.semicircle <- function(x, y, r, side="top", color="black", border=NA){
  # Adopted from https://stackoverflow.com/questions/31538534/plotting-half-circles-in-r
  rs <- seq(0, pi, len=100)
  xm <- r*cos(rs)
  ym <- r*sin(rs)
  if(side == "bottom"){
    sign <- -1
  }else{
    sign <- 1
  }
  xc <- x + (sign * xm)
  yc <- y + (sign * ym)
  polygon(xc, yc, col=color, border=border, lwd=1)
}

# Plot data for a single locus-phenotype pair
plot.locus.single_hpo <- function(loci, sumstats, region_ids, hpo, x, y, 
                                  stat.size="pvalue", max.stat.size=8,
                                  stat.color="lnor", max.stat.color=4){
  # Scale radii
  stat.size.idx <- which(colnames(sumstats) == stat.size)
  stat.color.idx <- which(colnames(sumstats) == stat.color)
  radii <- as.numeric(sapply(c("DUP", "DEL"), function(cnv){
    radius.x <- max(sumstats[which(sumstats$region_id %in% region_ids & sumstats$hpo==hpo & sumstats$cnv==cnv), stat.size.idx], na.rm=T)
    min(c(max(c(0, radius.x), na.rm=T), max.stat.size)) / (2 * max.stat.size)
  }))
  # Determine color based on stat.color
  colors <- sapply(c("DUP", "DEL"), function(cnv){
    color.x <- max(sumstats[which(sumstats$region_id %in% region_ids & sumstats$hpo==hpo & sumstats$cnv==cnv), stat.color.idx], na.rm=T)
    cnv.color.palettes[[cnv]][round(min(c(max(c(0, color.x), na.rm=T), max.stat.color)) * 100 / max.stat.color) + 1]
  })
  # Determine border based on significance
  sig <- sapply(c("DUP", "DEL"), function(cnv){
    if(any(hpo %in% unlist(loci$hpos[which(loci$region_id %in% region_ids & loci$cnv==cnv)]))){
      "black"
    }else{
      NA
    }
  })
  # Plot semicircles
  plot.semicircle(x, y, r=radii[1], side="top", color=colors[1], border=sig[1])
  plot.semicircle(x, y, r=radii[2], side="bottom", color=colors[2], border=sig[2])
}

# Plot a row of all hpos per region
plot.locus.row <- function(loci, sumstats, hpos, region_ids, y, 
                           stat.size="pvalue", max.stat.size=8, 
                           stat.color="lnor", max.stat.color=4,
                           outline=TRUE, outline.color="gray80", buffer=0.15){
  x.at <- sapply(1:length(hpos), function(x){x-0.5+(buffer * (x-1))})
  if(outline==TRUE){
    # segments(x0=0, x1=length(hpos), y0=y, y1=y, col=outline.color)
    sapply(1:length(hpos), function(x){
      draw.circle(x.at[x], y, radius=0.5, border=outline.color, col="white")
    })
  }
  sapply(1:length(hpos), function(x){
    plot.locus.single_hpo(loci, sumstats, region_ids, hpos[x], x.at[x], y, 
                          stat.size, max.stat.size, stat.color, max.stat.color)
  })
}

# Plot a grid of all regions
plot.all.loci <- function(loci, clusters, sumstats, hpos, 
                          stat.size="pvalue", max.stat.size=8, 
                          stat.color="lnor", max.stat.color=4,
                          outline=TRUE, outline.color="gray60", 
                          background=TRUE, background.color="gray80",
                          buffer=0.15){
  # Gather plot parameters
  nhpos <- length(hpos)
  x.cols <- sapply(1:length(hpos), function(i){i - 0.5 + (buffer * (i-1))})
  y.rows <- sapply(1:length(clusters), function(i){-i + 0.5 - (buffer * (i-1))})
  
  # Prep plot area
  par(mar=c(0.5, 6, 8, 0.5), bty="n")
  plot(x=c(0, ceiling(max(x.cols))), y=c(0, floor(min(y.rows))), type="n", asp=1, 
       xaxt="n", yaxt="n", xlab="", ylab="")
  axis(2, at=y.rows, tick=F, las=2, line=-0.9, labels=names(clusters))
  axis(3, at=x.cols, tick=F, las=2, line=-0.8,
       labels=sapply(hpos, function(hpo){hpo.abbrevs[which(names(hpo.abbrevs)==hpo)]}))
  # text(x=x.cols - 0.75, y=par("usr")[4], srt=45, pos=4, xpd=T, 
  #      labels=sapply(hpos, function(hpo){hpo.abbrevs[which(names(hpo.abbrevs)==hpo)]}))
  if(background==TRUE){
    bg.col <- adjustcolor(background.color, alpha=1/3)
    rect(xleft=min(x.cols - 0.35), xright=max(x.cols + 0.35),
         ybottom=min(y.rows - 0.35), ytop=max(y.rows + 0.35),
         border=NA, bty="n", col=bg.col)
    rect(xleft=x.cols - 0.35, xright=x.cols + 0.35,
         ybottom=min(y.rows), ytop=par("usr")[4],
         border=NA, bty="n", col=bg.col)
    rect(xleft=par("usr")[1], xright=max(x.cols),
         ybottom=y.rows-0.35, ytop=y.rows+0.35,
         border=NA, bty="n", col=bg.col)
  }
  
  # Plot each region
  sapply(1:length(clusters), function(ridx){
    region_ids <- clusters[[ridx]]
    plot.locus.row(loci, sumstats, hpos, region_ids, y.rows[ridx], 
                   stat.size, max.stat.size, stat.color, max.stat.color,
                   outline, outline.color, buffer)
  })
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(plotrix, quietly=T)
require(funr, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--clusters"), help="Comma-delimited list of regions to collapse."),
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog loci.bed sumstats.tsv hpos.tsv", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop(paste("Four positional arguments required: loci.bed, sumstats.tsv, hpos.tsv, output.pdf\n", sep=" "))
}

# Writes args & opts to vars
loci.in <- args$args[1]
sumstats.in <- args$args[2]
hpos.in <- args$args[3]
outfile <- args$args[4]
clusters.in <- opts$`clusters`
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# sumstats.in <- "~/scratch/rCNV.final_segments.loci.all_sumstats.tsv.gz"
# hpos.in <- "~/scratch/rCNV2_analysis_d1.reordered_hpos.txt"
# clusters.in <- "~/scratch/locus_clusters.txt"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# outfile <- "~/scratch/test.pdf"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Read data
loci <- load.loci(loci.in)
sumstats <- load.sumstats(sumstats.in)
hpos <- as.character(read.table(hpos.in, header=F, sep="\t")[, 1])
clusters <- load.clusters(clusters.in, loci)

# Plot locus grid
pdf(outfile, height=8, width=8)
plot.all.loci(loci, clusters, sumstats, hpos, 
              stat.size="pvalue", max.stat.size=6, 
              stat.color="lnor", max.stat.color=log(8),
              outline=T, outline.color="gray75", background=T, buffer=0.1)
dev.off()
