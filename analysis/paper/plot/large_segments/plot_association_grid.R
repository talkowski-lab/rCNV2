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
# Load summary dataframe for loci
load.loci <- function(loci.in){
  # Read data
  loci <- read.table(loci.in, sep="\t", header=T, comment.char="")
  colnames(loci)[1] <- "chr"
  
  # Split list-style columns
  loci$hpos <- strsplit(loci$hpos, split=";")
  loci$constituent_assocs <- strsplit(loci$constituent_assocs, split=";")
  loci$cred_interval_coords <- strsplit(loci$cred_interval_coords, split=";")
  loci$genes <- strsplit(loci$genes, split=";")
  
  # Convert numeric columns to numerics
  numeric.cols <- c("start_min", "end_max", "pooled_control_freq", "pooled_case_freq",
                    "pooled_ln_or", "pooled_ln_or_ci_lower", "pooled_ln_or_ci_upper", 
                    "min_ln_or", "max_ln_or", "n_hpos", "n_constituent_assocs", 
                    "n_cred_intervals", "cred_intervals_size", "n_genes")
  for(col in numeric.cols){
    cidx <- which(colnames(loci) == col)
    loci[, cidx] <- as.numeric(loci[, cidx])
  }
  
  # Add formatted locus names and sizes
  loci$name <- sapply(strsplit(loci$region_id, split="_"), function(parts){parts[4]})
  loci$size <- paste(prettyNum(round(loci$cred_intervals_size/1000, 0), big.mark=","), "kb", sep=" ")
  
  return(loci)
}

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
plot.semicircle <- function(x, y, r, side="top", color="black"){
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
  polygon(xc, yc, col=color, border=NA)
}

# Plot data for a single locus-phenotype pair
plot.locus.single_hpo <- function(loci, sumstats, region_ids, hpo, x, y, stat="pvalue", max.stat=10){
  # Scale radii
  stat.idx <- which(colnames(sumstats) == stat)
  radii <- as.numeric(sapply(c("DUP", "DEL"), function(cnv){
    x <- max(sumstats[which(sumstats$region_id %in% region_ids & sumstats$hpo==hpo & sumstats$cnv==cnv), stat.idx], na.rm=T)
    min(c(max(c(0, x), na.rm=T), max.stat)) / (2 * max.stat)
  }))
  # Determine color based on significance
  colors <- sapply(c("DUP", "DEL"), function(cnv){
    if(any(hpo %in% unlist(loci$hpos[which(loci$region_id %in% region_ids & loci$cnv==cnv)]))){
      as.character(cnv.colors[which(names(cnv.colors)==cnv)])
    }else{
      as.character(control.cnv.colors[which(names(control.cnv.colors)==cnv)])
    }
  })
  # Plot semicircles
  plot.semicircle(x, y, r=radii[1], side="top", color=colors[1])
  plot.semicircle(x, y, r=radii[2], side="bottom", color=colors[2])
}

# Plot a row of all hpos per region
plot.locus.row <- function(loci, sumstats, hpos, region_ids, y, 
                           stat="pvalue", max.stat=10, outline=TRUE,
                           outline.color="gray80", buffer=0.15){
  x.at <- sapply(1:length(hpos), function(x){x-0.5+(buffer * (x-1))})
  if(outline==TRUE){
    # segments(x0=0, x1=length(hpos), y0=y, y1=y, col=outline.color)
    sapply(1:length(hpos), function(x){
      draw.circle(x.at[x], y, radius=0.5, border=outline.color, col="white")
    })
  }
  sapply(1:length(hpos), function(x){
    plot.locus.single_hpo(loci, sumstats, region_ids, hpos[x], x.at[x], y, stat, max.stat)
  })
}

# Plot a grid of all regions
plot.all.loci <- function(loci, clusters, sumstats, hpos, stat="pvalue", max.stat=10, 
                          outline=TRUE, outline.color="gray60", 
                          background=TRUE, background.color="gray80",
                          buffer=0.15){
  # Gather plot parameters
  nhpos <- length(hpos)
  x.cols <- sapply(1:length(hpos), function(i){i - 0.5 + (buffer * (i-1))})
  y.rows <- sapply(1:length(clusters), function(i){-i + 0.5 - (buffer * (i-1))})
  
  # Prep plot area
  par(mar=c(0.5, 6, 6, 0.5), bty="n")
  plot(x=c(0, ceiling(max(x.cols))), y=c(0, floor(min(y.rows))), type="n", asp=1, 
       xaxt="n", yaxt="n", xlab="", ylab="")
  axis(2, at=y.rows, tick=F, las=2, line=-0.9, labels=names(clusters))
  text(x=x.cols - 0.75, y=par("usr")[4], srt=45, pos=4, xpd=T, 
       labels=sapply(hpos, function(hpo){hpo.abbrevs[which(names(hpo.abbrevs)==hpo)]}))
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
                   stat, max.stat, outline, outline.color, buffer)
  })
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse)
require(plotrix)

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
  stop(paste("Seven positional arguments: loci.bed, sumstats.tsv, hpos.tsv, output.pdf\n", sep=" "))
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
# hpos.in <- "~/scratch/hpos.txt"
# clusters.in <- "~/scratch/locus_clusters.txt"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# outfile <- "~/scratch/test.pdf"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Read data
loci <- load.loci(loci.in)
sumstats <- load.sumstats(sumstats.in)
hpos <- as.character(read.table(hpos.in, header=F, sep="\t")[, 1])
clusters <- load.clusters(clusters.in, loci)

# Plot locus grid
pdf(outfile, height=8, width=8)
plot.all.loci(loci, clusters, sumstats, hpos, stat="pvalue", max.stat=8, 
              outline=T, outline.color="gray75", background=T, buffer=0.1)
dev.off()