#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot correlations between rCNV-derived gene scores and existing gene scores


# Set parameters & define constants
options(stringsAsFactors=F, scipen=1000, family="sans")
scores.to.keep <- c("pHI", "pTS", "gnomad_pLI", "gnomad_pRec", "gnomad_pNull",
                    "gnomad_oe_lof_upper", "gnomad_oe_mis_upper", "rvis_pct",
                    "exac_cnv_z", "hurles_hi", "eds")
score.abbrevs <- c("pHI" = "pHI",
                   "pTS" = "pTS",
                   "gnomad_pLI" = "pLI",
                   "gnomad_pRec" = "pRec",
                   "gnomad_pNull" = "pNull",
                   "gnomad_oe_lof_upper" = "LOEUF",
                   "gnomad_oe_mis_upper" = "Mis. OEUF",
                   "rvis_pct" = "RVIS Pct.",
                   "exac_cnv_z" = "CNV Z-Score",
                   "hurles_hi" = "HI Index",
                   "eds" = "EDS")
plot.height <- 1.9 #global plot height, in inches
png.res <- 400 #dpi resolution for png


######################
### DATA FUNCTIONS ###
######################
# Merge all rCNV-derived and existing scores into a single dataframe
merge.scores <- function(rcnv.scores, features, scores.to.keep){
  merged <- merge(rcnv.scores[, which(colnames(rcnv.scores) %in% c("gene", scores.to.keep))],
        features[, which(colnames(features) %in% c("gene", scores.to.keep))], 
        by="gene", all=T, sort=T)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot basic distribution for a single score
plot.score.hist <- function(scores, score, n.bins=100, blue.bg=TRUE,
                            parmar=c(1.1, 2.4, 1.1, 0.25)){
  # Extract plotting values
  vals <- as.numeric(scores[, score])
  vals <- vals[which(!is.na(vals))]
  xlim <- range(vals)
  pal <- viridis(n.bins)
  
  # Bin values
  breaks <- seq(xlim[1], xlim[2]+10e-6, length.out=n.bins+1)
  counts <- sapply(2:length(breaks), function(i){
    length(which(vals>=breaks[i-1] & vals<breaks[i]))
  })
  ymax <- 1.05*max(counts, na.rm=T)
  
  # Set parameters
  if(blue.bg==TRUE){
    plot.bg <- bluewhite
    plot.border <- NA
    plot.bty <- "n"
    grid.col <- "white"
  }else{
    plot.bg <- "white"
    plot.border <- NA
    plot.bty <- "n"
    grid.col <- NA
  }
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=xlim, ylim=c(0, ymax),
       xaxt="n", yaxt="n", xlab="", ylab="", yaxs="i")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       col=plot.bg, border=plot.border, bty=plot.bty)
  y.ax.at <- axTicks(2)
  abline(h=y.ax.at, col=grid.col)
  
  # Add rectangles
  rect(xleft=breaks[-length(breaks)], xright=breaks[-1],
       ybottom=0, ytop=counts, col=pal, border=pal)
  segments(x0=breaks, x1=breaks, y0=c(0, counts), y1=c(counts, 0), col=blueblack)
  segments(x0=breaks[-length(breaks)], x1=breaks[-1], y0=counts, y1=counts, col=blueblack)
  
  # Add axes
  axis(1, at=c(-10e10, 10e10), tck=0, col=blueblack, labels=NA)
  mtext(1, line=0.1, text=score.abbrevs[score])
  y.ax.lab <- y.ax.at / 100
  axis(2, at=c(-10e10, 10e10), tck=0, col=blueblack, labels=NA)
  axis(2, at=y.ax.at, tck=-0.025, col=blueblack, labels=NA)
  axis(2, at=y.ax.at, tick=F, line=-0.7, las=2, labels=y.ax.lab)
  mtext(2, line=1.4, text="Genes (x100)")
}

# Plot correlation of two scores
score.corplot <- function(scores, x.score, y.score, palette, pt.cex=0.15,
                          blue.bg=TRUE, parmar=c(1.1, 1.1, 1.1, 1.1)){
  # Get plot values
  plot.df <- data.frame("x"=scores[, x.score], "y"=scores[, y.score])
  plot.df <- plot.df[which(complete.cases(plot.df)), ]
  plot.df <- color.points.by.density(plot.df$x, plot.df$y, palette=palette)
  if(blue.bg==TRUE){
    plot.bg <- bluewhite
    plot.border <- NA
    plot.bty <- "n"
    grid.col <- "white"
  }else{
    plot.bg <- "white"
    plot.border <- NA
    plot.bty <- "n"
    grid.col <- NA
  }
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(plot.df[, 1:2], type="n", xaxt="n", yaxt="n", xlab="", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       col=plot.bg, border=plot.border, bty=plot.bty)
  x.tick.at <- seq(par("usr")[1], par("usr")[2], length.out=6)
  y.tick.at <- seq(par("usr")[3], par("usr")[4], length.out=6)
  abline(v=x.tick.at, h=y.tick.at, col=grid.col)
  
  # Add points
  points(plot.df[, 1:2], pch=19, col=plot.df$col, cex=pt.cex)
  
  # Add axes
  axis(1, at=c(-10e10, 10e10), tck=0, col=blueblack, labels=NA)
  mtext(1, line=0.1, text=score.abbrevs[x.score])
  axis(2, at=c(-10e10, 10e10), tck=0, col=blueblack, labels=NA)
  mtext(2, line=0.1, text=score.abbrevs[y.score])
  
  # Add correlation coefficients
  r <- format(round(cor(plot.df$x, plot.df$y, method="pearson"), 2), nsmall=2)
  rho <- format(round(cor(plot.df$x, plot.df$y, method="spearman"), 2), nsmall=2)
  tau <- format(round(cor(plot.df$x, plot.df$y, method="kendall"), 2), nsmall=2)
  mtext(3, line=-0.1, cex=5/6,
        text=bquote("r" == .(r) * ";" ~ rho == .(rho) * ";" ~ tau == .(tau)))
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
args <- parse_args(OptionParser(usage=paste("%prog scores.tsv features.bed metadata.tsv out.prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop(paste("Four positional arguments required: scores.tsv, features.bed, feature_metadata.tsv, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
scores.in <- args$args[1]
features.in <- args$args[2]
feature.metadata.in <- args$args[3]
out.prefix <- args$args[4]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# scores.in <- "~/scratch/rCNV.gene_scores.tsv.gz"
# features.in <- "~/scratch/gencode.v19.canonical.pext_filtered.all_features.no_variation.bed.gz"
# feature.metadata.in <- "~/scratch/gene_feature_metadata.tsv"
# out.prefix <- "~/scratch/gene_score_correlations"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/gene_scores/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load scores, features, and feature metadata
rcnv.scores <- load.scores(scores.in)
features <- load.features(features.in, fill=NA, norm=F)
feat.meta <- load.gene.feature.metadata(feature.metadata.in)

# Merge all scores into single dataframe
scores <- merge.scores(rcnv.scores, features, scores.to.keep)

# Plot basic histogram of each score
sapply(scores.to.keep, function(score){
  pdf(paste(out.prefix, score, "hist.pdf", sep="."), height=plot.height, width=2.15)
  plot.score.hist(scores, score, n.bins=40, blue.bg=FALSE)
  dev.off()
})

# Set color palettes
del.pal <- colorRampPalette(c(cnv.blacks[1], cnv.colors[1], 
                              control.cnv.colors[1], cnv.whites[1]))(256)
dup.pal <- colorRampPalette(c(cnv.blacks[2], cnv.colors[2], 
                              control.cnv.colors[2], cnv.whites[2]))(256)

# Plot correlations of all scores vs pHI
sapply(scores.to.keep, function(score){
  # pdf(paste(out.prefix, score, "vs_pHI.pdf", sep="."), 
  #     height=plot.height, width=plot.height)
  png(paste(out.prefix, score, "vs_pHI.png", sep="."), 
      height=png.res*plot.height, width=png.res*plot.height, res=png.res)
  score.corplot(scores, score, "pHI", del.pal, blue.bg=FALSE)
  dev.off()
})

# Plot correlations of all scores vs pHI
sapply(scores.to.keep, function(score){
  # pdf(paste(out.prefix, score, "vs_pTS.pdf", sep="."), 
  #     height=plot.height, width=plot.height)
  png(paste(out.prefix, score, "vs_pTS.png", sep="."), 
      height=png.res*plot.height, width=png.res*plot.height, res=png.res)
  score.corplot(scores, score, "pTS", dup.pal, blue.bg=FALSE)
  dev.off()
})

