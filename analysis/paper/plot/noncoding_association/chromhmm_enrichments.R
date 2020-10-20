#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Compute & plot ChromHMM enrichments from genome annotation burden tests for rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load ChromHMM state manifest from .tsv
load.chromhmm.manifest <- function(chromhmm.manifest.in){
  mfst <- read.table(chromhmm.manifest.in, sep="\t", comment.char="", header=T)
  mfst$color <- sapply(mfst$COLOR.CODE, function(str){
    cvals <- as.numeric(unlist(strsplit(str, split=",")))
    rgb(cvals[1], cvals[2], cvals[3], maxColorValue=255)
  })
  mfst$COLOR.CODE <- NULL
  colnames(mfst) <- tolower(colnames(mfst))
  mfst$mnemonic <- gsub("/", "", mfst$mnemonic)
  mfst$description <- gsub("enhancer1", "Enhancer 1", mfst$description)
  mfst$description <- gsub("enhancer2", "Enhancer 2", mfst$description)
  mfst$description <- gsub("transcription", "Transcription", mfst$description)
  mfst$description <- gsub("genes & repeats", "Genes & Repeats", mfst$description)
  return(mfst)
}

# Compute proportion of tissues by significance/direction for a single ChromHMM state
calc.single.chmm.enrichment <- function(stats, cnv, state){
  elig.idx <- which(stats$cnv==cnv & grepl(state, stats$trackname, fixed=T))
  ctrl.idx <- which(stats$meta.lnOR < 0)
  case.idx <- which(stats$meta.lnOR >= 0)
  sig.idx <- which(stats$meta.phred_p >= -log10(0.05))
  n.tracks <- length(elig.idx)
  frac.ctrl.sig <- length(intersect(intersect(elig.idx, ctrl.idx), sig.idx)) / n.tracks
  frac.ctrl.ns <- length(setdiff(intersect(elig.idx, ctrl.idx), sig.idx)) / n.tracks
  frac.case.ns <- length(setdiff(intersect(elig.idx, case.idx), sig.idx)) / n.tracks
  frac.case.sig <- length(intersect(intersect(elig.idx, case.idx), sig.idx)) / n.tracks
  return(c(frac.ctrl.sig, frac.ctrl.ns, frac.case.ns, frac.case.sig))
}

# Compute all case:control enrichments for all ChromHMM states for a single CNV type
calc.all.chmm.enrichments <- function(stats, cnv, states){
  enrich <- as.data.frame(t(sapply(1:nrow(states), function(i){
    state <- paste(states[i, 1:2], collapse="_")
    calc.single.chmm.enrichment(stats, cnv, state)
  })))
  colnames(enrich) <- c("ctrl.sig", "ctrl.ns", "case.ns", "case.sig")
  neworder <- with(enrich, order(case.sig + case.ns))
  enrich <- enrich[neworder, ]
  rownames(enrich) <- states$mnemonic[neworder]
  return(enrich)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot ordered enrichment of rCNV burden for ChromHMM tracks across tissues
plot.enrich <- function(stats, cnv, states, parmar=c(0.25, 8.75, 2.25, 1)){
  # Get plot data
  enrich <- calc.all.chmm.enrichments(stats, cnv, states)
  colors <- sapply(rownames(enrich), function(x){states$color[which(states$mnemonic==x)]})
  n.rows <- nrow(enrich)
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(-0.1, 1), ylim=c(0, n.rows),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  
  # Add background markings
  x.ax.at <- seq(0, 1, 0.25)
  abline(v=x.ax.at, col=bluewhite)
  abline(v=0.5, col=blueblack)
  
  # Add bars per state
  sapply(1:nrow(enrich), function(i){
    # Get state-specific values
    xvals <- c(0, cumsum(as.numeric(enrich[i, ])))
    yvals <- c(i-0.85, i-0.15)
    s.color <- colors[i]
    rcols <- as.character(c(ns.color.dark, ns.color, control.cnv.colors[cnv], cnv.colors[cnv]))
    border.col <- cnv.blacks[cnv]
    
    # Add Y-axis legend & label
    axis(2, at=i-0.5, tick=F, line=-1, las=2, cex.axis=5/6,
         labels=states$description[which(states$mnemonic==rownames(enrich)[i])])
    rect(xleft=-0.075, xright=-0.025, ybottom=yvals[1], ytop=yvals[2],
         col=s.color, border=blueblack)
    
    # Plot bars
    rect(xleft=xvals[1:4], xright=xvals[2:5],
         ybottom=yvals[1], ytop=yvals[2],
         col=rcols, border=rcols, lwd=0.5)
    segments(x0=xvals[3], x1=xvals[3], y0=yvals[1], y1=yvals[2],
             lend="butt", col=cnv.blacks[cnv])
    rect(xleft=0, xright=par("usr")[2],
         ybottom=yvals[1], ytop=yvals[2],
         col=NA, border=blueblack)
  })
  
  # Add axes
  abline(v=0, col=blueblack)
  axis(4, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(3, at=c(0, 10e10), col=blueblack, tck=0, labels=NA)
  axis(3, at=x.ax.at, tck=-0.025, labels=NA, col=blueblack)
  sapply(x.ax.at, function(x){
    axis(3, at=x, tick=F, line=-0.75, cex.axis=5/6, 
         labels=paste(round(100*x), "%", sep=""))
  })
  axis(3, at=0.5, tick=F, line=0.2, labels="Tissues in Roadmap")
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
args <- parse_args(OptionParser(usage=paste("%prog stats.tsv  REP_states.tsv out.prefix"),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop(paste("Three positional arguments required: stats.tsv, REP_states.tsv, and out.prefix\n"))
}

# Writes args & opts to vars
stats.in <- args$args[1]
chromhmm.manifest.in <- args$args[2]
out.prefix <- args$args[3]
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# stats.in <- "~/scratch/rCNV.burden_stats.tsv.gz"
# chromhmm.manifest.in <- "~/scratch/REP_state_manifest.tsv"
# out.prefix <- "~/scratch/chromhmm_enrich_test"
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

# Load ChromHMM manifest
states <- load.chromhmm.manifest(chromhmm.manifest.in)

# Plot enrichment separately for DEL & DUP
sapply(c("DEL", "DUP"), function(cnv){
  pdf(paste(out.prefix, "chromhmm_enrichment", cnv, "pdf", sep="."),
      height=3.1, width=3.7)
  plot.enrich(stats, cnv, states)
  dev.off()
})
