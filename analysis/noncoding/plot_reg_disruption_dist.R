#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot distributions of regulatory disruption scores for noncoding CNVs by cohort


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load & parse regulatory disruption scoring results
load.reg.scores <- function(in.dir, cohorts, cnv){
  scores <- lapply(cohorts, function(cohort){
    dat <- lapply(c("developmental", "adult", "controls"), function(pheno){
      in.path <- paste(in.dir, paste(cohort, "reg_info", cnv, pheno, "tsv.gz", sep="."), sep="/")
      read.table(in.path, header=T, sep="\t")
    })
    names(dat) <- c("dev", "adult", "control")
    return(dat)
  })
  names(scores) <- cohorts
  return(scores)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot predicted regulatory effects
plot.reg.scores <- function(scores, cnv, variable, title=NULL, cnv.prefix="Noncoding",
                            parmar=c(0.5, 3, 3, 0.5)){
  # Get plot data
  pdat <- lapply(scores, function(l){lapply(l, function(l){
    if(nrow(l) > 0){
      as.numeric(l[, variable])
    }else{
      return(0)
    }
    })})
  ylims <- range(unlist(unlist(pdat)), na.rm=T)
  xlims <- c(0, length(pdat))
  group.colors <- c(cnv.colors[cnv], control.cnv.colors[cnv], ns.color)

  # Prep plot area
  par(bty="n", mar=parmar)
  plot(NA, xlim=xlims, ylim=ylims,
       xaxt="n", xaxs="i", xlab="",
       yaxt="n", ylab="")

  # Add vioplots, lines, and labels for each cohort
  sapply(1:length(pdat), function(i){
    # Determine cohort-specific outlier cutoff
    outlier.cutoff <- quantile(abs(unlist(pdat[[i]])), probs=0.95)
    # Violins for the middle 99% of data
    vioplot(lapply(pdat[[i]], function(v){v[which(abs(v) < outlier.cutoff)]}),
            at=i-c(0.75, 0.5, 0.25), add=T,
            col=group.colors, drawRect=FALSE, wex=0.25)
    # Points for the top 1% of all data
    sapply(1:3, function(k){
      kvals <- pdat[[i]][[k]]
      points(x=rep(i-(1-(0.25*k)), length(which(abs(kvals) >= outlier.cutoff))),
             y=kvals[which(abs(kvals) >= outlier.cutoff)],
             pch=1, cex=0.2, lwd=0.75, xpd=T,
             col=adjustcolor(group.colors[k], alpha=0.4))
      # beeswarm(kvals[which(abs(kvals) >= outlier.cutoff)],
      #          at=i-(1-(0.25*k)), add=T, corral="wrap", corralWidth=0.25,
      #          pch=19, col=group.colors[k], cex=0.3)
    })
    # Axis and title
    axis(3, at=i-c(0.1, 0.9), tck=0, labels=NA)
    axis(3, at=i-0.5, line=-0.9, cex.axis=0.9, tick=F, labels=cohort.names[i])
  })

  # Add axes & labels
  axis(2, tck=0, labels=NA, at=c(-10e10, 10e10))
  y.ax.at <- axTicks(2)
  axis(2, at=y.ax.at, tck=-0.015, labels=NA)
  axis(2, at=y.ax.at, tick=F, las=2, cex.axis=0.85, line=-0.65)
  mtext(2, line=2, text=title)
  mtext(3, line=2, font=2,
        text=paste(title, "for", cnv.prefix,
                   c("DEL" = "Deletions", "DUP" = "Duplications")[cnv]))

  # Add legend
  lside <- if(cnv=="DEL"){"bottomright"}else{"topright"}
  legend(lside, cex=0.75, bg="white",
         legend=c("Developmental Case", "Adult-Onset Case", "Control"),
         fill=group.colors, xpd=T)
}



#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(rCNV2, quietly=T)
require(vioplot, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--cnv"), default='DEL', metavar='character',
              help="CNV type to plot [default: DEL]"),
  make_option(c("--prefix"), default='Noncoding', metavar='character',
              help="CNV prefix for title of plot [default: Noncoding]")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog in_dir out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 2){
  stop(paste("Three positional arguments required: in_dir, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
in.dir <- args$args[1]
out.prefix <- args$args[2]
cnv <- opts$cnv
cnv.prefix <- opts$prefix

# # DEV PARAMETERS
# in.dir <- "~/scratch/reg_scoring_results"
# out.prefix <- "~/scratch/reg_test"
# cnv <- "DEL"

# Load reg scores per cohort split by phenotype
scores <- load.reg.scores(in.dir, names(cohort.abbrevs), cnv)

# Plot predicted expression
pdf(paste(out.prefix, "pred_exp.pdf", sep="."),
    height=4, width=8)
plot.reg.scores(scores, cnv, "pred_exp",
                title="Predicted Expression Z-scores")
dev.off()

# Plot regulatory disruption scores
pdf(paste(out.prefix, "reg_dist.pdf", sep="."),
    height=4, width=8)
plot.reg.scores(scores, cnv, "reg_dist",
                title="Regulatory Disruption Score",
                cnv.prefix=cnv.prefix)
dev.off()
