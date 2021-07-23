#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Optimize CDS overlap from precomputed meta-analysis stats
# (Subroutine of optimize_cds_overlap.sh)


# Load libraries
require(rCNV2, quietly=TRUE)
require(optparse, quietly=T)


# Compute CDS optimization data and determine optimal cutoff
optimize.cds <- function(stats){
  # Scale primary P-value and odds ratios as a fraction of maximum possible
  stats$frac_p <- stats$meta_phred_p / max(stats$meta_phred_p, na.rm=T)
  stats$frac_or <- stats$meta_lnOR / max(stats$meta_lnOR, na.rm=T)

  # Compute Euclidean distance of scaled P & OR from best possible value (1, 1)
  stats$pareto_dist <- sqrt(((stats$frac_p-1)^2) + ((stats$frac_or-1)^2))

  # Prepare optimization data frame
  opt.df <- data.frame("min.cds" = stats$min_cds,
                       "power.frac" = stats$frac_p,
                       "or.frac" = stats$frac_or,
                       "pareto.dist" = stats$pareto_dist)

  # Print optimal cutoff
  best.idx <- which(opt.df$pareto.dist == min(opt.df$pareto.dist, na.rm=T))
  cat(paste("\nOptimal cutoff: CDS >= ", 100 * opt.df$min.cds[best.idx], "%\n", sep=""))

  # Return data for plotting
  return(opt.df)
}


# Plot optimization data
plot.opt <- function(opt.df){
  # Get plot parameters
  min.val <- min(c(opt.df$power.frac, opt.df$or.frac), na.rm=T)
  best.cds <- opt.df$min.cds[which(opt.df$pareto.dist == min(opt.df$pareto.dist, na.rm=T))]

  # Prep plot area
  par(bty="n", mar=c(2.5, 3, 1, 0.2))
  plot(NA, xlim=c(0, 1), ylim=c(min.val, 1),
       xaxt="n", yaxt="n", xlab="", ylab="", yaxs="i")
  rect(xleft=best.cds-0.005, xright=best.cds+0.005,
       ybottom=par("usr")[3], ytop=par("usr")[4],
       col=adjustcolor(highlight.color, alpha=0.3),
       border=highlight.color, xpd=T)
  points(x=best.cds, y=par("usr")[4]+(0.02*diff(par("usr")[3:4])), xpd=T, pch=25, bg="black")

  # Plot power & effect size
  points(opt.df$min.cds, opt.df$or.frac, lwd=3, col="#F95700FF", type="l")
  points(opt.df$min.cds, opt.df$power.frac, lwd=3, col="#00A4CCFF", type="l")

  # Add axes
  axis(1, labels=NA, tck=-0.02)
  axis(1, tick=F, line=-0.65, cex.axis=0.85)
  mtext(1, text="Minimum CDS Cutoff", line=1.25)
  axis(2, labels=NA, tck=-0.02)
  axis(2, tick=F, line=-0.65, cex.axis=0.85, las=2)
  mtext(2, text="Fraction of Maximum Possible", line=2)

  # Add legend
  legend("bottomright", cex=0.85, lwd=3, col=c("#F95700FF", "#00A4CCFF"),
         legend=c("Effect Size (lnOR)", "Power (-log10[P])"), bty="n")
}


# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog meta.stats out.prefix",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Writes args & opts to vars
meta.in <- args$args[1]
out.prefix <- args$args[2]

# Load meta stats
stats <- load.meta.stats(meta.in, keep.n.cols=1)

# Compute optimization data
opt.df <- optimize.cds(stats)

# Plot optimization data
pdf(paste(out.prefix, "cds_optimization_results.pdf", sep="."),
    height=3, width=4)
plot.opt(opt.df)
dev.off()

# Write optimization data to .tsv
colnames(opt.df)[1] <- paste("#", colnames(opt.df)[1], sep="")
write.table(opt.df, paste(out.prefix, "cds_optimization_results.tsv", sep="."),
            col.names=T, row.names=F, sep="\t", quote=F)
