#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot summary distributions of outcome of gene-based fine-mapping for rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load table of PIPs
load.pips <- function(path){
  pips <- read.table(path, header=T, sep="\t", comment.char="")
  colnames(pips)[1] <- gsub("X.", "", colnames(pips)[1], fixed=T)

  # Joint model is already averaged across phenotypes,
  # but need to do the same for posterior and prior
  for(gene in unique(pips$gene)){
    gene.idxs <- which(pips$gene == gene)
    pips$PIP[gene.idxs] <- mean(pips$PIP[gene.idxs], na.rm=T)
  }

  # Drop HPO and deduplicate
  pips <- pips[, c("gene", "PIP", "cnv")]

  return(pips[which(!duplicated(pips)), ])
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot histogram of number of genes per credible set
plot.cs.size.hist <- function(credsets, sig.only=TRUE,  bin.width=1,
                              parmar=c(2.15, 2.85, 1, 0.6)){
  # Get plotting data
  if(sig.only){
    col <- "n_sig_in_credset"
    x.title <- "Sig. Genes per Credible Set"
  }else{
    col <- "n_genes"
    x.title <- "Genes per Credible Set"
  }
  pdat <- lapply(c("DEL", "DUP"), function(cnv){
    credsets[which(credsets$cnv == cnv), col]
  })
  xmax <- max(unlist(pdat))
  breaks <- seq(0, xmax+1, bin.width)
  counts.del <- hist(pdat[[1]], breaks=breaks, plot=F)$counts
  counts.dup <- hist(pdat[[2]], breaks=breaks, plot=F)$counts
  counts <- counts.del + counts.dup
  bar.colors <- c(rep(cnv.colors[1], length(counts)),
                  rep(cnv.colors[2], length(counts)))

  # Prepare plot area
  par(bty="n", mar=parmar)
  plot(NA, xlim=c(0, xmax), ylim=c(0, 1.025*max(counts)),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  y.ax.at <- axTicks(2)

  # Add stacked histogram
  rect(xleft=breaks[-length(breaks)], xright=breaks[-1],
       ybottom=c(rep(0, length(counts.del)), counts.del),
       ytop=c(counts.del, counts.del+counts.dup),
       col=bar.colors, border=bar.colors, lwd=0.5, xpd=T)
  rect(xleft=breaks[-length(breaks)], xright=breaks[-1],
       ybottom=rep(0, length(counts.del)),
       ytop=counts.del+counts.dup,
       col=NA, border=blueblack, xpd=T)

  # Add axes
  x.ax.at <- axTicks(1)
  axis(1, x.ax.at, tck=-0.025, labels=NA, col=blueblack)
  sapply(x.ax.at, function(x){axis(1, x, tick=F, line=-0.7)})
  mtext(1, line=1.25, text=x.title)
  axis(2, c(0, 10e10), tck=0, labels=NA, col=blueblack)
  axis(2, y.ax.at, tck=-0.025, labels=NA, col=blueblack)
  axis(2, y.ax.at, tick=F, las=2, line=-0.6)
  mtext(2, line=1.9, text="Credible Sets")
}


# Scatterplot of two sets of PIPs matched on gene, CNV, and phenotype
pip.scatter <- function(stats1, stats2, label1, label2,
                        pt.cex=0.75, pt.lwd=1, blue.bg=TRUE,
                        parmar=c(2.75, 2.75, 0.25, 0.25)){
  # Join PIP sets & assign colors
  x <- merge(stats1, stats2, by=c("gene", "cnv"), suffix=c(".1", ".2"))
  set.seed(2021)
  x <- x[sample(1:nrow(x), nrow(x), replace=F), ]
  pt.pch <- apply(x[, c("PIP.1", "PIP.2")], 1, function(vals){
    if(as.numeric(vals[1]) >= as.numeric(vals)[2]){25}else{24}
  })
  pt.col <- cnv.colors[x$cnv]
  pt.bcol <- cnv.blacks[x$cnv]
  x <- x[, grep("PIP", colnames(x), fixed=T)]
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
  plot(NA, xlim=c(0, 1), ylim=c(0, 1), xaxt="n", yaxt="n", xlab="", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=par("usr")[4],
       col=plot.bg, border=plot.border, bty=plot.bty)
  ax.at <- seq(0, 1, 0.2)
  abline(h=ax.at, v=ax.at, col=grid.col)
  abline(0, 1, lty=5, col=blueblack)

  # Add axes
  axis(1, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(1, at=ax.at, labels=NA, tck=-0.025, col=blueblack)
  sapply(ax.at, function(x){axis(1, at=x, tick=F, line=-0.75)})
  mtext(1, text=label1, line=1.25)
  axis(2, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(2, at=ax.at, labels=NA, tck=-0.025, col=blueblack)
  sapply(ax.at, function(y){axis(2, at=y, tick=F, line=-0.65, las=2)})
  mtext(2, text=label2, line=1.7)

  # Add points
  points(x, pch=pt.pch, cex=pt.cex, col=pt.bcol, bg=pt.col, lwd=pt.lwd)
}


# Histogram of numeric difference in PIPs matched on gene, CNV, and phenotype
pip.split.hist <- function(stats1, stats2, label1, label2, bin.width=0.05,
                           parmar=c(2.1, 2.75, 0.25, 0.25)){
  # Join PIP sets
  x <- merge(stats1, stats2, by=c("gene", "cnv"), suffix=c(".1", ".2"))

  # Compute PIP differences split by CNV type
  d <- lapply(c("DEL", "DUP"), function(cnv){
    x <- x[which(x$cnv==cnv), grep("PIP", colnames(x), fixed=T)]
    x$PIP.2 - x$PIP.1
  })

  # Compute bar counts
  breaks <- seq(-1, 1, by=bin.width)
  counts.del <- list(hist(d[[1]], breaks=breaks, plot=F)$counts[which(breaks<0)],
                     hist(d[[1]], breaks=breaks, plot=F)$counts[setdiff(which(breaks>=0), length(breaks))])
  counts.dup <- list(hist(d[[2]], breaks=breaks, plot=F)$counts[which(breaks<0)],
                     hist(d[[2]], breaks=breaks, plot=F)$counts[setdiff(which(breaks>=0), length(breaks))])
  counts <- lapply(1:2, function(i){counts.del[[i]] + counts.dup[[i]]})

  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(-1, 1), ylim=c(0, 1.03*max(unlist(counts))),
       xaxt="n", yaxt="n", xlab="", ylab="", yaxs="i")
  x.ax.at <- seq(-1, 1, 0.5)
  y.ax.at <- axTicks(2)

  # Add bars
  bar.colors <- lapply(1:2, function(i){
    c(rep(control.cnv.colors[i], length(list(counts.del, counts.dup)[[i]][[1]])-1),
      rep(cnv.whites[i], 2),
      rep(cnv.colors[i], length(list(counts.del, counts.dup)[[i]][[1]])-1))
  })
  rect(xleft=breaks[-length(breaks)], xright=breaks[-1],
       ybottom=0, ytop=unlist(counts.del),
       col=bar.colors[[1]], border=bar.colors[[1]], lwd=0.5)
  rect(xleft=breaks[-length(breaks)], xright=breaks[-1],
       ybottom=unlist(counts.del), ytop=unlist(counts),
       col=bar.colors[[2]], border=bar.colors[[2]], lwd=0.5)
  rect(xleft=breaks[-length(breaks)], xright=breaks[-1],
       ybottom=0, ytop=unlist(counts), col=NA,
       border=c(rep(cnv.colors[2], length(counts[[1]])-1),
                NA, NA,
                rep(cnv.blacks[2], length(counts[[2]])-1)))

  # Add axes
  axis(1, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(1, at=x.ax.at, labels=NA, tck=-0.03, col=blueblack)
  sapply(x.ax.at, function(x){axis(1, at=x, tick=F, line=-0.75)})
  mtext(1, text=bquote(.(label2) - .(label1)), line=1.2)
  axis(2, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(2, at=y.ax.at, labels=NA, tck=-0.03, col=blueblack)
  sapply(y.ax.at, function(y){
    axis(2, at=y, tick=F, line=-0.65, las=2, cex.axis=5/6, labels=prettyNum(y, big.mark=","))})
  mtext(2, text="Genes", line=1.7)

  # Add annotations
  abline(v=c(-1, 1)*bin.width, lty=5, col=blueblack)
  n.left <- sum(counts[[1]][-length(counts[[1]])])
  text(x=par("usr")[1]-(2*bin.width), y=0.6*par("usr")[4], pos=4,
       labels=paste(label1, "\ngreater by\n", bin.width, " for\n",
                    prettyNum(n.left, big.mark=","), " genes", sep=""),
       cex=5/6, font=3, col=control.cnv.colors[2])
  n.right <- sum(counts[[2]][-1])
  text(x=par("usr")[2]+(2*bin.width), y=0.6*par("usr")[4], pos=2,
       labels=paste(label2, "\ngreater by\n", bin.width, " for\n",
                    prettyNum(n.right, big.mark=","), " genes", sep=""),
       cex=5/6, font=3, col=cnv.colors[2])
}


# Plot histogram of PIPs for all genes in all credible sets
plot.pip.hist <- function(allgenes, conf.cutoff=0.2, vconf.cutoff=0.8,
                          breaks=seq(0, 1, 0.025), parmar=c(2.15, 2.85, 1, 0.6)){
  # Get plotting data
  pips <- lapply(c("DEL", "DUP"), function(cnv){
    all.hits <- allgenes[which(allgenes$cnv==cnv & !is.na(allgenes$credible_set)),
                         c("gene", "PIP")]
    all.hits$PIP[!duplicated(all.hits)]
  })
  counts.del <- hist(pips[[1]], breaks=breaks, plot=F)$counts
  counts.dup <- hist(pips[[2]], breaks=breaks, plot=F)$counts
  counts <- counts.del + counts.dup
  bar.colors <- lapply(c("DEL", "DUP"), function(cnv){
    sapply(breaks[-1], function(x){
      x <- round(x, 2)
      if(x <= conf.cutoff){
        cnv.whites[cnv]
      }else if(x <= vconf.cutoff){
        control.cnv.colors[cnv]
      }else{
        cnv.colors[cnv]
      }
    })
  })

  # Prepare plot area
  par(bty="n", mar=parmar)
  plot(NA, xlim=c(0, 1), ylim=c(0, 1.025*max(counts)),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  y.ax.at <- axTicks(2)
  rect(xleft=c(conf.cutoff, vconf.cutoff), xright=rep(par("usr")[2], 2),
       ybottom=par("usr")[3], ytop=par("usr")[4], border=NA, bty="n",
       col=adjustcolor(bluewhite, alpha=0.5))
  abline(v=c(conf.cutoff, vconf.cutoff), lty=2, col=blueblack)
  text(x=c(conf.cutoff, conf.cutoff, vconf.cutoff)+c(-0.05, 0.05, 0.05),
       y=sum(par("usr")[3:4])/2, srt=90, col=c(ns.color, rep(blueblack, 2)),
       labels=c("Unlikely", "Confident", "Highly Confident"), cex=5/6)

  # Add stacked histogram
  rect(xleft=breaks[-length(breaks)], xright=breaks[-1],
       ybottom=c(rep(0, length(counts.del)), counts.del),
       ytop=c(counts.del, counts.del+counts.dup),
       col=unlist(bar.colors), border=unlist(bar.colors), lwd=0.5, xpd=T)
  rect(xleft=breaks[-length(breaks)], xright=breaks[-1],
       ybottom=rep(0, length(counts.del)),
       ytop=counts.del+counts.dup,
       col=NA, border=blueblack, xpd=T)

  # Add axes
  x.ax.at <- axTicks(1)
  axis(1, x.ax.at, tck=-0.025, labels=NA, col=blueblack)
  sapply(x.ax.at, function(x){axis(1, x, tick=F, line=-0.7)})
  mtext(1, line=1.25, text="Posterior Inclusion Prob. (PIP)")
  axis(2, c(0, 10e10), tck=0, labels=NA, col=blueblack)
  axis(2, y.ax.at, tck=-0.025, labels=NA, col=blueblack)
  axis(2, y.ax.at, tick=F, las=2, line=-0.6)
  mtext(2, line=1.9, text="Genes in Credible Sets")
}


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(optparse, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog credsets.bed assocs.bed prior.pips.tsv",
                                            "posterior.pips.tsv fullmodel.pips.tsv",
                                            "out.prefix"),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 6){
  stop(paste("Four positional arguments required: credsets.bed, assocs.bed,",
             "prior.pips.tsv, posterior.pips.tsv, fullmodel.pips.tsv, and out.prefix\n"))
}

# Writes args & opts to vars
credsets.in <- args$args[1]
assocs.in <- args$args[2]
prior.pips.in <- args$args[3]
posterior.pips.in <- args$args[4]
final.pips.in <- args$args[5]
out.prefix <- args$args[6]

# # DEV PARAMETERS
# credsets.in <- "~/scratch/rCNV.final_genes.credible_sets.bed.gz"
# assocs.in <- "~/scratch/rCNV.final_genes.associations.bed.gz"
# prior.pips.in <- "~/scratch/all_PIPs.prior.tsv"
# posterior.pips.in <- "~/scratch/all_PIPs.posterior.tsv"
# final.pips.in <- "~/scratch/all_PIPs.full_model.tsv"
# out.prefix <- "~/scratch/finemap_distribs_test"

# Load credible sets and associations
credsets <- load.credsets(credsets.in)
assocs <- load.gene.associations(assocs.in)

# Load PIPs for all genes
pips <- list("Prior" = load.pips(prior.pips.in),
             "Posterior" = load.pips(posterior.pips.in),
             "Full Model" = load.pips(final.pips.in),
             "full_verbose" = read.table(final.pips.in, header=T, sep="\t",
                                         comment.char="", check.names=F))
colnames(pips$full_verbose)[1] <- "HPO"

# Histogram of number of significant genes per credible set
sapply(c(TRUE, FALSE), function(sig.only){
  if(sig.only){
    suffix <- "sig_genes"
  }else{
    suffix <- "all_genes"
  }
  pdf(paste(out.prefix, "credset_size_hist", suffix, ".pdf", sep="."),
      height=2.1, width=2.5)
  plot.cs.size.hist(credsets, sig.only, parmar=c(2.25, 2.85, 0.3, 1))
  dev.off()
})

# Compute average decrease in # of significant genes per credset due to fine-mapping
finemap.decrease.pct <- 100*(mean(credsets$n_sig_genes) - mean(credsets$n_sig_in_credset)) / mean(credsets$n_sig_genes)
cat(paste("Fine-mapping reduced the average number of significant genes per block by ",
          round(finemap.decrease.pct), "%\n", sep=""))

# Scatterplot of effect of fine-mapping on number of genes per association
n_gene_range <- c(0, max(credsets[, c("n_sig_genes", "n_sig_in_credset")]))
pdf(paste(out.prefix, "genes_per_credset_before_vs_after.scatter.pdf", sep="."),
    height=2.5, width=2.6)
credsets.scatter(credsets, credsets$n_sig_genes, credsets$n_sig_in_credset, add.lm=T,
                 xlims=n_gene_range, ylims=n_gene_range,
                 abline.a=0, abline.b=1, abline.lty=5, blue.bg=FALSE, pt.cex=0.7,
                 xtitle="Genes Before Fine-Mapping", x.title.line=1.25,
                 ytitle="Genes After Fine-Mapping",
                 parmar = c(2.4, 2.7, 0.25, 0.5))
dev.off()

# Plot pairwise comparisons of pips
sapply(list(c(1, 2), c(1, 3), c(2, 3)), function(idxs){
  set1 <- pips[[idxs[1]]]
  name1 <- names(pips)[idxs[1]]
  set2 <- pips[[idxs[2]]]
  name2 <- names(pips)[idxs[2]]
  # Scatterplot
  pdf(paste(out.prefix, "pip_scatter",
            gsub(" ", "_", name1, fixed=T),
            gsub(" ", "_", name2, fixed=T),
            "pdf", sep="."),
      height=2.4, width=2.4)
  pip.scatter(set1, set2, name1, name2, blue.bg=FALSE, pt.cex=0.4)
  dev.off()
  # Histogram of residuals
  pdf(paste(out.prefix, "residual_hist",
            gsub(" ", "_", name1, fixed=T),
            gsub(" ", "_", name2, fixed=T),
            "pdf", sep="."),
      height=1.8, width=2.4)
  pip.split.hist(set1, set2, name1, name2)
  dev.off()
})

# Plot histograms of PIPs for all genes in credible sets
pdf(paste(out.prefix, "finemapped_distribs.credset_pip.pdf", sep="."),
    height=2.1, width=2.7)
plot.pip.hist(pips[[4]])
dev.off()
