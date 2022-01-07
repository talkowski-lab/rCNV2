#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot geneset enrichments of fine-mapped genes for rCNV2 manuscript


######################
### DATA FUNCTIONS ###
######################
# Load a list of gene lists
load.gene.lists <- function(genelists.in){
  x <- read.table(genelists.in, header=F, sep="\t")
  genelists <- lapply(x[, 2], function(path){unique(read.table(path, header=F)[, 1])})
  names(genelists) <- as.character(x[, 1])
  return(genelists)
}

# Append one-hot features of membership per gene for a list of gene lists
append.glist.feats <- function(feats, genelists){
  for(i in 1:length(genelists)){
    feats[, names(genelists)[i]] <- as.numeric(feats$gene %in% genelists[[i]])
  }
  return(feats)
}

# Calculate excess PTVs for DDD
calc.ddd.excess <- function(feats){
  # Scale LoF mutation rate by total number of PTVs in DDD
  lof.exp <- feats$gnomad_mu_lof * (sum(feats$ddd_dn_lof)/sum(feats$gnomad_mu_lof))
  feats$ddd_dn_lof - lof.exp
}

# Compute enrichment for a single feature for a pre-specified list of genes
calc.enrichment <- function(genes, feats, feat, ci="bootstrap", add.sig=FALSE,
                            ref.genes=NULL){
  # Get values for genes in question
  vals <- feats[which(feats$gene %in% genes), feat]

  # Compute statistics depending on confidence interval specified
  if(ci=="bootstrap"){
    # Bootstrap feature mean & 95% CI
    mean.for.boot <- function(vals, indexes){mean(vals[indexes], na.rm=T)}
    estimate <- mean(vals, na.rm=T)
    ci <- boot.ci(boot(data=vals, statistic=mean.for.boot, R=1000),
                  conf=0.95, type="norm")$normal[, -1]
    res <- as.numeric(c(estimate, ci))
  }else if(ci=="binomial"){
    # Compute binomial mean & 95% CI
    baseline <- mean(feats[, feat], na.rm=T)
    binom.res <- binom.test(x=sum(vals), n=length(vals), p=baseline)
    res <- as.numeric(c(binom.res$estimate, binom.res$conf.int))
  }

  # Add Fisher's exact significance test results, if optioned
  if(add.sig==TRUE){
    if(is.null(ref.genes)){
      ref.vals <- feats[which(!(feats$gene %in% genes)), feat]
    }else{
      ref.vals <- feats[which(feats$gene %in% ref.genes), feat]
    }
    fisher.res <- fisher.test(matrix(c(length(which(ref.vals==0)), length(which(vals==0)),
             length(which(ref.vals>0)), length(which(vals>0))),
           nrow=2, byrow=T))
    res <- c(res, as.numeric(c(fisher.res$estimate, fisher.res$p.value)))
  }

  return(res)
}

# Compute enrichment for a single feature across all groups of genes for plotting
get.plot.data <- function(feats, feat, gene.groups, not.credset.genes, ci="bootstrap", add.sig=FALSE){
  # Group 1: PIP ≥ 0.2, top gene
  # Group 2: PIP ≥ 0.2, not top gene
  # Group 3: PIP < 0.2, not top gene
  # Group 4: all genes outside of credible sets
  lapply(c("DEL", "DUP", "CNV"), function(cnv){
    g1 <- unique(c(gene.groups[[cnv]]$top.vconf, gene.groups[[cnv]]$top.conf))
    g2 <- unique(c(gene.groups[[cnv]]$nottop.vconf, gene.groups[[cnv]]$nottop.conf))
    g3 <- unique(gene.groups[[cnv]]$nottop.notconf)
    g4 <- not.credset.genes
    pdat <- as.data.frame(do.call("rbind", lapply(list(g1, g2, g3, g4), calc.enrichment,
                                                  feats=feats, feat=feat, ci=ci,
                                                  add.sig=add.sig, ref.genes=g4)))
    if(add.sig==FALSE){
      colnames(pdat) <- c("mean", "lower", "upper")
    }else{
      colnames(pdat) <- c("mean", "lower", "upper", "odds.ratio", "p.value")
    }
    return(pdat)
  })
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot enrichment of gene groups for a single feature
plot.enrichment <- function(feats, feat, gene.groups, ci="bootstrap",
                            print.fisher=FALSE, bars=FALSE, pt.cex=1,
                            y.ax.title="Proportion of Genes", y.ax.title.line=0.75,
                            add.y.labels=TRUE, blue.bg=TRUE,
                            parmar=c(0.2, 4, 0.2, 0.2)){
  # Collect data for plotting
  not.credset.genes <- feats$gene[which(!(feats$gene %in% all.credset.genes))]
  pdat <- get.plot.data(feats, feat, gene.groups, not.credset.genes, ci, add.sig=print.fisher)

  # Print Fisher exact test results, if optioned
  if(print.fisher==TRUE){
    or <- round(pdat[[3]]$odds.ratio[1], 2)
    pval <- format(pdat[[3]]$p.value[1], scientific=TRUE, nsmall=2, digits=4)
    cat(paste("Two-sided Fisher's exact test of top-ranked confident genes\nvs.",
              "genes not present in any credible set:",
              "\nOdds ratio =", or, "\nP-value =", pval, "\n\n"))
  }

  # Get plot parameters
  ymax <- min(c(max(unlist(pdat), na.rm=T),
                1.3*max(do.call("rbind", pdat)[, 1], na.rm=T)))
  min.means <- min(do.call("rbind", pdat)[, 1], na.rm=T)
  if(bars==TRUE){
    ymin <- 0
  }else{
    ymin <- max(c(min(unlist(pdat), na.rm=T),
                  min.means - 0.5*abs(min.means)))
  }
  if(ci=="binomial"){
    ymax <- min(c(ymax, 1))
    ymin <- max(c(0, ymin))
  }
  cell.hex <- 0.15
  cell.height <- cell.hex*diff(c(ymin, ymax))
  cell.buffer <- 0.3*cell.height
  ylims <- c(ymin-(2*cell.height)-cell.buffer, ymax)
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
  plot(NA, xlim=c(0, 4), ylim=ylims,
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")

  # Set plot parameters
  x.at <- (1:4)-0.5
  x.at.cnv <- list("DEL"=x.at-0.1, "DUP"=x.at+0.1)
  y.ax.at <- sort(unique(round(10*axTicks(2)[which(axTicks(2) > ymin)])))/10

  # Add background shading
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       border=plot.border, bty=plot.bty, col=plot.bg)
  abline(h=y.ax.at, col=grid.col)
  abline(h=pdat[[1]]$mean[4], lty=5, col=ns.color)

  # Add points (or bars, if bars==TRUE)
  point.colors <- lapply(1:2, function(i){c(cnv.colors[i], rep(control.cnv.colors[i], 2))})
  if(bars==TRUE){
    bar.border <- lapply(1:2, function(i){c(cnv.blacks[i], rep(cnv.colors[i], 2))})
    ci.colors <- list(c(redblack, rep(cnv.colors[1], 2)),
                      c(blueblack, rep(cnv.colors[2], 2)))
    ci.lwd <- 1
    ci.x.at <- list(x.at-0.2, x.at+0.2)
  }else{
    ci.colors <- point.colors
    ci.lwd <- 1.5
    ci.x.at <- x.at.cnv
  }
  if(bars==TRUE){
    lapply(1:2, function(i){
      rect(xleft=x.at[1:3]+(0.4*(i-2)), xright=x.at[1:3]+(0.4*(i-1)),
           ybottom=rep(0, 3), ytop=pdat[[i]]$mean[1:3],
           col=point.colors[[i]], border=bar.border[[i]])
      segments(x0=rep(ci.x.at[[i]][1:3]-0.1, 2), x1=rep(ci.x.at[[i]][1:3]+0.1, 2),
               y0=c(pdat[[i]]$lower[1:3], pdat[[i]]$upper[1:3]),
               y1=c(pdat[[i]]$lower[1:3], pdat[[i]]$upper[1:3]),
               lwd=ci.lwd, col=ci.colors[[i]], lend="butt")
    })
    rect(xleft=x.at[4]-0.2, xright=x.at[4]+0.2,
         ybottom=0, ytop=pdat[[1]]$mean[4], col=ns.color, border=ns.color.dark)
    segments(x0=x.at[4], x1=x.at[4], y0=pdat[[1]]$lower[4], y1=pdat[[1]]$upper[4],
             lwd=ci.lwd, col=ns.color.dark, lend="round")
    segments(x0=rep(x.at[4]-0.1, 2), x1=rep(x.at[4]+0.1, 2),
             y0=c(pdat[[1]]$lower[4], pdat[[1]]$upper[4]),
             y1=c(pdat[[1]]$lower[4], pdat[[1]]$upper[4]),
             lwd=ci.lwd, col=ns.color.dark, lend="butt")
  }
  lapply(1:2, function(i){
    segments(x0=ci.x.at[[i]][1:3], x1=ci.x.at[[i]][1:3],
             y0=pdat[[i]]$lower[1:3], y1=pdat[[i]]$upper[1:3],
             lwd=ci.lwd, col=ci.colors[[i]], lend="round")
  })
  if(bars==F){
    segments(x0=x.at[4], x1=x.at[4], y0=pdat[[1]]$lower[4], y1=pdat[[1]]$upper[4],
             lwd=ci.lwd, col=ns.color, lend="round")
    lapply(1:2, function(i){
      points(x=x.at.cnv[[i]][1:3], y=pdat[[i]]$mean[1:3],
             col=point.colors[[i]], pch=19, cex=pt.cex)
    })
    points(x=x.at[4], y=pdat[[1]]$mean[4], col=ns.color, pch=19, cex=pt.cex)
  }

  # Add lower grid
  rect(xleft=par("usr")[1], xright=par("usr")[2], ybottom=par("usr")[3], ytop=ymin,
       border=NA, bty="n", xpd=T, col="white")
  rect(xleft=par("usr")[1], xright=1, ybottom=par("usr")[3], ytop=ymin-cell.buffer,
       border=NA, bty="n", xpd=T, col=adjustcolor(highlight.color, alpha=0.3))
  sapply(1:2, function(i){
    rect(xleft=0:2, xright=1:3,
         ybottom=ymin-((i-1)*cell.height)-cell.buffer,
         ytop=ymin-(i*cell.height)-cell.buffer,
         col=NA, border=blueblack, xpd=T)
  })
  text(x=x.at[1:3], y=ymin-(1.5*cell.height)-cell.buffer, cex=5/6,
       labels=c("Yes", "No", "No"), font=c(2, 1, 1))
  text(x=x.at[1:2], y=ymin-(0.5*cell.height)-cell.buffer, cex=5/6,
       labels=bquote("" >= 0.2))
  text(x=x.at[3], y=ymin-(0.5*cell.height)-cell.buffer, cex=5/6,
       labels=bquote("" < 0.2))
  text(x=x.at[4], y=ymin-cell.height-cell.buffer, cex=5/6, font=3,
       labels="All\nOther\nGenes", col=ns.color, xpd=T)
  if(add.y.labels==TRUE){
    axis(2, at=ymin-(c(0.5, 1.5)*cell.height)-cell.buffer, tick=F, line=-0.85, cex=5/6,
         labels=c("PIP", "Top"), las=2)
  }

  # Add remaining axes
  axis(2, at=c(ymin, 10e10), tck=0, labels=NA, col=blueblack)
  axis(2, at=y.ax.at, labels=NA, tck=-0.025, col=blueblack)
  axis(2, at=y.ax.at, tick=F, line=-0.65, las=2)
  if(add.y.labels==TRUE){
    axis(2, at=mean(c(ymin, ymax)), tick=F, line=y.ax.title.line, labels=y.ax.title)
  }
  if(bars==TRUE){
    abline(h=0, col=blueblack)
  }
}


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(optparse, quietly=T)
require(boot, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog credsets.bed assocs.bed",
                                            "features.bed genelists.tsv",
                                            "out.prefix"),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 5){
  stop(paste("Five positional arguments required: credsets.bed, assocs.bed,",
             "features.bed, genelists.tsv, and out.prefix\n"))
}

# Writes args & opts to vars
credsets.in <- args$args[1]
assocs.in <- args$args[2]
features.in <- args$args[3]
genelists.in <- args$args[4]
out.prefix <- args$args[5]

# # DEV PARAMETERS
# credsets.in <- "~/scratch/rCNV.final_genes.credible_sets.bed.gz"
# assocs.in <- "~/scratch/rCNV.final_genes.associations.bed.gz"
# features.in <- "~/scratch/gencode.v19.canonical.pext_filtered.all_features.bed.gz"
# genelists.in <- "~/scratch/gset_genelists.tmp.tsv"
# out.prefix <- "~/scratch/test_finemapped_enrichments"

# Load credible sets and associations
credsets <- load.credsets(credsets.in)
assocs <- load.gene.associations(assocs.in)

# Split fine-mapped genes into categories based on top/not top status and conf/vconf
gene.groups <- categorize.genes(credsets)
summarize.finemapping.categories(gene.groups)
all.credset.genes <- as.character(unlist(unlist(gene.groups)))

# Load gene-level features & extract list of all genes _not_ in any credible set
feats <- load.features(features.in)
genelists <- load.gene.lists(genelists.in)
feats <- append.glist.feats(feats, genelists)
feats$excess_ddd_dn_lof <- calc.ddd.excess(feats)

# Plot gene list-based panels
sapply(names(genelists), function(glist){
  pdf(paste(out.prefix, glist, "enrichments.nolabel.pdf", sep="."),
      height=2.1, width=2.6)
  cat(paste(glist, "\n"))
  plot.enrichment(feats, glist, gene.groups, ci="binomial", add.y.labels=FALSE,
                  bars=TRUE, print.fisher=TRUE, blue.bg=FALSE)
  dev.off()
  pdf(paste(out.prefix, glist, "enrichments.withlabel.pdf", sep="."),
      height=2.1, width=2.6)
  plot.enrichment(feats, glist, gene.groups, ci="binomial", add.y.labels=TRUE,
                  bars=TRUE, blue.bg=FALSE)
  dev.off()
})

# Plot DNM counts from DDD
pdf(paste(out.prefix, "ddd_dn_lof", "enrichments.pdf", sep="."),
    height=2.1, width=2.6)
plot.enrichment(feats, "excess_ddd_dn_lof", gene.groups, ci="bootstrap", pt.cex=1.25,
                y.ax.title="Excess De Novo\nPTVs per Gene", y.ax.title.line=0.9,
                blue.bg=FALSE)
dev.off()
