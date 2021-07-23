#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot quality assessment analyses of various gene fine-mapping approaches


options(scipen=100000, stringsAsFactors=F)
require(rCNV2)


#################
### FUNCTIONS ###
#################
# Load truth set
load.truth <- function(truth.in){
  truth.list <- read.table(truth.in, header=F, sep="\t")
  truth.sets <- lapply(1:nrow(truth.list), function(i){
    x <- read.table(truth.list[i, 2], sep="\t", header=T, comment.char="")
    colnames(x)[1] <- "HPO"
    return(list("truth.genes"=x, "name"=truth.list[i, 1]))
  })
  names(truth.sets) <- truth.list[, 1]
  return(truth.sets)
}

# Compute ROC
roc <- function(stats, steps=seq(1, 0, -0.001)){
  roc_res <- as.data.frame(t(sapply(steps, function(minPIP){
    idx.all <- which(stats$PIP > minPIP)
    fall <- length(idx.all) / nrow(stats)
    # Note: performance statistics only defined for known condition positives or negatives
    stats <- stats[which(stats$true | stats$false), ]
    idxs <- which(stats$PIP > minPIP)
    ftrue <- length(which(stats$true[idxs])) / length(which(stats$true))
    ffalse <- length(which(stats$false[idxs])) / length(which(stats$false))
    return(c(minPIP, fall, ffalse, ftrue))
  })))
  colnames(roc_res) <- c("minPIP", "frac_all", "frac_false", "frac_true")
  return(roc_res)
}

# Compute PRC
prc <- function(stats, steps=seq(1, 0, -0.001)){
  prc_res <- as.data.frame(t(sapply(steps, function(minPIP){
    idx.all <- which(stats$PIP > minPIP)
    fall <- length(idx.all) / nrow(stats)
    # Note: performance statistics only defined for known condition positives or negatives
    stats <- stats[which(stats$true | stats$false), ]
    idxs <- which(stats$PIP > minPIP)
    prec <- length(which(stats$true[idxs])) / length(idxs)
    recall <- length(which(stats$true[idxs])) / length(which(stats$true))
    return(c(minPIP, fall, prec, recall))
  })))
  colnames(prc_res) <- c("minPIP", "frac_all", "precision", "recall")
  return(prc_res)
}

# Load a single dataset and annotate vs truth sets
load.data.single <- function(path, truth, false.genes){
  x <- read.table(path, header=T, sep="\t", comment.char="")[, 1:4]
  colnames(x)[1] <- "HPO"
  truth$true <- TRUE
  x <- merge(x, truth, all.x=T, all.y=F, by=c("HPO", "gene"))
  x$true[which(is.na(x$true))] <- FALSE
  x$false <- FALSE
  x$false[which(x$gene %in% false.genes & x$true == FALSE)] <- TRUE
  return(x)
}

# Find ROC-optimal point
optimize.roc <- function(roc.res, return.dist=F){
  # Compute Euclidean distance from (0, 1)
  x2 <- (roc.res$frac_false - 0) ^ 2
  y2 <- (roc.res$frac_true - 1) ^ 2
  d <- sqrt(x2 + y2)
  d.best <- min(d, na.rm=T)
  best.idx <- sort(which(d == d.best))[1]
  if(return.dist==F){
    return(roc.res[best.idx, ])
  }else{
    return(d)
  }
}

# Find PRC-optimal point
optimize.prc <- function(prc.res, return.dist=F){
  # Compute Euclidean distance from (1, 1)
  x2 <- (prc.res$recall - 1) ^ 2
  y2 <- (prc.res$precision - 1) ^ 2
  d <- sqrt(x2 + y2)
  d.best <- min(d, na.rm=T)
  best.idx <- sort(which(d == d.best))[1]
  if(return.dist==F){
    return(prc.res[best.idx, ])
  }else{
    return(d)
  }
}

# Gather calibration
get.calibration <- function(stats, n.pip.bins=20){
  pip.breaks <- quantile(stats$PIP, seq(0, 1, by=1/n.pip.bins), na.rm=T)
  cal <- as.data.frame(t(sapply(1:(length(pip.breaks)-1), function(i){
    idx <- which(stats$PIP > pip.breaks[i] & stats$PIP <= pip.breaks[i+1])
    true <- length(which(stats$true[idx]))
    binconf(true, length(idx))
  })))
  colnames(cal) <- c("est", "lower", "upper")
  return(cal)
}

# Wrapper to load all datasets
load.datasets <- function(data.in, truth, false.genes, subgroup="all"){
  datlist <- read.table(data.in, header=F, sep="\t")
  colnames(datlist) <- c("name", "color", "lty", "path")
  data <- lapply(1:nrow(datlist), function(i){
    stats <- load.data.single(datlist$path[i], truth, false.genes)
    if(subgroup == "developmental"){
      stats <- stats[which(stats$HPO %in% developmental.hpos), ]
    }else if(subgroup == "adult"){
      stats <- stats[which(stats$HPO %in% adult.hpos), ]
    }
    roc.res <- roc(stats)
    roc.opt <- optimize.roc(roc.res)
    prc.res <- prc(stats)
    prc.opt <- optimize.prc(prc.res)
    if(any(stats$true)){
      roc.auc <- flux::auc(roc.res$frac_false, roc.res$frac_true)
      prc.auc <- flux::auc(prc.res$recall, prc.res$precision)
    }else{
      roc.auc <- NA
      prc.auc <- NA
    }
    return(list("stats"=stats,
                "roc"=roc.res,
                "roc.opt"=roc.opt,
                "roc.auc"=roc.auc,
                "roc.PIP_0.1"=roc.res[which(round(roc.res$minPIP, 3) == 0.100), ],
                "roc.PIP_0.3"=roc.res[which(round(roc.res$minPIP, 3) == 0.300), ],
                "roc.PIP_0.5"=roc.res[which(round(roc.res$minPIP, 3) == 0.500), ],
                "roc.PIP_0.7"=roc.res[which(round(roc.res$minPIP, 3) == 0.700), ],
                "roc.PIP_0.9"=roc.res[which(round(roc.res$minPIP, 3) == 0.900), ],
                "prc"=prc.res,
                "prc.opt"=prc.opt,
                "prc.auc"=prc.auc,
                "prc.PIP_0.1"=prc.res[which(round(prc.res$minPIP, 3) == 0.100), ],
                "prc.PIP_0.3"=prc.res[which(round(prc.res$minPIP, 3) == 0.300), ],
                "prc.PIP_0.5"=prc.res[which(round(prc.res$minPIP, 3) == 0.500), ],
                "prc.PIP_0.7"=prc.res[which(round(prc.res$minPIP, 3) == 0.700), ],
                "prc.PIP_0.9"=prc.res[which(round(prc.res$minPIP, 3) == 0.900), ],
                "calibration"=get.calibration(stats),
                "color"=datlist$color[i],
                "lty"=datlist$lty[i]))
  })
  names(data) <- datlist$name
  return(data)
}

# PIP comparison scatterplot
pip.scatter <- function(stats1, stats2, color1, color2, label1, label2,
                        title.cex=0.75, lpos="right"){
  # Join PIPs
  x <- merge(stats1, stats2, by=c("HPO", "gene"), suffix=c(".1", ".2"))
  x <- x[, grep("PIP", colnames(x), fixed=T)]
  colors <- apply(x, 1, function(p){if(p[1] > p[2]){color1}else{color2}})
  # Plot
  par(mar=c(4, 4, 2, 2))
  plot(x, xlim=c(0, 1), ylim=c(0, 1), col=colors, pch=19, cex=0.75,
       panel.first=c(abline(0, 1, lty=2)),
       xlab=paste("PIP (", label1, ")", sep=""),
       ylab=paste("PIP (", label2, ")", sep=""))
  mtext(3, text=paste(label1, "vs.", label2), font=2, line=0.2, cex=title.cex)
  legend(lpos, cex=0.7, col=c(color1, color2), pch=19,
         legend=c(paste(label1 , ">", label2),
                  paste(label2 , ">", label1)))
}

# PIP comparison histogram
pip.split.hist <- function(stats1, stats2, color1, color2, label1, label2,
                           title.cex=0.75){
  # Join PIPs
  x <- merge(stats1, stats2, by=c("HPO", "gene"), suffix=c(".1", ".2"))
  x <- x[, grep("PIP", colnames(x), fixed=T)]
  d <- x$PIP.2 - x$PIP.1
  # Plot
  par(mar=c(4, 4, 2, 2), bty="n")
  h <- hist(d, breaks=seq(-1, 1, 0.05), plot=F)
  plot(x=c(-1, 1), y=c(0, max(h$counts)), type="n",
       xlab=paste("PIP (", label2 , " - ", label1, ")", sep=""),
       ylab="Genes (count)")
  rect(xleft=h$breaks[which(h$mids<0)],
       xright=h$breaks[which(h$mids<0)+1],
       ybottom=0, ytop=h$counts[which(h$mids<0)], col=color1)
  rect(xleft=h$breaks[which(h$mids>0)],
       xright=h$breaks[which(h$mids>0)+1],
       ybottom=0, ytop=h$counts[which(h$mids>0)], col=color2)
  mtext(3, text=paste(label1, "vs.", label2), font=2, line=0.2, cex=title.cex)
  legend("topleft", cex=0.7, text.col=color1, bty="n",
         legend=paste(label1 , ">\n", label2))
  legend("topright", cex=0.7, text.col=color2, bty="n",
         legend=paste(label2 , ">\n", label1))
}

# ROC plot
plot.roc <- function(data, title=NULL){
  par(mar=c(3, 3, 1.5, 1))
  plot(x=c(0, 1), y=c(0, 1), type="n",
       xaxs="i", yaxs="i", xlab="", ylab="")
  abline(0, 1, col="gray70", lty=2)
  lorder <- order(-sapply(data, function(x){x$roc.auc}))
  sapply(rev(lorder), function(i){
    x <- data[[i]]
    points(x$roc$frac_false, x$roc$frac_true,
           type="l", col=x$color, lwd=2, lty=x$lty)
  })
  sapply(rev(lorder), function(i){
    x <- data[[i]]
    points(x$roc.opt$frac_false, x$roc.opt$frac_true,
           pch=21, bg=x$color)
    sapply(c(1, 3, 5, 7, 9), function(d){
      fname <- paste("roc.PIP_0", d, sep=".")
      vals <- unlist(x[which(names(x) == fname)])
      names(vals) <- colnames(x$roc)
      points(vals[3], vals[4],
             pch=23, bg="white", cex=1.2)
      text(x=vals[3], y=vals[4],
           labels=d, col=x$color, cex=0.65)
    })
  })
  legend("bottomright", lwd=5, cex=0.75, bty="n",
         col=sapply(data, function(x){x$color})[lorder],
         legend=paste(names(data), " (",
                      sapply(data, function(x){format(round(x$roc.auc, 2), nsmall=2)}),
                      ")", sep="")[lorder])
  legend("topleft", cex=0.75, bty="n", pch=c(21, 23),
         legend=c("ROC-optimal cutoff", "PIP > 0.N"),
         pt.bg=c("gray50", "white"), pt.cex=c(1, 1.4))
  legend("topleft", cex=0.75, bty="n", pch=c(NA, "N"),
         legend=c("", ""), pt.cex=0.5)
  mtext(1, line=2, text="Fraction of false genes retained")
  mtext(2, line=2, text="Fraction of known associations retained")
  mtext(3, line=0.1, text=title)
}


# PRC plot
plot.prc <- function(data, title=NULL){
  par(mar=c(3, 3, 1.5, 1))
  prec.max <- max(sapply(data, function(l){l$prc$precision}), na.rm=T)
  prec.baseline <- data[[1]]$prc$precision[which(data[[1]]$prc$minPIP==0)]
  plot(x=c(0, 1), y=c(0, prec.max), type="n",
       xaxs="i", yaxs="i", xlab="", ylab="")
  abline(h=prec.baseline, col="gray70", lty=2)
  lorder <- order(-sapply(data, function(x){x$prc.auc}))
  sapply(rev(lorder), function(i){
    x <- data[[i]]
    points(x$prc$recall, x$prc$precision,
           type="l", col=x$color, lwd=2, lty=x$lty)
  })
  sapply(rev(lorder), function(i){
    x <- data[[i]]
    points(x$prc.opt$recall, x$prc.opt$precision,
           pch=21, bg=x$color)
    sapply(c(1, 3, 5, 7, 9), function(d){
      fname <- paste("prc.PIP_0", d, sep=".")
      vals <- unlist(x[which(names(x) == fname)])
      names(vals) <- colnames(x$prc)
      points(vals[4], vals[3],
             pch=23, bg="white", cex=1.2)
      text(x=vals[4], y=vals[3],
           labels=d, col=x$color, cex=0.65)
    })
  })
  legend("topright", lwd=5, cex=0.75, bty="n",
         col=sapply(data, function(x){x$color})[lorder],
         legend=paste(names(data), " (",
                      sapply(data, function(x){format(round(x$prc.auc, 2), nsmall=2)}),
                      ")", sep="")[lorder])
  legend("bottomleft", cex=0.75, bty="n", pch=c(21, 23),
         legend=c("PRC-optimal cutoff", "PIP > 0.N"),
         pt.bg=c("gray50", "white"), pt.cex=c(1, 1.4))
  legend("bottomleft", cex=0.75, bty="n", pch=c(NA, "N"),
         legend=c("", ""), pt.cex=0.5)
  mtext(1, line=2, text="Recall (known associations)")
  mtext(2, line=2, text="Precision (known associations)")
  mtext(3, line=0.1, text=title)
}


# Calibration plot
plot.calibration <- function(data, title=NULL){
  plot.dat <- do.call("cbind", lapply(data, function(d){d$roc$frac_true/d$roc$frac_all}))
  xvals <- data[[1]]$roc$minPIP
  # keep.idx <- which(sapply(xvals, round, 3) %in% seq(0, 1, 0.05))
  # plot.dat <- plot.dat[keep.idx, ]
  # xvals <- xvals[keep.idx]
  par(mar=c(2.75, 4, 1.5, 1))
  plot(x=range(xvals), y=c(0, max(plot.dat, na.rm=T)), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="")
  abline(h=1, lty=2, col="gray80")
  sapply(1:ncol(plot.dat), function(i){
    # points(x=xvals, y=plot.dat[, i], col=data[[i]]$color, pch=19, cex=0.7)
    points(x=xvals, y=plot.dat[, i], col=data[[i]]$color, type="l", lwd=2)
  })
  axis(1, at=seq(0, 1, 0.2), labels=NA)
  axis(1, at=seq(0, 1, 0.2), tick=F, line=-0.5)
  mtext(1, line=1.75, text="Min PIP")
  axis(2, las=2)
  mtext(2, line=3, text="Fold-Enrichment")
  mtext(3, line=0.2, text=title, font=2)
  legend("bottomright", cex=0.6, pch=19, bg="white",
         col=sapply(data, function(x){x$color}),
         legend=names(data))
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse, quietly=T)
require(flux, quietly=T)
require(Hmisc, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog roc_input.tsv truth_genes.tsv false_genes.list subgroup out.prefix",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 5){
  stop("Five positional arguments: data.tsv, truth_sets.tsv, false_genes.list, subgroup, and out_prefix\n")
}

# Writes args & opts to vars
data.in <- args$args[1]
truth.in <- args$args[2]
false.in <- args$args[3]
pheno.group <- args$args[4]
out.prefix <- args$args[5]

# # DEV PARAMTERS
# setwd("~/scratch")
# data.in <- "finemap_roc_input.tsv"
# truth.in <- "finemap_roc_truth_sets.tsv"
# false.in <- "true_negatives.genes.tsv"
# pheno.group <- "developmental"
# out.prefix <- "finemap_roc"

# Read truth sets & negative genes
truth <- load.truth(truth.in)
false.genes <- sort(unique(as.character(read.table(false.in, header=F)[, 1])))

# Load each input, annotate vs. truth, and store as list
data <- lapply(truth, function(tlist){
  load.datasets(data.in, tlist$truth.genes, false.genes, pheno.group)
})
names(data) <- names(truth)

# Plot histograms of PIPs
pdf(paste(out.prefix, "PIP_distributions.pdf", sep="."),
    height=2.25, width=12)
par(mfrow=c(1, length(data[[1]])),
    mar=c(4, 4, 1.5, 1.5))
hmax <- max(sapply(1:length(data[[1]]), function(i){
  max(hist(data[[1]][[i]]$stats$PIP, breaks=seq(0, 1, 0.02), plot=F)$counts)
}))
sapply(1:length(data[[1]]), function(i){
  color <- data[[1]][[i]]$color
  hist(data[[1]][[i]]$stats$PIP, breaks=seq(0, 1, 0.02),
       col=color, border=color, xlab="PIP", ylab="Associations",
       main=names(data[[1]])[i], ylim=c(0, hmax))
})
dev.off()

# Plot PIP comparison scatterplots
if(names(data[[1]])[1] %in% c("Prior", "Posterior")){
  png(paste(out.prefix, "PIP_comparisons", tolower(names(data[[1]])[1]), "scatter.png", sep="."),
      height=2.5*300, width=2.5*(length(data[[1]])-1)*300, res=300)
  par(mfrow=c(1, length(data[[1]]) - 1))
  sapply(2:length(data[[1]]), function(i){
    pip.scatter(stats1=data[[1]][[1]]$stats,
                stats2=data[[1]][[i]]$stats,
                color1=data[[1]][[1]]$color,
                color2=data[[1]][[i]]$color,
                label1=names(data[[1]])[1],
                label2=names(data[[1]])[i],
                lpos="bottomright")
  })
  dev.off()
}
if(names(data[[1]])[2] %in% c("Prior", "Posterior")){
  png(paste(out.prefix, "PIP_comparisons", tolower(names(data[[1]])[2]), "scatter.png", sep="."),
      height=2.5*300, width=2.5*(length(data[[1]])-1)*300, res=300)
  par(mfrow=c(1, length(data[[1]]) - 1))
  sapply(1:length(data[[1]]), function(i){
    if(i != 2){
      pip.scatter(stats1=data[[1]][[2]]$stats,
                  stats2=data[[1]][[i]]$stats,
                  color1=data[[1]][[2]]$color,
                  color2=data[[1]][[i]]$color,
                  label1=names(data[[1]])[2],
                  label2=names(data[[1]])[i],
                  lpos="topleft")
    }
  })
  dev.off()
}

# Plot PIP comparison histograms
if(names(data[[1]])[1] %in% c("Prior", "Posterior")){
  png(paste(out.prefix, "PIP_comparisons", tolower(names(data[[1]])[1]), "hist.png", sep="."),
      height=2.5*300, width=2.5*(length(data[[1]])-1)*300, res=300)
  par(mfrow=c(1, length(data[[1]]) - 1))
  sapply(2:length(data[[1]]), function(i){
    pip.split.hist(stats1=data[[1]][[1]]$stats,
                   stats2=data[[1]][[i]]$stats,
                   color1=data[[1]][[1]]$color,
                   color2=data[[1]][[i]]$color,
                   label1=names(data[[1]])[1],
                   label2=names(data[[1]])[i])
  })
  dev.off()
}
if(names(data[[1]])[2] %in% c("Prior", "Posterior")){
  png(paste(out.prefix, "PIP_comparisons", tolower(names(data[[1]])[2]), "hist.png", sep="."),
      height=2.5*300, width=2.5*(length(data[[1]])-1)*300, res=300)
  par(mfrow=c(1, length(data[[1]]) - 1))
  sapply(1:length(data[[1]]), function(i){
    if(i != 2){
      pip.split.hist(stats1=data[[1]][[2]]$stats,
                     stats2=data[[1]][[i]]$stats,
                     color1=data[[1]][[2]]$color,
                     color2=data[[1]][[i]]$color,
                     label1=names(data[[1]])[2],
                     label2=names(data[[1]])[i])
    }
  })
  dev.off()
}

# Plot ROCs
sapply(1:length(data), function(i){
  sanitized.name <- gsub('_$', '', gsub('[_]+', "_", gsub('[\ |(|)|&|/|-]', "_", tolower(names(data)[i]))))
  pdf(paste(out.prefix, sanitized.name, "roc.pdf", sep="."), height=4, width=4)
  plot.roc(data[[i]], title=names(data)[i])
  dev.off()
})

# Plot PRCs
sapply(1:length(data), function(i){
  sanitized.name <- gsub('_$', '', gsub('[_]+', "_", gsub('[\ |(|)|&|/|-]', "_", tolower(names(data)[i]))))
  pdf(paste(out.prefix, sanitized.name, "prc.pdf", sep="."), height=4, width=4)
  plot.prc(data[[i]], title=names(data)[i])
  dev.off()
})

# Plot enrichments
sapply(1:length(data), function(i){
  sanitized.name <- gsub('_$', '', gsub('[_]+', "_", gsub('[\ |(|)|&|/|-]', "_", tolower(names(data)[i]))))
  pdf(paste(out.prefix, sanitized.name, "enrichment.pdf", sep="."), height=3, width=4.5)
  plot.calibration(data[[i]], title=names(data)[i])
  dev.off()
})

