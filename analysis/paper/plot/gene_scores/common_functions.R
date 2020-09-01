#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Common functions for gene dosage sensitivity score analyses in rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load dataframe of all gene scores
load.scores <- function(scores.in){
  scores <- read.table(scores.in, header=T, sep="\t", comment.char="", check.names=F)
  colnames(scores) <- c("gene", "pHI", "pTS")
  scores[, -1] <- apply(scores[, -1], 2, as.numeric)
  return(scores)
}

# Load gene metadata from .bed
load.gene.metadata <- function(cov.in){
  cov <- read.table(cov.in, header=T, sep="\t", comment.char="")[, -c(1:3)]
  cov[, -1] <- apply(cov[, -1], 2, as.numeric)
  return(cov)
}

# Load gene feature category metadata from .tsv
load.gene.feature.metadata <- function(meta.in){
  feat.meta <- read.table(meta.in, header=T, sep="\t", comment.char="")
  colnames(feat.meta)[1] <- gsub("X.", "", colnames(feat.meta)[1], fixed=T)
  return(feat.meta)
}

# Partition genes into subgroups based on scores
classify.genes <- function(scores, hc.cutoff=0.9, lc.cutoff=0.5){
  ds.hc <- scores$gene[which(scores$pHI>=hc.cutoff & scores$pTS>=hc.cutoff)]
  ds.lc <- scores$gene[which(scores$pHI>=lc.cutoff & scores$pTS>=lc.cutoff & !(scores$gene %in% ds.hc))]
  hi.hc <- scores$gene[which(scores$pHI>=hc.cutoff & scores$pTS<lc.cutoff)]
  hi.lc <- scores$gene[which(scores$pHI>=lc.cutoff & scores$pTS<lc.cutoff & !(scores$gene %in% hi.hc))]
  ts.hc <- scores$gene[which(scores$pHI<lc.cutoff & scores$pTS>=hc.cutoff)]
  ts.lc <- scores$gene[which(scores$pHI<lc.cutoff & scores$pTS>=lc.cutoff & !(scores$gene %in% ts.hc))]
  ns <- scores$gene[which(scores$pHI<lc.cutoff & scores$pTS<lc.cutoff)]
  list("ds.hc" = ds.hc, "ds.lc" = ds.lc,
       "hi.hc" = hi.hc, "hi.lc" = hi.lc,
       "ts.hc" = ts.hc, "ts.lc" = ts.lc,
       "ns" = ns)
}

# Get gene color based on membership in ds.groups
get.gene.color.byscore <- function(gene, ds.groups){
  if(gene %in% ds.groups$ds.hc){
    cnv.colors[3]
  }else if(gene %in% ds.groups$ds.lc){
    control.cnv.colors[3]
  }else if(gene %in% ds.groups$hi.hc){
    cnv.colors[1]
  }else if(gene %in% ds.groups$hi.lc){
    control.cnv.colors[1]
  }else if(gene %in% ds.groups$ts.hc){
    cnv.colors[2]
  }else if(gene %in% ds.groups$ts.lc){
    control.cnv.colors[2]
  }else{
    ns.color
  }
}

# Compute ROC of a single score vs. predefined true/false genes
roc <- function(stats, score, true.genes, false.genes, steps=seq(1, 0, -0.001)){
  x <- data.frame("score" = stats[, which(colnames(stats) == score)],
                  "true" = stats$gene %in% true.genes,
                  "false" = stats$gene %in% false.genes)
  roc_res <- as.data.frame(t(sapply(steps, function(k){
    idxs <- which(x$score >= k)
    ftrue <- length(which(x$true[idxs])) / length(which(x$true))
    ffalse <- length(which(x$false[idxs])) / length(which(x$false))
    fother <- length(which(!x$true[idxs])) / length(which(!x$true))
    fall <- length(idxs) / nrow(x)
    return(c(k, fall, fother, ftrue, ffalse))
  })))
  roc_res <- rbind(c(-Inf, 0, 0, 0, 0),
                   roc_res,
                   c(Inf, 1, 1, 1, 1))
  colnames(roc_res) <- c("min_score", "frac_all", "frac_other", "frac_true", "frac_false")
  return(roc_res)
}

# Compute PRC of a single score vs. predefined true/false genes
prc <- function(stats, score, true.genes, false.genes, steps=seq(1, 0, -0.001)){
  x <- data.frame("score" = stats[, which(colnames(stats) == score)],
                  "true" = stats$gene %in% true.genes,
                  "false" = stats$gene %in% false.genes)
  prc_res <- as.data.frame(t(sapply(steps, function(k){
    idxs <- which(x$score >= k)
    prec <- length(which(x$true[idxs])) / (length(which(x$true[idxs])) + length(which(x$false[idxs])))
    recall <- length(which(x$true[idxs])) / length(which(x$true))
    fall <- length(idxs) / nrow(x)
    return(c(k, fall, prec, recall))
  })))
  prc_res <- rbind(c(-Inf, 0, 1, 0),
                   prc_res)
  colnames(prc_res) <- c("min_score", "frac_all", "precision", "recall")
  return(prc_res)
}

# Wrapper to calculate all performance stats a single score
evaluate.score <- function(stats, score, true.genes, false.genes){
  roc.res <- roc(stats, score, true.genes, false.genes)
  roc.auc <- flux::auc(roc.res$frac_false, roc.res$frac_true)
  prc.res <- prc(stats, score, true.genes, false.genes)
  prc.auc <- flux::auc(prc.res$recall, prc.res$precision)
  return(list("roc"=roc.res,
              "roc.auc"=roc.auc,
              "prc"=prc.res,
              "prc.auc"=prc.auc))
}

# Load & normalize other (non-rCNV) scores
load.other.scores <- function(meta.in){
  other.scores <- c("gnomad_pLI", "gnomad_oe_mis_upper", "gnomad_oe_lof_upper",
                    "exac_cnv_z", "rvis_pct", "eds", "hurles_hi")
  meta <- read.table(meta.in, header=T, sep="\t", comment.char="")[, c("gene", other.scores)]
  oeuf.idxs <- grep("_upper", colnames(meta))
  meta[, oeuf.idxs] <- apply(meta[, oeuf.idxs], 2, function(vals){
    vals <- vals - min(vals, na.rm=T)
    vals <- vals / max(vals, na.rm=T)
    vals <- 1 - vals
  })
  meta$exac_cnv_z <- meta$exac_cnv_z - min(meta$exac_cnv_z, na.rm=T)
  meta$exac_cnv_z <- meta$exac_cnv_z / max(meta$exac_cnv_z, na.rm=T)
  meta$rvis_pct <- (100 - meta$rvis_pct) / 100
  meta$eds <- meta$eds - min(meta$eds, na.rm=T)
  meta$eds <- meta$eds / max(meta$eds, na.rm=T)
  return(meta)
}

# Modify gene features to improve interpretability of regression coefficients
mod.features <- function(feats){
  # Remove PCA-based features (uninterpretable)
  remove.feats <- colnames(feats)[grep("_component_", colnames(feats), fixed=T)]
  if(length(remove.feats) > 0){
    feats <- feats[, -which(colnames(feats) %in% remove.feats)]
  }
  
  # Mirror some features where a smaller value = more important
  reverse.feats <- c(colnames(feats)[grep("_oe_", colnames(feats), fixed=T)],
                     "rvis", "rvis_pct")
  feats[, reverse.feats] <- -feats[, reverse.feats]
  
  # Helper function to return list of highly correlated features (optional)
  get.correlated.feat.pairs <- function(feats, min.cor=0.7){
    feat.names <- colnames(feats)[-1]
    cor.mat <- cor(feats[, -1])
    cor.pairs <- data.frame("feat1"=character(), "feat2"=character(), "r"=numeric())
    for(ridx in 1:nrow(cor.mat)){
      for(cidx in 1:ncol(cor.mat)){
        fname.r <- rownames(cor.mat)[ridx]
        fname.c <- colnames(cor.mat)[cidx]
        cor.stat <- as.numeric(cor.mat[ridx, cidx])
        if(cor.stat >= min.cor &
           fname.r != fname.c &
           cidx > ridx){
          cor.pairs <- rbind(cor.pairs, c(fname.r, fname.c, cor.stat))
        }
      }
    }
    colnames(cor.pairs) <- c("feat1", "feat2", "r")
    return(cor.pairs)
  }
  
  # Selectively remove highly correlated features
  cor.feats.to.prune <- c("median_expression_q1", "median_expression_q3",
                          "expression_mad_min", "expression_mad_q1", 
                          "expression_mad_q3", "expression_mad_max", 
                          "expression_mad_sd", 
                          intersect(colnames(feats)[grep("chromhmm_", colnames(feats), fixed=T)],
                                    colnames(feats)[grep("_sd", colnames(feats), fixed=T)]),
                          colnames(feats)[grep("gnomad_mu", colnames(feats), fixed=T)],
                          "cen_dist", "tel_dist",
                          "gnomad_oe_mis", "gnomad_oe_lof", "rvis_pct")
  feats <- feats[, which(!colnames(feats) %in% cor.feats.to.prune)]
  
  return(feats)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot distributions of values by gene score category
plot.feature.bydsgroup <- function(feats, ds.groups, feat.idx=2, title=NULL, 
                                   swarm.max=5000, max.ylim=6, 
                                   parmar=c(2.25, 2.75, 1.25, 0.5)){
  # Get plot values
  vals <- lapply(ds.groups, function(genes){
    as.numeric(feats[which(feats$gene %in% genes), feat.idx])
  })
  names(vals) <- names(ds.groups)
  vals <- list("hi.hc"=vals$hi.hc, "ds.hc"=vals$ds.hc, "ts.hc"=vals$ts.hc, 
               "hi.lc"=vals$hi.lc, "ds.lc"=vals$ds.lc, "ts.lc"=vals$ts.lc,
               "ns"=vals$ns)
  ylims <- range(unlist(vals), na.rm=T)
  if(ylims[1] < -max.ylim){ylims[1] <- -max.ylim}
  if(ylims[2] > max.ylim){ylims[2] <- max.ylim}
  x.at <- c(1:3, 5:7, 9)-0.5
  plot.colors <- c(cnv.colors[c(1, 3, 2)], control.cnv.colors[c(1, 3, 2)], ns.color)
  border.colors <- c(cnv.blacks[c(1, 3, 2)], cnv.colors[c(1, 3, 2)], "gray30")
  black.colors <- c(rep(cnv.blacks[c(1, 3, 2)], 2), "gray30")
  lab.colors <- c(plot.colors[1:6], "black")
  pt.cex <- c(rep(0.12, 3), rep(0.075, 4))
  if(is.null(title)){
    title <- colnames(feats)[feat.idx]
  }
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, 9), ylim=ylims,
       xaxt="n", xlab="", yaxt="n", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       bty="n", border=NA, col=bluewhite)
  y.ax.at <- sort(unique(round(axTicks(2))))
  if(length(y.ax.at) > 6){
    y.ax.at <- y.ax.at[seq(1, length(y.ax.at), 2)]
  }
  abline(h=0, col="white", lwd=2)
  # rect(xleft=c(par("usr")[1], 3.1, 7.1, 9.1), 
  #      xright=c(-0.1, 3.9, 7.9, par("usr")[2]), 
  #      ybottom=par("usr")[3], ytop=par("usr")[4],
  #      bty="n", border=NA, col="white")
  
  # Add axes
  x.labs <- c(rep(c("HI", "DS", "TS"), 2), "NS")
  sapply(1:7, function(i){
    axis(1, at=x.at[i], tick=F, line=-0.9,
         labels=x.labs[i], col.axis=lab.colors[i])
  })
  axis(1, at=c(0.1, 2.9), tck=0, col=blueblack, labels=NA, line=1.1)
  axis(1, at=c(4.1, 6.9), tck=0, col=blueblack, labels=NA, line=1.1)
  axis(1, at=1.5, line=0.1, labels="High Conf.", tick=F)
  axis(1, at=5.5, line=0.1, labels="Low Conf.", tick=F)
  axis(2, at=c(-10e10, 10e10), tck=0, col=blueblack)
  axis(2, at=y.ax.at, tck=-0.03, labels=NA, col=blueblack)
  axis(2, at=y.ax.at, tick=F, las=2, line=-0.65)
  mtext(2, text="Feature Value (Z-Score)", line=1.7)
  mtext(3, text=title)
  
  # Add values
  sapply(1:7, function(i){
    ivals <- vals[[i]][which(!is.na(vals[[i]]))]
    if(length(unique(ivals)) > 2){
      if(length(ivals) <= swarm.max){
        beeswarm(ivals, at=x.at[i], col=plot.colors[i], cex=pt.cex[i], add=T,
                 corral="wrap", corralWidth=0.8)
      }else{
        vioplot(ivals, at=x.at[i], col=plot.colors[i], border=border.colors[i],
                add=T, drawRect=F, h=0.1)
      }
      vmed <- median(ivals)
      segments(x0=x.at[i]+0.2, x1=x.at[i]-0.2, y0=vmed, y1=vmed, 
               lend="round", col=black.colors[i], lwd=2)
    }
  })
}

# Superimposed barplot of one statistic per category from evaluate.scores()
superimposed.barplot <- function(values, colors, xleft, xright, ybottom, ytop, 
                                 min.value=NULL, max.value=NULL, title=NULL, 
                                 add.labels=TRUE, lab.cex=0.6, lab.colors=NULL, 
                                 buffer=0.025){
  # Ådd subpanel
  rect(xleft=xleft+buffer, xright=xright-buffer, ybottom=ybottom+buffer, ytop=ytop-buffer, 
       col=bluewhite, border=blueblack)
  text(x=mean(c(xleft, xright)), y=ytop-(2*buffer), labels=title, pos=3)
  
  # Set parameters & get plot values
  inner.xleft <- xleft + (2*buffer)
  inner.xright <- xright - (2*buffer)
  inner.ybottom <- ybottom + (2*buffer)
  inner.ytop <- ytop - (2*buffer)
  n.bars <- length(values)
  bar.buffer <- 0.1
  bar.height <- (inner.ytop - inner.ybottom) / (n.bars + 1)
  bar.y.breaks <- seq(inner.ybottom, inner.ytop, length.out=n.bars + 1)
  bar.ybottom <- bar.y.breaks[1:n.bars] + (bar.buffer * bar.height)
  bar.ytop <- bar.y.breaks[2:(n.bars+1)] - (bar.buffer * bar.height)
  max.bar.length <- inner.xright - inner.xleft
  bar.mid <- mean(c(inner.xright, inner.xleft))
  centered.values <- values - min.value
  norm.values <- centered.values * (max.value / (max.value - min.value))
  scaled.values <- norm.values * max.bar.length
  bar.xleft <- rep(inner.xleft, times=n.bars)
  bar.xright <- inner.xleft + scaled.values
  
  # Add background gridlines
  gridlines.x.at <- seq(inner.xleft, inner.xright, length.out=6)
  segments(x0=gridlines.x.at, x1=gridlines.x.at,
           y0=inner.ybottom, y1=inner.ytop,
           col="white")
  
  # Add bars & labels
  rect(xleft=bar.xleft, xright=bar.xright, ybottom=bar.ybottom, ytop=bar.ytop, border=NA, col=colors)
  if(is.null(lab.colors)){
    lab.colors <- rep("black", n.bars)
    lab.colors[1:ceiling(n.bars/2)] <- "white"
  }
  lab.pos <- sapply(norm.values, function(x){if(x>=0.5){2}else{4}})
  if(any(lab.pos==4)){
    lab.colors[which(lab.pos==4)] <- "black"
  }
  lab.x.adj <- sapply(norm.values, function(x){if(x>=0.5){2*buffer}else{-2*buffer}})
  lab.y.at <- (bar.y.breaks[1:n.bars]+bar.y.breaks[-1])/2
  text(x=bar.xright+lab.x.adj, y=lab.y.at - (bar.buffer * bar.height), pos=lab.pos, 
       col=lab.colors, labels=round(values, 3), cex=lab.cex)
  
  # Add cleanup line
  segments(x0=inner.xleft, x1=inner.xleft, y0=inner.ybottom, y1=inner.ytop, col=blueblack)
}

# Plot ROC curves from a list of evaluate.score() outputs
plot.roc <- function(data, colors=NULL, nested.auc=TRUE, auc.text.colors=NULL, 
                     ax.tick=-0.025, parmar=c(2.5, 2.5, 0.75, 0.75)){
  # Get plot data
  if(is.null(colors)){
    colors <- rev(viridis(length(data)))
  }
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, 1), ylim=c(0, 1), type="n",
       xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlab="", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2], 
       ybottom=par("usr")[3], ytop=par("usr")[4],
       border=NA, bty="n", col="white")
  abline(h=axTicks(2), v=axTicks(1), col=bluewhite)
  abline(0, 1, col=bluewhite, lwd=2, lty=2)
  box(col=bluewhite, bty="o", xpd=T)
  
  # Add curves
  lorder <- order(-sapply(data, function(x){x$roc.auc}))
  sapply(rev(lorder), function(i){
    x <- data[[i]]
    points(x$roc$frac_false, x$roc$frac_true,
           type="l", lwd=3, col=colors[i], xpd=T)
  })
  
  # Add nested AUC subpanel
  if(nested.auc==TRUE){
    superimposed.barplot(rev(sapply(data, function(l){l$roc.auc})), rev(colors),
                         xleft=0.5, xright=1, ybottom=0, ytop=0.5, 
                         min.value=0, max.value=1, title="AUC",
                         lab.colors=auc.text.colors, buffer=0.02)
  }
  
  # Add axes
  axis(1, labels=NA, col=blueblack, tck=ax.tick)
  axis(1, tick=F, line=-0.6)
  mtext(1, line=1.25, text="False positive rate")
  axis(2, labels=NA, col=blueblack, tck=ax.tick)
  axis(2, tick=F, line=-0.6, las=2)
  mtext(2, line=1.65, text="True positive rate")
}

# Plot PRC curves from a list of evaluate.score() outputs
plot.prc <- function(data, colors=NULL, nested.auc=TRUE, auc.text.colors=NULL,
                     ax.tick=-0.025, parmar=c(2.5, 2.5, 0.75, 0.75)){
  # Get plot data
  if(is.null(colors)){
    colors <- rev(viridis(length(data)))
  }
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, 1), ylim=c(0, 1), type="n",
       xaxt="n", yaxt="n", xaxs="i", yaxs="i", xlab="", ylab="")
  rect(xleft=par("usr")[1], xright=par("usr")[2], 
       ybottom=par("usr")[3], ytop=par("usr")[4],
       border=NA, bty="n", col="white")
  abline(h=axTicks(2), v=axTicks(1), col=bluewhite)
  box(col=bluewhite, bty="o", xpd=T)
  
  # Add curves
  lorder <- order(-sapply(data, function(x){x$prc.auc}))
  sapply(rev(lorder), function(i){
    x <- data[[i]]
    points(x$prc$recall, x$prc$precision,
           type="l", lwd=4, col=colors[i], xpd=T)
  })
  
  # Add nested AUC subpanel
  if(nested.auc==TRUE){
    superimposed.barplot(rev(sapply(data, function(l){l$prc.auc})), rev(colors),
                         xleft=0, xright=0.5, ybottom=0, ytop=0.5, 
                         min.value=0, max.value=1, title="AUC",
                         lab.colors=auc.text.colors, buffer=0.02)
  }
  
  # Add axes & cleanup
  axis(1, labels=NA, col=blueblack, tck=ax.tick)
  axis(1, tick=F, line=-0.6)
  mtext(1, line=1.25, text="Precision")
  axis(2, labels=NA, col=blueblack, tck=ax.tick)
  axis(2, tick=F, line=-0.6, las=2)
  mtext(2, line=1.65, text="Recall")
}

# Plot simple color legend
simple.legend <- function(labels, colors){
  par(mar=rep(0.25, 4), bty="n")
  plot(NA, xlim=c(0, 1), ylim=c(0, length(labels)),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  points(x=rep(0.1, length(labels)), y=(1:length(labels))-0.5, pch=22, 
         col=blueblack, bg=colors, cex=1.8)
  text(x=rep(0.1, length(labels)), y=(1:length(labels))-0.58, pos=4, labels=labels, xpd=T)
}


