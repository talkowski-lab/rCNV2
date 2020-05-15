#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Common functions for large segment analyses in rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load summary dataframe for locus-level association stats
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
  loci$formatted_name <- sapply(strsplit(loci$region_id, split="_"), function(parts){parts[4]})
  loci$formatted_size <- paste(prettyNum(round(loci$cred_intervals_size/1000, 0), big.mark=","), "kb", sep=" ")
  
  return(loci)
}

# Load summary table for all segments (including lit. GDs & NAHR regions)
load.segment.table <- function(segs.in){
  # Read data
  segs <- read.table(segs.in, header=T, sep="\t", comment.char="")
  colnames(segs)[1] <- "chr"
  
  # Split list-style columns
  listcol.idxs <- c(which(colnames(segs) %in% c("coords", "genes")),
                    intersect(grep("^n_", colnames(segs), invert=T, fixed=F),
                            grep("_genes$", colnames(segs), fixed=F)))
  segs[, listcol.idxs] <- apply(segs[, listcol.idxs], 2, strsplit, split=";")
  
  # Convert numeric columns to numerics
  numcol.idxs <- c(which(colnames(segs) %in% c("start", "end", "size")),
                   grep("^n_", colnames(segs), fixed=F))
  segs[, numcol.idxs] <- apply(segs[, numcol.idxs], 2, as.numeric)
  
  # Convert boolean dummy columns to logicals
  boolcol.idxs <- which(colnames(segs) %in% c("gw_sig", "hc_gd", "mc_gd", "lc_gd", "any_gd",
                                              "pathogenic", "benign", "nahr", "pleiotropic"))
  segs[, boolcol.idxs] <- apply(segs[, boolcol.idxs], 2, function(vals){
    sapply(vals, function(val){if(val==1){TRUE}else{FALSE}})})
  
  # Add formatted sizes
  segs$formatted_size <- paste(prettyNum(round(segs$size/1000, 0), big.mark=","), "kb", sep=" ")
  
  # Return cleaned dataframe
  return(segs)
}

# Merge locus association data with master segments table
merge.loci.segs <- function(loci, segs){
  shared.columns <- intersect(colnames(loci), colnames(segs))
  gw <- merge(loci, segs, by=shared.columns, all.x=T, all.y=F, 
              sort=F, suffixes=c(".l", ".s"))
  gw$color <- cnv.colors[sapply(gw$cnv, function(cnv){which(names(cnv.colors)==cnv)})]
  gw$black <- cnv.blacks[sapply(gw$cnv, function(cnv){which(names(cnv.blacks)==cnv)})]
  return(gw)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Generic segment scatterplot function
gw.scatter <- function(gw, x, y, xlims=NULL, ylims=NULL, add.lm=T,
                       xtitle=NULL, x.at=NULL, x.labs=NULL, x.labs.at=NULL, parse.x.labs=FALSE,
                       ytitle=NULL, y.at=NULL, y.labs=NULL, y.labs.at=NULL, parse.y.labs=FALSE,
                       parmar=c(3, 3, 0.8, 0.8)){
  # Get plot values
  if(is.null(xlims)){
    xlims <- range(x[which(!is.infinite(x))], na.rm=T)
  }
  if(is.null(ylims)){
    ylims <- range(y[which(!is.infinite(y))], na.rm=T)
  }
  del.idx <- which(gw$cnv=="DEL")
  dup.idx <- which(gw$cnv=="DUP")
  
  # Prep plot area
  par(mar=parmar)
  plot(NA, xlim=xlims, ylim=ylims, xlab="", ylab="", xaxt="n", yaxt="n")
  
  # Add linear fits, if optioned
  if(add.lm==T){
    del.fit <- robust.lm(x=x[del.idx], y=y[del.idx], conf=0.95)
    dup.fit <- robust.lm(x=x[dup.idx], y=y[dup.idx], conf=0.95)
    polygon(x=c(del.fit$ci$x, rev(del.fit$ci$x)),
            y=c(del.fit$ci$lower, rev(del.fit$ci$upper)),
            border=NA, col=adjustcolor(cnv.colors[1], alpha=0.2))
    polygon(x=c(dup.fit$ci$x, rev(dup.fit$ci$x)),
            y=c(dup.fit$ci$lower, rev(dup.fit$ci$upper)),
            border=NA, col=adjustcolor(cnv.colors[2], alpha=0.2))
    abline(del.fit$fit, lwd=2, col=cnv.colors[1])
    abline(dup.fit$fit, lwd=2, col=cnv.colors[2])
  }
  
  # Add points
  points(x, y, pch=21, bg=gw$color, col=gw$black)
  
  # Add axis ticks
  if(is.null(x.at)){
    x.at <- axTicks(1)
  }
  if(is.null(y.at)){
    y.at <- axTicks(2)
  }
  axis(1, at=x.at, labels=NA, tck=-0.03, col=blueblack)
  axis(2, at=y.at, labels=NA, tck=-0.03, col=blueblack)
  
  # Add axis labels
  if(is.null(x.labs)){
    x.labs <- x.at
  }
  if(is.null(x.labs.at)){
    x.labs.at <- x.at
  }
  if(is.null(y.labs)){
    y.labs <- y.at
  }
  if(is.null(y.labs.at)){
    y.labs.at <- y.at
  }
  sapply(1:length(x.labs.at), function(i){
    if(parse.x.labs==TRUE){
      axis(1, at=x.labs.at[i], labels=parse(text=x.labs[i]), tick=F, line=-0.6)
    }else{
      axis(1, at=x.labs.at[i], labels=x.labs[i], tick=F, line=-0.6)
    }
  })
  sapply(1:length(y.labs.at), function(i){
    if(parse.y.labs==TRUE){
      axis(2, at=y.labs.at[i], labels=parse(text=y.labs[i]), tick=F, line=-0.6, las=2)
    }else{
      axis(2, at=y.labs.at[i], labels=y.labs[i], tick=F, line=-0.6, las=2)
    }
  })
  
  # Add axis titles
  mtext(1, text=xtitle, line=1.75)
  mtext(2, text=ytitle, line=1.75)
  
  # Add cleanup box
  box(col=blueblack, bty="o")
}

# Generic swarm/boxplot function
gw.swarm <- function(gw, x.bool, y, cnv.split=TRUE, ylims=NULL, 
                     add.pvalue=FALSE, stat.test="wilcoxon", alternative="two.sided",
                     xtitle=NULL, x.labs=c("FALSE", "TRUE"),
                     ytitle=NULL, y.at=NULL, y.labs=NULL, y.labs.at=NULL, parse.y.labs=FALSE,
                     parmar=c(2.3, 3, 0.5, 0.5)){
  
  require(beeswarm, quietly=T)
  
  # Get plot values
  if(is.null(ylims)){
    ylims <- range(y[which(!is.infinite(y))], na.rm=T)
  }
  if(cnv.split==TRUE){
    del.idx <- which(gw$cnv=="DEL")
    dup.idx <- which(gw$cnv=="DUP")
    x.at <- c(0.3, 0.7, 1.3, 1.7)
    width <- 0.2
    y.vals <- list(y[intersect(which(!x.bool), del.idx)],
                   y[intersect(which(!x.bool), dup.idx)],
                   y[intersect(which(x.bool), del.idx)],
                   y[intersect(which(x.bool), dup.idx)])
    pt.color.list <- lapply(1:4, function(i){
      rep(rep(cnv.colors, 2)[i], length(y.vals[[i]]))
    })
    pt.border.list <- lapply(1:4, function(i){
      rep(rep(cnv.blacks, 2)[i], length(y.vals[[i]]))
    })
    boxplot.colors <- rep(cnv.colors, 2)
    boxplot.fill <- rep(control.cnv.colors, 2)
  }else{
    x.at <- c(0.5, 1.5)
    width <- 0.4
    y.vals <- list(y[which(!x.bool)],
                   y[which(x.bool)])
    pt.color.list <- list(gw$color[which(!x.bool)],
                          gw$color[which(x.bool)])
    pt.border.list <- list(gw$black[which(!x.bool)],
                          gw$black[which(x.bool)])
    boxplot.colors <- rep(blueblack, 2)
    boxplot.fill <- rep(bluewhite, 2)
  }
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, 2), ylim=ylims, type="n", xlab="", ylab="", xaxt="n", yaxt="n")
  
  # Add boxplots
  boxplot(y.vals, at=x.at, outline=F, lty=1, add=T,
          outwex=width, staplewex=width, boxwex=width,
          border=boxplot.colors, col=boxplot.fill, 
          xaxt="n", yaxt="n")
  
  # Add swarms
  sapply(1:length(y.vals), function(i){
    beeswarm(y.vals[[i]], add=T, at=x.at[i], pch=21,
             pwbg=pt.color.list[[i]], pwcol=pt.border.list[[i]],
             corral="wrap", corralWidth=width)
  })
  
  # Add x-axis
  sapply(1:2, function(x){
    axis(1, at=x-c(0.1, 0.9), tck=0, labels=NA, col=blueblack)
    axis(1, at=x-0.5, line=-0.9, labels=x.labs[x], tick=F)
  })
  axis(1, at=c(0.1, 1.9), tck=0, labels=NA, line=1.2, col=blueblack)
  mtext(1, line=1.3, text=xtitle)
  
  # Add y-axis ticks
  if(is.null(y.at)){
    y.at <- axTicks(2)
  }
  axis(2, at=y.at, labels=NA, tck=-0.03, col=blueblack)
  if(is.null(y.labs)){
    y.labs <- y.at
  }
  if(is.null(y.labs.at)){
    y.labs.at <- y.at
  }
  sapply(1:length(y.labs.at), function(i){
    if(parse.y.labs==TRUE){
      axis(2, at=y.labs.at[i], labels=parse(text=y.labs[i]), tick=F, line=-0.6, las=2)
    }else{
      axis(2, at=y.labs.at[i], labels=y.labs[i], tick=F, line=-0.6, las=2)
    }
  })
  mtext(2, text=ytitle, line=1.75)
  
  # Add P-values, if optioned
  if(stat.test=="wilcoxon"){
    stat.res <- wilcox.test(y ~ x.bool, alternative=alternative)
    pval <- stat.res$p.value
    print(stat.res)
  }
  if(add.pvalue==T){
    if(cnv.split==T){
      warning("P-value labeling with cnv.split=T is currently unsupported")
    }else{
      axis(3, at=c(0.5, 1.5), tck=0.03, col=blueblack, labels=NA, line=0.5)
      mtext(3, text=format.pval(pval), line=0.5)
    }
  }
}
