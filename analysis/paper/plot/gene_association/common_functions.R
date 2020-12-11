#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Common functions for gene association and fine-mapping analyses in rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load credible sets from BED
load.credsets <- function(credsets.in){
  # Read data & clean columns
  credsets <- read.table(credsets.in, header=T, sep="\t", comment.char="")
  colnames(credsets)[1] <- gsub("X.", "", colnames(credsets)[1], fixed=T)
  
  # Ensure numerics
  numeric.cols <- c("start", "end", "mean_control_freq", "mean_case_freq",
                    "pooled_ln_or", "pooled_ln_or_ci_lower", "pooled_ln_or_ci_upper",
                    "best_pvalue", "n_genes")
  credsets[, numeric.cols] <- apply(credsets[, numeric.cols], 2, as.numeric)
  
  # Split list-style columns
  credsets$all_genes <- strsplit(credsets$all_genes, split=";")
  credsets$vconf_genes <- strsplit(credsets$vconf_genes, split=";")
  credsets$conf_genes <- strsplit(credsets$conf_genes, split=";")
  
  return(credsets)
}

# Load all individual gene associations from BED
load.associations <- function(assocs.in){
  # Read data & clean columns
  assocs <- read.table(assocs.in, header=T, sep="\t", comment.char="")
  colnames(assocs)[1] <- gsub("X.", "", colnames(assocs)[1], fixed=T)
  
  # Ensure numerics
  numeric.cols <- c("start", "end", "control_freq", "case_freq",
                    "ln_or", "ln_or_ci_lower", "ln_or_ci_upper",
                    "pvalue", "pip")
  assocs[, numeric.cols] <- apply(assocs[, numeric.cols], 2, as.numeric)
  
  return(assocs)
}

# Categorize genes based on their rank in credible set & PIP
categorize.genes <- function(credsets){
  g.cats <- lapply(c("DEL", "DUP"), function(cnv){
    cnv.idxs <- which(credsets$cnv==cnv)
    g.top <- unique(sort(credsets$top_gene[cnv.idxs]))
    g.nottop <- setdiff(unique(sort(unlist(credsets$all_genes[cnv.idxs]))),
                        g.top)
    g.vconf <- unique(sort(unlist(credsets$vconf_genes[cnv.idxs])))
    g.conf <- setdiff(unique(sort(unlist(credsets$conf_genes[cnv.idxs]))),
                      g.vconf)
    g.notconf <- setdiff(setdiff(unique(sort(unlist(credsets$all_genes[cnv.idxs]))),
                                 g.vconf), g.conf)
    g.top.vconf <- intersect(g.top, g.vconf)
    g.nottop.vconf <- intersect(g.nottop, g.vconf)
    g.top.conf <- intersect(g.top, g.conf)
    g.nottop.conf <- intersect(g.nottop, g.conf)
    g.top.notconf <- intersect(g.top, g.notconf)
    g.nottop.notconf <- intersect(g.nottop, g.notconf)
    list("top.vconf"=g.top.vconf, "nottop.vconf"=g.nottop.vconf,
         "top.conf"=g.top.conf, "nottop.conf"=g.nottop.conf,
         "top.notconf"=g.top.notconf, "nottop.notconf"=g.nottop.notconf)
  })
  names(g.cats) <- c("DEL", "DUP")
  return(g.cats)
}

# Get all genes for a single quadrant of 2x2 grid
get.quadrant.genes <- function(gene.groups, top, conf){
  group.name <- paste(top, conf, sep=".")
  sort(unique(c(gene.groups$DEL[[group.name]],
                gene.groups$DUP[[group.name]])))
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Generic credsets scatterplot function
credsets.scatter <- function(credsets, x, y, subset_to_regions=NULL,
                         xlims=NULL, ylims=NULL, add.lm=T, pt.cex=1,
                         horiz.lines.at=NULL, horiz.lines.lty=1, horiz.lines.color=NULL,
                         abline.a=NULL, abline.b=NULL, abline.lty=1,
                         xtitle=NULL, x.title.line=1.5, x.at=NULL, x.labs=NULL, x.labs.at=NULL, parse.x.labs=FALSE, x.title.cex=1,
                         ytitle=NULL, y.title.line=1.7, y.at=NULL, y.labs=NULL, y.labs.at=NULL, parse.y.labs=FALSE, y.title.cex=1,
                         parmar=c(2.7, 2.7, 0.25, 0.25)){
  # Get plot values
  if(!is.null(subset_to_regions)){
    keep.idx <- which(credsets$region_id %in% subset_to_regions)
    credsets <- credsets[keep.idx, ]
    x <- x[keep.idx]
    y <- y[keep.idx]
  }
  if(is.null(xlims)){
    xlims <- range(x[which(!is.infinite(x))], na.rm=T)
  }
  if(is.null(ylims)){
    ylims <- range(y[which(!is.infinite(y))], na.rm=T)
  }
  del.idx <- which(credsets$cnv=="DEL")
  dup.idx <- which(credsets$cnv=="DUP")
  
  # Prep plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=xlims, ylim=ylims, xlab="", ylab="", xaxt="n", yaxt="n")
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=par("usr")[3], ytop=par("usr")[4],
       border=NA, bty="n", col=bluewhite)
  
  # Add gridlines
  if(is.null(x.at)){
    x.at <- axTicks(1)
  }
  if(is.null(y.at)){
    y.at <- axTicks(2)
  }
  abline(v=x.at, h=y.at, col="white")
  if(!is.null(horiz.lines.at)){
    n.horiz.lines <- length(horiz.lines.at)
    if(is.null(horiz.lines.color) | length(horiz.lines.color) != n.horiz.lines){
      horiz.lines.color <- rep(blueblack, n.horiz.lines)
    }
    sapply(1:length(horiz.lines.at), function(i){
      abline(h=horiz.lines.at[i], lty=horiz.lines.lty[i], col=horiz.lines.color[i])
    })
  }
  if(!is.null(abline.a)){
    sapply(1:length(abline.a), function(i){
      abline(a=abline.a[i], b=abline.b[i], lty=abline.lty[i], col=blueblack)
    })
  }
  
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
  
  # Add points (always add gw-sig last)
  points(x, y, pch=21, cex=pt.cex,
         bg=cnv.colors[credsets$cnv], col=cnv.blacks[credsets$cnv])
  
  # Add axis ticks
  axis(1, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(1, at=x.at, labels=NA, tck=-0.03, col=blueblack)
  axis(2, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
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
  mtext(1, text=xtitle, line=x.title.line, cex=x.title.cex)
  mtext(2, text=ytitle, line=y.title.line, cex=y.title.cex)
}

