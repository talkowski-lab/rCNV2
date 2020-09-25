#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Common functions for locus highlight plots in rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Extract gene components from a GTF
load.genes <- function(gtf.in, region){
  # Required for bedr in local Rstudio only:
  # Sys.setenv(PATH = paste(Sys.getenv("PATH"),
  #                         "/Users/collins/anaconda3/envs/py3/bin",
  #                         sep = ":"))
  
  # Tabix region of interest
  if(!file.exists(paste(gtf.in, "tbi", sep="."))){
    stop(paste("tabix index not found for input file", gtf.in))
  }
  require(bedr, quietly=T)
  gtf <- bedr::tabix(region, gtf.in, check.chr=FALSE)
  
  # Reformat entries
  if(!is.null(gtf)){
    colnames(gtf) <- c("chr", "source", "feature", "start", "end", "score", 
                       "strand", "frame", "attribute")
    gtf$gene <- sapply(gtf$attribute, function(atrs.str){
      atrs <- unlist(strsplit(atrs.str, split=";"))
      gsub("\"", "", unlist(strsplit(atrs[grep("gene_name", atrs)], split=" "))[[2]])
    })
    gtf <- gtf[, c("chr", "start", "end", "gene", "strand", "feature")]
    gtf[, c("start", "end")] <- apply(gtf[, c("start", "end")], 2, as.numeric)
  }else{
    gtf <- data.frame("chr"=character(), "start"=numeric(), "end"=numeric(),
                      "gene"=character(), "strand"=character(), "feature"=character())
  }
  
  return(gtf)
}

# Extract summary statistics from a single BED
load.sumstats <- function(bedpath, region){
  # Required for bedr in local Rstudio only:
  # Sys.setenv(PATH = paste(Sys.getenv("PATH"),
  #                         "/Users/collins/anaconda3/envs/py3/bin",
  #                         sep = ":"))
  
  # Tabix region of interest
  if(!file.exists(paste(bedpath, "tbi", sep="."))){
    stop(paste("tabix index not found for input file", bedpath))
  }
  require(bedr, quietly=T)
  ss <- bedr::tabix(region, bedpath, check.chr=FALSE)
  
  # Add midpoint
  ss$pos <- (ss$start + ss$stop)/2
  
  # Return columns of interest
  ss[, c("chr", "pos", "meta_phred_p", "meta_lnOR", "meta_lnOR_lower", "meta_lnOR_upper")]
}

# Extract PIPs from BED for specified genes for a single HPO
load.pips <- function(pips.in, genes, hpo){
  pips <- read.table(pips.in, header=T, sep="\t", comment.char="")
  pips[which(pips[, 1] == hpo & pips$gene %in% genes), c("gene", "PIP")]
}

# Load sample sizes for case/control contrast directly from .tsv
get.sample.sizes <- function(table.in, case.hpos=c("HP:0000118"), 
                             ctrl.hpo="HEALTHY_CONTROL"){
  n <- read.table(table.in, header=T, sep="\t", comment.char="")
  n.case <- apply(n[which(n[, 1] %in% case.hpos), grep("meta", colnames(n), fixed=T)], 2, max, na.rm=T)
  n.ctrl <- n[which(n[, 1] == ctrl.hpo), grep("meta", colnames(n), fixed=T)]
  return(list("case"=as.numeric(as.vector(n.case)), "ctrl"=as.numeric(as.vector(n.ctrl))))
}

# Load CNVs from a single BED, and split by case/control
load.cnvs <- function(bedpath, region, cnv=NULL, case.hpos=c("HP:0000118"), 
                      ctrl.hpo="HEALTHY_CONTROL"){
  # Required for bedr in local Rstudio only:
  # Sys.setenv(PATH = paste(Sys.getenv("PATH"),
  #                         "/Users/collins/anaconda3/envs/py3/bin",
  #                         sep = ":"))
  
  # Tabix region of interest
  if(!file.exists(paste(bedpath, "tbi", sep="."))){
    stop(paste("tabix index not found for input file", bedpath))
  }
  require(bedr, quietly=T)
  cnvs <- bedr::tabix(region, bedpath, check.chr=FALSE)
  if(!is.null(cnvs)){
    # Ensure consistent column names
    colnames(cnvs) <- c("chr", "start", "end", "cnv_id", "cnv", "pheno")
    
    # Sort & filter CNVs
    cnvs <- cnvs[with(cnvs, order(start, end)), ]
    if(!is.null(cnv)){
      cnvs <- cnvs[which(cnvs$cnv==cnv), ]
    }
    case.cnv.idxs <- which(sapply(cnvs$pheno, function(pstr){
      any(case.hpos %in% unlist(strsplit(pstr, split=";", fixed=T)))
    }))
    case.cnvs <- cnvs[case.cnv.idxs, ]
    ctrl.cnvs <- cnvs[grep(ctrl.hpo, cnvs$pheno, fixed=T), ]
  }else{
    empty.df <- data.frame("chr"=character(), "start"=numeric(), "end"=numeric(), 
                           "cnv_id"=character(), "cnv"=character(), "pheno"=character())
    case.cnvs <- empty.df
    ctrl.cnvs <- empty.df
  }
  
  return(list("case"=case.cnvs, "ctrl"=ctrl.cnvs))
}

# Load all CNVs from a list of BEDs
load.cnvs.multi <- function(cnvlist, region, cnv=NULL, case.hpo="HP:0000118", 
                            ctrl.hpo="HEALTHY_CONTROL"){
  # cnvlist should be a tsv of (cohort, path) pairs
  cnv.list <- read.table(cnvlist, header=F, sep="\t")
  cnvs <- lapply(1:nrow(cnv.list), function(i){
    load.cnvs(cnv.list[i, 2], region, cnv, case.hpo, ctrl.hpo)
  })
  names(cnvs) <- cnv.list[, 1]
  return(cnvs)
}

# Transform CNV coordinates to plotting values (with colors)
pileup.cnvs <- function(cnvs, start=NULL, end=NULL, dx=100, 
                        cnv.height=1, cnv.buffer=0.15, bevel.switch.pct=0.025,
                        col=blueblack, highlight.hpo=NA, highlight.col=NULL){
  # Set range of values to evaluate
  if(is.null(start)){
    start <- min(cnvs$start)
  }
  if(is.null(end)){
    end <- max(cnvs$end)
  }
  x <- seq(start, end, length.out=dx+1)
  
  # Set other scaling parameters
  cnv.y.buf <- cnv.height * cnv.buffer
  
  # Build empty dataframe of CNV pileup
  counts <- data.frame("pos"=x, "count"=0, "idx"=1:length(x))
  
  # Create plotting values for each CNV
  # Note: must use for loop to increment counts after each CNV is added
  cnv.plot.values <- list()
  if(nrow(cnvs) > 0){
    for(i in 1:nrow(cnvs)){
      # Get CNV info and increment counts
      cnv.id <- cnvs$cnv_id[i]
      cnv.start <- cnvs$start[i]
      cnv.end <- cnvs$end[i]
      cnv.x.idxs <- which(x >= cnv.start & x <= cnv.end)
      counts$count[cnv.x.idxs] <- counts$count[cnv.x.idxs] + cnv.height
      cnv.hpos <- unlist(strsplit(cnvs$pheno[i], split=";", fixed=T))
      
      # Create plotting vectors for each CNV
      cnv.x <- c(x[cnv.x.idxs], rev(x[cnv.x.idxs]))
      cnv.y <- c(counts$count[cnv.x.idxs] - cnv.y.buf,
                 rev(counts$count[cnv.x.idxs] - cnv.height + cnv.y.buf))
      
      # Bevel edges by single dx
      # Bevel left edge according to the CNV that came before
      max.change.dist <- ceiling(bevel.switch.pct * dx)
      if(i == 1 | cnv.x.idxs[1] == 1){
        dist.to.nearest.change <- 10e10
      }else{
        same.height.idxs <- which(prev.counts$count == max(counts$count[cnv.x.idxs]))
        if(length(same.height.idxs) > 0){
          dist.to.nearest.change <- min(cnv.x.idxs) - max(same.height.idxs)
        }else{
          dist.to.nearest.change <- 10e10
        }
      }
      if(dist.to.nearest.change > max.change.dist){
        cnv.y[1] <- counts$count[cnv.x.idxs[1]] - cnv.height + cnv.y.buf
      }else{
        cnv.y[length(cnv.y)] <- counts$count[cnv.x.idxs[1]] - cnv.y.buf
      }
      # Always bevel right edge sloping outward
      cnv.y[length(cnv.x.idxs)] <- counts$count[cnv.x.idxs[length(cnv.x.idxs)]] - cnv.height + cnv.y.buf
      
      # Assign color
      if(highlight.hpo %in% cnv.hpos){
        cnv.color <- highlight.col
      }else{
        cnv.color <- col
      }
      
      # Add CNV plotting values to output list
      cnv.plot.values[[cnv.id]] <- list("x"=cnv.x, "y"=cnv.y, "color"=cnv.color)
      
      # Save previous CNV's x indexes and counts for comparisons
      prev.x.idxs <- cnv.x.idxs
      prev.counts <- counts
    }
  }
  
  # Increment final counts by cnv.y.buf (for plotting)
  counts$counts  <- counts$count + cnv.y.buf
  
  return(list("cnvs"=cnv.plot.values, "counts"=counts))
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Add stick figure for idiogram
add.idio.stick <- function(genome.in, chrom, start, end, y.at, tick.height=0.1, 
                           wex=0.925, y.label=TRUE){
  # Load & scale genomic coordinates for chromosome of interest
  g <- read.table(genome.in, sep="\t", header=F)
  len <- g[which(g[, 1] == chrom), 2]
  pbuf <- ((1-wex)/2) * diff(par("usr")[1:2])
  pstart <- par("usr")[1] + pbuf
  pend <- par("usr")[2] - pbuf
  
  # Scale highlight coordinates
  iwidth <- pend - pstart
  h.start <- pstart + ((start / len) * iwidth)
  h.end <- pstart + ((end / len) * iwidth)
  
  # Plot idiogram stick + highlight box
  segments(x0=pstart, x1=pend, y0=y.at, y1=y.at, col="gray50")
  rect(xleft=h.start, xright=h.end, ybottom=y.at-tick.height, ytop=y.at+tick.height,
       col="red")
  if(y.label==TRUE){
    axis(2, at=y.at, tick=F, line=-0.85, labels=paste("chr", chrom, sep=""), las=2)
  }
}

# Add rounded coordinate line
add.coord.line <- function(start, end, y0, highlight.start, highlight.end, highlight.col="red", 
                           tick.height=0.1, max.tick=6, vlines=FALSE, vlines.bottom=NULL, lab.cex=4.5/6){
  require(shape, quietly=TRUE)
  
  # Get plot coordinates
  log.base <- floor(log10((end-start)/6))
  tick.space <- ceiling((end-start)/(10^log.base))*(10^(log.base-1))
  tick.at <- seq(0, 10e8, by=tick.space)
  tick.at <- tick.at[which(tick.at >= start & tick.at <= end)]
  if(length(tick.at) > max.tick){
    tick.at <- tick.at[seq(1, length(tick.at), 2)]
    tick.space <- tick.at[2] - tick.at[1]
  }
  n.ticks <- length(tick.at)
  
  # Set tick labels
  if(log.base < 3){
    suffix <- "kb"
    denom <- 1000
    nsmall <- (3 - log.base) + 1
  }else{
    suffix <- "Mb"
    denom <- 1000000
    nsmall <- (6 - log.base) + 1
  }
  tick.labels <- paste(format(round(tick.at / denom, nsmall), nsmall=nsmall), suffix)
  tick.labels.y.at <- rep(y0, n.ticks)
  
  # Draw vertical lines for ticks, if optioned
  if(vlines==TRUE){
    if(is.null(vlines.bottom)){
      vlines.bottom <- par("usr")[3]
    }
    segments(x0=tick.at, x1=tick.at, y0=y0, y1=vlines.bottom, col=bluewhite)
  }
  
  # Draw coordinate line & highlight box
  abline(h=y0, col=blueblack, lwd=2)
  rect(xleft=highlight.start, xright=highlight.end, 
       ybottom=y0 - (0.5 * tick.height), ytop=y0 + (0.5 * tick.height), 
       bty="n", border=NA, col=highlight.col)
  shape::Arrowhead(x0=par("usr")[1:2], y0=rep(y0, 2), angle=c(180, 0), 
                   arr.length=0.3, arr.width=0.2, arr.col=blueblack, arr.type="triangle", 
                   xpd=T, lcol=NA)
  segments(x0=tick.at, x1=tick.at, 
           y0=rep(y0-tick.height, n.ticks), y1=rep(y0+tick.height, n.ticks),
           col=blueblack)
  
  # Add tick labels
  text(x=tick.at, y=tick.labels.y.at, cex=lab.cex, labels=tick.labels, pos=3, xpd=T)
}

# Add gene bodies to plot
add.genes <- function(genes, y0, n.rows=1, panel.height=0.2, col=ns.color,
                      y.axis.title=NULL, y.axis.title.col=NULL, label.genes=c()){
  gene.names <- unique(genes$gene)
  n.genes <- length(gene.names)
  
  # Get y scaling
  panel.bottom <- y0 - (0.5*panel.height)
  panel.top <- y0 + (0.5*panel.height)
  row.height <- panel.height / n.rows
  exon.height <- (2/3) * row.height
  utr.height <- (1/3) * row.height
  row.mids <- seq(panel.bottom + (0.5*row.height), 
                  panel.top - (0.5*row.height), 
                  length.out=n.rows)
  
  # Iterate over genes and plot each one at a time
  sapply(1:n.genes, function(i){
    gene <- gene.names[i]
    row.idx <- (i %% n.rows) + 1
    tx.coords <- genes[which(genes$gene==gene & genes$feature=="transcript"), 2:3]
    tx.coords$end[which(tx.coords$end > par("usr")[2])] <- par("usr")[2]
    tx.coords$start[which(tx.coords$start < par("usr")[1])] <- par("usr")[1]
    ex.coords <- genes[which(genes$gene==gene & genes$feature=="exon"), 2:3]
    utr.coords <- genes[which(genes$gene==gene & genes$feature=="UTR"), 2:3]
    if(nrow(tx.coords) > 0){
      segments(x0=tx.coords[, 1], x1=tx.coords[, 2], 
               y0=row.mids[row.idx], y1=row.mids[row.idx],
               lend="round", col=col)
    }
    if(nrow(ex.coords) > 0){
      rect(xleft=ex.coords[, 1], xright=ex.coords[, 2],
           ybottom=row.mids[row.idx]-(0.5*exon.height),
           ytop=row.mids[row.idx]+(0.5*exon.height),
           border=col, lwd=0.3, col=col)
    }
    if(nrow(utr.coords) > 0){
      rect(xleft=utr.coords[, 1], xright=utr.coords[, 2],
           ybottom=row.mids[row.idx]-(0.5*utr.height),
           ytop=row.mids[row.idx]+(0.5*utr.height),
           border=col, lwd=0.3, col=col)
    }
    if(gene %in% label.genes){
      if(row.mids[row.idx] >= y0){
        glabel.pos <- 3
        glabel.y <- row.mids[row.idx]-(0.25*row.height)
      }else{
        glabel.pos <- 1
        glabel.y <- row.mids[row.idx]+(0.25*row.height)
      }
      text(x=mean(as.numeric(tx.coords[1, 1:2])), y=glabel.y,
           cex=5/6, font=3, labels=gene, pos=glabel.pos)
    }
  })
  
  # Add y-axis title (if optioned)
  if(is.null(y.axis.title.col)){
    y.axis.title.col <- col
  }
  axis(2, at=y0, tick=F, line=-0.9, las=2, labels=y.axis.title, col.axis=y.axis.title.col)
}

# Add square bracket to plot
add.bracket <- function(xleft, xright, y0, height, col=blueblack, staple.wex=0.025, 
                        lwd=1, left.label=NULL){
  # Get various dimensions for plotting
  staple.width <- staple.wex*diff(par("usr")[1:2])
  half.height <- height / 2
  ytop <- y0 - half.height
  ybottom <- y0 + half.height
  
  # Plot staples
  segments(x0=rep(c(xleft, xright), 2), x1=rep(c(xleft-staple.width, xright+staple.width), 2), 
           y0=c(ytop, ytop, ybottom, ybottom), y1=c(ytop, ytop, ybottom, ybottom), 
           col=col, lwd=lwd, lend="round")
  
  # Plot vertical connectors
  segments(x0=c(xleft-staple.width, xright+staple.width), x1=c(xleft-staple.width, xright+staple.width), 
           y0=rep(ybottom, 2), y1=rep(ytop, 2),
           col=col, lwd=lwd, lend="round")
  
  # Add labels, if optioned
  text(x=xleft-staple.width, y=y0, pos=2, col=col, labels=left.label)
}

# Add panel of phred-scaled P-values 
add.pvalues <- function(ss, y0, cnv.type, panel.height=0.2, pt.cex=0.6){
  # Get panel parameters
  half.height <- 0.5*panel.height
  ybottom <- y0 - half.height
  ytop <- y0 + half.height
  
  # Set CNV-based plotting values
  if(cnv.type=="DEL"){
    pt.col <- cnv.blacks[1]
  }else if(cnv.type=="DUP"){
    pt.col <- cnv.blacks[2]
  }
  
  # Scale p-values according to y0 and panel.height
  pos <- as.numeric(ss$pos)
  pvals.orig <- as.numeric(ss$meta_phred_p)
  max.pval.orig <- max(pvals.orig, na.rm=T)
  pvals.scaled <- (panel.height / (ceiling(max.pval.orig) + 1)) * pvals.orig
  pvals <- pvals.scaled + y0 - half.height
  
  # Add horizontal gridlines
  y.ax.tick.spacing <- seq(-half.height, half.height, length.out=6)
  abline(h=c(y0 + y.ax.tick.spacing), col="white")
  
  # Add points
  points(x=pos, y=pvals, pch=19, col=pt.col, cex=pt.cex)
  
  # Add Y-axis
  y.ax.label.cex <- 5/6
  axis(2, at=c(ybottom, ytop), tick=0, labels=NA, col=blueblack)
  axis(2, at=y0 + y.ax.tick.spacing, tck=-0.0075, col=blueblack, labels=NA)
  axis(2, at=y0+c(-half.height, half.height), tick=F, las=2, line=-0.65,
       labels=c(0, ceiling(max.pval.orig)+1), cex.axis=y.ax.label.cex)
  axis(2, at=y0, line=-0.2, tick=F, labels=bquote(-log[10](italic("P")[.(cnv.type)])), las=2)
  
  # Add cleanup top & bottom lines
  abline(h=c(ytop, ybottom), col=blueblack)
  segments(x0=par("usr")[2], x1=par("usr")[2], y0=ybottom, y1=ytop, col=blueblack, xpd=T)
}

# Add panel of natural log-scaled odds ratios 
add.ors <- function(ss, y0, cnv.type, panel.height=0.2, pt.cex=0.7){
  # Get panel parameters
  half.height <- 0.5*panel.height
  ybottom <- y0 - half.height
  ytop <- y0 + half.height
  
  # Set CNV-based plotting values
  if(cnv.type=="DEL"){
    line.col <- cnv.colors[1]
    ci.col <- control.cnv.colors[1]
  }else if(cnv.type=="DUP"){
    line.col <- cnv.colors[2]
    ci.col <- control.cnv.colors[2]
  }
  
  # Scale odds ratios according to y0 and panel.height
  pos <- as.numeric(ss$pos)
  ors.orig <- apply(ss[, grep("meta_lnOR", colnames(ss), fixed=T)], 2, as.numeric)
  ors.scalar <- max(ors.orig[, 1], na.rm=T)
  ors.scaled <- (panel.height / (ceiling(ors.scalar) + 1)) * ors.orig
  ors.scaled[which(ors.scaled > panel.height)] <- panel.height
  ors.scaled[which(ors.scaled < 0)] <- 0
  ors <- ors.scaled + y0 - half.height
  
  # Only keep points with non-NA ORs
  keep.idx <- which(apply(ors, 1, function(vals){all(!is.na(vals))}))
  pos <- pos[keep.idx]
  ors <- ors[keep.idx, ]
  
  # Add horizontal gridlines
  y.ax.tick.spacing <- seq(-half.height, half.height, length.out=6)
  abline(h=c(y0 + y.ax.tick.spacing), col="white")
  
  # Add points & shading
  polygon(x=c(pos, rev(pos)), y=c(ors[, 2], rev(ors[, 3])),
          col=adjustcolor(ci.col, alpha=0.5), border=NA, bty="n")
  points(x=pos, y=ors[, 1], lwd=3, col=line.col, type="l")
  
  # Add Y-axis
  y.ax.label.cex <- 5/6
  axis(2, at=c(ybottom, ytop), tick=0, labels=NA, col=blueblack)
  axis(2, at=y0 + y.ax.tick.spacing, tck=-0.0075, col=blueblack, labels=NA)
  axis(2, at=y0+c(-half.height, half.height), tick=F, las=2, line=-0.65,
       labels=c(0, ceiling(ors.scalar) + 1), cex.axis=y.ax.label.cex)
  axis(2, at=y0, line=-0.2, tick=F, labels=bquote("ln" * ("OR"[.(cnv.type)])), las=2)
  
  # Add cleanup top & bottom lines
  abline(h=c(ytop, ybottom), col=blueblack)
  segments(x0=par("usr")[2], x1=par("usr")[2], y0=ybottom, y1=ytop, col=blueblack, xpd=T)
}

# Add plot of PIPs for all genes
add.pips <- function(pips, genes, y0, panel.height=0.2, label.genes=NULL,
                     col=blueblack, highlight.col="red", highlight.genes=NULL){
  # Get panel parameters
  half.height <- 0.5*panel.height
  ybottom <- y0 - half.height
  ytop <- y0 + half.height
  
  # Merge gene coordinates & PIPs
  gcoords <- genes[which(genes$gene %in% pips$gene & genes$feature=="transcript"), 
                   c("start", "end", "gene")]
  pdat <- merge(gcoords, pips, all=F, sort=F, by="gene")
  
  # Scale PIPs according to y0 and panel.height
  max.pip <- ceiling(11*max(pdat$PIP, na.rm=T)) / 10
  pdat$y <- (pdat$PIP * panel.height / min(c(1, max.pip))) + y0 - half.height
  
  # Add horizontal gridlines
  y.ax.tick.spacing <- seq(-half.height, half.height, length.out=6)
  abline(h=c(y0 + y.ax.tick.spacing), col="white")
  
  # Add Y-axis
  y.ax.label.cex <- 5/6
  axis(2, at=c(ybottom, ytop), tick=0, labels=NA, col=blueblack)
  axis(2, at=y0 + y.ax.tick.spacing, tck=-0.0075, col=blueblack, labels=NA)
  axis(2, at=y0+c(-half.height, half.height), tick=F, las=2, line=-0.65,
       labels=c(0, max.pip), cex.axis=y.ax.label.cex)
  axis(2, at=y0, line=-0.2, tick=F, labels="PIP", las=2)
  
  # Add cleanup top & bottom lines
  abline(h=c(ytop, ybottom), col=blueblack)
  segments(x0=par("usr")[2], x1=par("usr")[2], y0=ybottom, y1=ytop, col=blueblack, xpd=T)
  
  # Add segments for each gene
  if(!is.null(highlight.genes)){
    highlight.idxs <- which(pdat$gene %in% highlight.genes)
    other.idxs <- which(!(pdat$gene %in% highlight.genes))
    segments(x0=pdat$start[highlight.idxs], x1=pdat$end[highlight.idxs], 
             y0=pdat$y[highlight.idxs], y1=pdat$y[highlight.idxs], 
             lwd=4, lend="butt", col=highlight.col)
  }else{
    other.idxs <- 1:nrow(pdat)
  }
  segments(x0=pdat$start[other.idxs], x1=pdat$end[other.idxs], 
           y0=pdat$y[other.idxs], y1=pdat$y[other.idxs], 
           lwd=4, lend="butt", col=col)
  
  # Label any genes, if optioned
  if(!is.na(label.genes)){
    sapply(label.genes, function(gene){
      best.pip <- pdat$PIP[which(pdat$gene==gene)]
      best.y <- pdat$y[which(pdat$gene==gene)]
      if(best.pip > 0.6 * max.pip){
        best.pos <- 1
        best.y <- best.y + 0.1*panel.height
      }else{
        best.pos <- 3
        best.y <- best.y - 0.1*panel.height
      }
      text(x=mean(as.numeric(pdat[which(pdat$gene==gene), c("start", "end")])), 
           y=best.y, labels=gene, font=3, cex=5/6, pos=best.pos)
    })
  }
}

# Add mirrored case & control CNV pileups for a single cohort to an existing coordinate plot
add.cnv.panel <- function(cnvs, n.case, n.ctrl, y0, cnv.type, highlight.hpo=NULL,
                          max.freq=NULL, start=NULL, end=NULL, 
                          y.axis.title="CNV\nFreq.", expand.pheno.label=TRUE,
                          case.legend.side="left", case.legend.topbottom="top", 
                          ctrl.legend.side="left", cc.legend.colors=rep(blueblack, 2),
                          add.cohort.label=FALSE, cohort.label=NULL,
                          panel.height=2, dx=100){
  # Standardize inputs
  n.case <- as.numeric(n.case)
  n.ctrl <- as.numeric(n.ctrl)
  
  # Set horizontal dimensions
  if(is.null(start)){
    start <- par("usr")[1]
  }
  if(is.null(end)){
    end <- par("usr")[2]
  }
  
  # Set CNV-based plotting values
  if(cnv.type=="DEL"){
    col.case.highlight <- cnv.blacks[1]
    col.case.other <- cnv.colors[1]
    col.ctrl <- control.cnv.colors[1]
    col.midline <- cnv.blacks[1]
  }else if(cnv.type=="DUP"){
    col.case.highlight <- cnv.blacks[2]
    col.case.other <- cnv.colors[2]
    col.ctrl <- control.cnv.colors[2]
    col.midline <- cnv.blacks[2]
  }
  
  # Collect raw (unscaled) CNV pileups
  raw.pileups <- lapply(cnvs, pileup.cnvs, start=start, end=end, dx=dx, cnv.buffer=0)
  max.case.n <- max(unlist(lapply(raw.pileups$case$cnvs, function(l){max(l$y)})))
  max.ctrl.n <- max(unlist(lapply(raw.pileups$ctrl$cnvs, function(l){max(l$y)})))
  
  # Set scaling
  half.height <- panel.height / 2 
  ytop <- y0 + half.height
  ybottom <- y0 - half.height
  if(is.null(max.freq)){
    max.case.freq <- (max.case.n + 1) / n.case
    max.ctrl.freq <- (max.ctrl.n + 1) / n.ctrl
    max.freq <- max(c(max.case.freq, max.ctrl.freq))
  }
  case.cnv.height <- half.height / ceiling(n.case * max.freq)
  ctrl.cnv.height <- half.height / ceiling(n.ctrl * max.freq)
  
  # Gather scaled CNV pileups
  case.pileup <- pileup.cnvs(cnvs$case, start=start, end=end, dx=dx, cnv.height=case.cnv.height,
                             col=col.case.other, highlight.hpo=highlight.hpo, highlight.col=col.case.highlight)
  ctrl.pileup <- pileup.cnvs(cnvs$ctrl, start=start, end=end, dx=dx, cnv.height=ctrl.cnv.height,
                             col=col.ctrl)
  
  # Add horizontal gridlines
  y.ax.tick.spacing <- seq(-half.height, half.height, length.out=7)
  abline(h=c(y0, y0 + y.ax.tick.spacing, y0 - y.ax.tick.spacing), col="white")
  
  # Plot midline, pileups, and outlines
  lapply(case.pileup$cnvs, function(l){polygon(l$x, y0 + l$y, border=NA, col=l$color)})
  points(case.pileup$counts[, 1], case.pileup$counts[, 2] + y0, type="l", col=col.case.other)
  lapply(ctrl.pileup$cnvs, function(l){polygon(l$x, y0 - l$y, border=NA, col=l$color)})
  points(ctrl.pileup$counts[, 1], -ctrl.pileup$counts[, 2] + y0, type="l", col=col.ctrl)
  abline(h=y0, col=col.midline)
  
  # Add Y axes
  y.ax.label.cex <- 5/6
  axis(2, at=y0 + y.ax.tick.spacing, labels=NA, col=blueblack, tck=-0.0075)
  max.freq.fmt <- format.scientific(max.freq, nsmall=0, max.decimal=0)
  axis(2, at=y0-half.height, tick=F, las=2, line=-0.65, labels=max.freq.fmt, 
       cex.axis=y.ax.label.cex)
  axis(2, at=y0, tick=F, las=2, line=-0.65, labels=0,
       cex.axis=y.ax.label.cex)
  axis(2, at=y0+half.height, tick=F, las=2, line=-0.65, labels=max.freq.fmt,
       cex.axis=y.ax.label.cex)
  axis(2, at=y0, line=0.25, tick=F, labels=y.axis.title)
  
  # Add case legend
  legend.y.cex <- 0.7
  legend.text.cex <- 5/6
  if(!is.na(case.legend.side)){
    case.legend.y <- y0 + (legend.y.cex * half.height)
    if(case.legend.side == "left"){
      case.legend.x <- start
      case.legend.pos <- 4
    }else{
      case.legend.x <- end
      case.legend.pos <- 2
    }
    if(expand.pheno.label == TRUE){
      case.legend.label <- paste(hpo.abbrevs[case.hpo], " cases (N=", prettyNum(n.case, big.mark=","), ")", sep="")
    }else{
      case.legend.label <- paste("Cases (N=", prettyNum(n.case, big.mark=","), ")", sep="")
    }
    
    text(x=case.legend.x, y=case.legend.y, labels=case.legend.label, 
         pos=case.legend.pos, col=cc.legend.colors[1], cex=legend.text.cex)
  }
  
  # Add control legend
  if(!is.na(ctrl.legend.side)){
    ctrl.legend.y <- y0 - (legend.y.cex * half.height)
    if(ctrl.legend.side == "left"){
      ctrl.legend.x <- start
      ctrl.legend.pos <- 4
    }else{
      ctrl.legend.x <- end
      ctrl.legend.pos <- 2
    }
    ctrl.legend.label <- paste("Controls (N=", prettyNum(n.ctrl, big.mark=","), ")", sep="")
    text(x=ctrl.legend.x, y=ctrl.legend.y, labels=ctrl.legend.label, 
         pos=ctrl.legend.pos, col=cc.legend.colors[2], cex=legend.text.cex)
  }
  
  # Add cohort label (upper panel, opposite case.legend.side)
  if(add.cohort.label==TRUE){
    if(case.legend.topbottom=="top"){
      cohort.label.y.adj <- 1
    }else{
      cohort.label.y.adj <- -1
    }
    cohort.label.y <- y0 + (cohort.label.y.adj * legend.y.cex * half.height)
    if(case.legend.side=="left"){
      cohort.label.x <- end
      cohort.label.pos <- 2
    }else{
      cohort.label.x <- start
      cohort.label.pos <- 4
    }
    text(x=cohort.label.x, y=cohort.label.y, labels=cohort.label, 
         pos=cohort.label.pos, col=blueblack, cex=legend.text.cex)
  }
  
  # Add cleanup top & bottom lines
  abline(h=c(ytop, ybottom), col=blueblack)
  segments(x0=par("usr")[2], x1=par("usr")[2], y0=ybottom, y1=ytop, col=blueblack, xpd=T)
}

# Plot CNV key
add.cnv.key <- function(cnv.type, y0, total.n.ctrls, all.case.hpos, total.n.cases,
                        highlight.case.hpo=NA, total.n.cases.highlight=NA,
                        panel.height=0.2, text.cex=5/6, pt.cex=1.3){
  # Get plot data
  x.at <- par("usr")[1] + (c(0.075, 0.65) * diff(par("usr")[1:2]))
  y.buf <- panel.height/3
  
  # Set CNV-based plotting values
  if(cnv.type=="DEL"){
    col.case.highlight <- cnv.blacks[1]
    col.case.other <- cnv.colors[1]
    col.ctrl <- control.cnv.colors[1]
    cnv.label <- "Deletions:"
  }else if(cnv.type=="DUP"){
    col.case.highlight <- cnv.blacks[2]
    col.case.other <- cnv.colors[2]
    col.ctrl <- control.cnv.colors[2]
    cnv.label <- "Duplications:"
  }
  
  # Allow everything to plot beyond the boundaries
  par(xpd=T)
  
  # Add left-most text
  axis(2, at=y0, tick=F, line=-2, las=2, labels=cnv.label, cex=text.cex)
  
  # Add case labels
  all.case.label <- paste(hpo.abbrevs[all.case.hpos[1]], " cases (total N=", 
                          prettyNum(total.n.cases, big.mark=","), ")", sep="")
  all.case.label <- gsub("cases cases", "cases", all.case.label, fixed=T)
  if(is.na(highlight.case.hpo)){
    points(x=x.at[1], y=y0, pch=15, col=col.case.other, cex=pt.cex)
    text(x.at[1], y=y0, pos=4, xpd=T, cex=text.cex, labels=all.case.label)
  }else{
    points(x=x.at[1], y=y0+y.buf, pch=15, col=col.case.other, cex=pt.cex)
    text(x.at[1], y=y0+y.buf, pos=4, xpd=T, cex=text.cex, labels=all.case.label)
    points(x=x.at[1], y=y0-y.buf, pch=15, col=col.case.highlight, cex=pt.cex)
    text(x.at[1], y=y0-y.buf, pos=4, xpd=T, cex=text.cex, 
         labels=paste(hpo.abbrevs[highlight.case.hpo], " cases (total N=", 
                      prettyNum(total.n.cases.highlight, big.mark=","), ")", sep=""))
  }
  
  # Add control label
  ctrl.label <- paste("Controls (total N=", 
                      prettyNum(total.n.ctrls, big.mark=","), 
                      ")", sep="")
  points(x=x.at[2], y=y0, pch=15, col=col.ctrl, cex=pt.cex)
  text(x.at[2], y=y0, pos=4, xpd=T, cex=text.cex, labels=ctrl.label)
}

