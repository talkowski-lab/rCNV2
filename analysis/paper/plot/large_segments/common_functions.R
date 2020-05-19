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
  shared.to.drop.from.segs <- shared.columns[which(shared.columns != "region_id")]
  segs.sub <- segs[, -which(colnames(segs) %in% shared.to.drop.from.segs)]
  gw <- merge(loci, segs.sub, by="region_id", all.x=T, all.y=F, 
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
             corral="random", corralWidth=width)
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

# Function to plot segment permutation test results
plot.seg.perms <- function(gw, perms, feature, measure, n.bins=100, 
                           x.title=NULL, xlims=NULL,
                           diamond.cex=1.5, parmar=c(2.25, 2, 0.5, 0.5)){
  # Get plot data
  perm.dat <- do.call("rbind", lapply(perms, function(df){
    if(measure == "mean"){
      c("ALL" = mean(df[, which(colnames(df)==feature)], na.rm=T),
        "DEL" = mean(df[which(df$cnv=="DEL"), 
                        which(colnames(df)==feature)], na.rm=T),
        "DUP" = mean(df[which(df$cnv=="DUP"), 
                        which(colnames(df)==feature)], na.rm=T))
    }else if(measure == "median"){
      c("ALL" = median(df[, which(colnames(df)==feature)], na.rm=T),
        "DEL" = median(df[which(df$cnv=="DEL"), 
                        which(colnames(df)==feature)], na.rm=T),
        "DUP" = median(df[which(df$cnv=="DUP"), 
                        which(colnames(df)==feature)], na.rm=T))
    }else if(measure == "sum"){
      c("ALL" = sum(df[, which(colnames(df)==feature)]),
        "DEL" = sum(df[which(df$cnv=="DEL"), 
                       which(colnames(df)==feature)]),
        "DUP" = sum(df[which(df$cnv=="DUP"), 
                       which(colnames(df)==feature)]))
    }else if(measure == "frac.any"){
      c("ALL" = length(df[, which(colnames(df)==feature)] > 0) / nrow(df),
        "DEL" = length(df[which(df$cnv=="DEL"), 
                       which(colnames(df)==feature)] > 0) / length(which(df$cnv=="DEL")),
        "DUP" = length(df[which(df$cnv=="DUP"), 
                       which(colnames(df)==feature)] > 0) / length(which(df$cnv=="DUP")))
    }
  }))
  # Convert boolean gw.dat column back to numeric, if needed
  gw.vals <- gw[, which(colnames(gw)==feature)]
  if(all(sapply(gw.vals, function(x){is.logical(x) | is.na(x)}))){
    gw.vals <- sapply(gw.vals, function(x){if(is.na(x)){NA}else if(x==T){1}else if(x==F){0}else{NA}})
  }
  if(measure == "mean"){
    gw.dat <- c("ALL" = mean(gw.vals, na.rm=T),
                "DEL" = mean(gw.vals[which(gw$cnv=="DEL")], na.rm=T),
                "DUP" = mean(gw.vals[which(gw$cnv=="DUP")], na.rm=T))
  }else if(measure == "median"){
    gw.dat <- c("ALL" = median(gw.vals, na.rm=T),
                "DEL" = median(gw.vals[which(gw$cnv=="DEL")], na.rm=T),
                "DUP" = median(gw.vals[which(gw$cnv=="DUP")], na.rm=T))
  }else if(measure == "sum"){
    gw.dat <- c("ALL" = sum(gw.vals, na.rm=T),
                "DEL" = sum(gw.vals[which(gw$cnv=="DEL")], na.rm=T),
                "DUP" = sum(gw.vals[which(gw$cnv=="DUP")], na.rm=T))
  }else if(measure == "frac.any"){
    gw.dat <- c("ALL" = length(which(gw.vals > 0)) / length(gw.vals),
                "DEL" = length(which(gw.vals[which(gw$cnv=="DEL")] > 0)) / length(which(gw$cnv=="DEL")),
                "DUP" = length(which(gw.vals[which(gw$cnv=="DUP")] > 0)) / length(which(gw$cnv=="DUP")))
  }
  val.range <- range(rbind(perm.dat, gw.dat), na.rm=T)
  val.range <- c(floor(val.range[1]), ceiling(1.25*val.range[2]))
  if(val.range[2] - val.range[1] <= n.bins & 
     length(unique(as.numeric(perm.dat))) <= n.bins){
    bins <- val.range[1]:val.range[2]
  }else{
    bins <- seq(val.range[1], val.range[2], length.out=n.bins)
  }
  bin.width <- bins[2]-bins[1]
  if(is.null(xlims)){
    xlims <- range(bins)
  }
  perm.means <- apply(perm.dat, 2, mean, na.rm=T)
  perm.hists <- as.data.frame(apply(perm.dat, 2, function(vals){hist(vals, breaks=bins, plot=F)$counts}))
  perm.pvals <- sapply(1:3, function(i){calc.perm.p(perm.vals=perm.dat[, i], obs.val=gw.dat[i])})
  
  # Helper function to add a single mirrored violin-histogram hybrid of values to an existing plot
  plot.viohist <- function(values, bins, y.at, width=0.8, 
                           obs.val=NA, obs.color=NA, obs.border=NA, 
                           color=bluewhite, border=blueblack,
                           diamond.cex=4, y.title=NULL){
    values[which(values==0)] <- NA
    values <- values / (max(values, na.rm=T) * 2/width)
    rect(xleft=bins[-length(bins)], xright=bins[-1],
         ybottom=-values + y.at, ytop=values + y.at, 
         col=color, border="white")
    segments(x0=rep(bins[-length(bins)], 2),
             x1=rep(bins[-1], 2),
             y0=c(values + y.at, -values + y.at), 
             y1=c(values + y.at, -values + y.at),
             col=border)
    segments(x0=rep(bins, 2), x1=rep(bins, 2),
             y0=c(y.at, values + y.at, y.at, -values + y.at),
             y1=c(values + y.at, y.at, -values + y.at, y.at),
             col=border)
    segments(x0=obs.val, x1=obs.val, 
             y0=y.at - 0.2, y1=y.at + 0.2, 
             col=obs.color, lwd=3, lend="round")
    points(x=obs.val, y=y.at, pch=23, bg=obs.color, 
           col=obs.border, cex=diamond.cex)
    if(!is.null(y.title)){
      axis(2, at=y.at, line=-0.9, tick=F, las=2, labels=y.title)
      # axis(2, at=c(y.at-(width/2), y.at+(width/2)), tck=0, labels=NA, col=blueblack)
    }
  }
  vio.colors <- rep(bluewhite, 3)
  vio.borders <- rep(blueblack, 3)
  # vio.colors <- c(purplewhite, redwhite, bluewhite)
  # vio.borders <- c(purpleblack, redblack, blueblack)
  row.colors <- c(cnv.colors[c(3, 1:2)])
  row.borders <- rep("black", 3)
  # row.borders <- c(purpleblack, redblack, blueblack)
  row.labels <- c("All\nCNV", "DEL", "DUP")
  stats.cex <- 0.85
  
  # Prep global plot area
  par(bty="n", mar=parmar)
  plot(NA, xlim=xlims, ylim=c(3, 0), type="n",
       xaxt="n", yaxt="n", xlab="", ylab="")
  sapply(1:3, function(i){
    plot.viohist(perm.hists[, i], bins, i-0.5,
                 color=vio.colors[i], border=vio.borders[i],
                 y.title=row.labels[i], diamond.cex=diamond.cex, obs.val=gw.dat[i], 
                 obs.color=row.colors[i], obs.border=row.borders[i])
    segments(x0=perm.means[i], x1=perm.means[i],
             y0=i-0.7, y1=i-0.3, lwd=3, 
             col=vio.borders[i], lend="round")
    text(x=gw.dat[i]-(0.03*(par("usr")[2]-par("usr")[1])), y=i-0.7, pos=4, 
         labels=perm.pvals[2, ][[i]], xpd=T, cex=stats.cex)
    # text(x=gw.dat[i]-(0.03*(par("usr")[2]-par("usr")[1])), y=i-0.7, pos=4, 
    #      labels=paste(prettyNum(round(1.73233, 1), small.interval=1), "fold", sep="-"),
    #      xpd=T, cex=stats.cex)
    print(paste(prettyNum(round(gw.dat[i]/mean(perm.dat[, i], na.rm=T), 1), small.interval=1), "fold", sep="-"))
  })
  axis(1, at=unique(c(0, axTicks(1))), labels=NA, col=blueblack, tck=-0.03)
  sapply(1:length(axTicks(1)), function(i){
    axis(1, at=axTicks(1)[i], labels=prettyNum(axTicks(1)[i], big.mark=","), line=-0.65, tick=F)
  })
  mtext(1, text=x.title, line=1.2)
}
