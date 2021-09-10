#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot summary of CNV filtering for rCNV2 formal analyses


options(stringsAsFactors=F, scipen=1000)
require(rCNV2, quietly=T)


######################
### DATA FUNCTIONS ###
######################
# Load CNV sizes and split by DEL/DUP and case/control
load.sizes <- function(paths.in, control.hpo="HEALTHY_CONTROL"){
  # Load table of cohorts and paths
  paths <- read.table(paths.in, header=F, comment.char="", sep="\t")
  paths <- paths[, c(1, ncol(paths))]
  colnames(paths) <- c("cohort", "path")

  # Process cohorts one at a time
  sizes <- lapply(1:nrow(paths), function(i){
    bed <- read.table(paths$path[i], header=T, sep="\t", comment.char="")
    sizes <- bed$end - bed$start
    x <- lapply(c("DEL", "DUP"), function(cnv){
      list("ctrl" = sizes[intersect(which(bed$cnv==cnv),
                                    grep(control.hpo, bed$pheno, fixed=T))],
           "case" = sizes[intersect(which(bed$cnv==cnv),
                                    grep(control.hpo, bed$pheno, fixed=T, invert=T))])
    })
    names(x) <- c("DEL", "DUP")
    return(x)
  })
  names(sizes) <- paths$cohort
  return(sizes)
}

# Load table of CNV summary metadata
load.cnv.stats <- function(stats.in){
  read.table(stats.in, header=T, sep="\t", comment.char="")
}

# Load list of metacohort members
load.metacohorts <- function(metacohorts.in){
  meta <- read.table(metacohorts.in, sep="\t", comment.char="", header=F)
  meta <- meta[grep("meta", meta[, 1], fixed=T), ]
  mnames <- meta[, 1]
  members <- lapply(strsplit(meta[, 2], split=";"), sort)
  names(members) <- mnames
  return(members)
}

# Calculate Poisson confidence interval
pois.ci <- function(lambda, n, conf=0.95){
  z <- qnorm(1 - ((1-conf)/2))
  ci.d <- z * sqrt(lambda / n)
  c(lambda - ci.d, lambda + ci.d)
}

# Collect flattened df of counts of CNVs per sample per cohort
get.cnv.persample <- function(stats, cohorts, cnv, log=F){
  # Get column indexes
  case_cnv.idx <- which(colnames(stats) == paste(cnv, "per_case", sep="_"))
  ctrl_cnv.idx <- which(colnames(stats) == paste(cnv, "per_ctrl", sep="_"))

  # Iterate over cohorts, gather point estimate, and compute 95% Poisson CI
  res <- do.call("rbind", lapply(cohorts, function(cohort){
    cohort.idx <- which(stats$cohort==cohort)
    n_case <- stats$n_case[cohort.idx]
    n_ctrl <- stats$n_ctrl[cohort.idx]
    case_cnv <- stats[cohort.idx, case_cnv.idx]
    ctrl_cnv <- stats[cohort.idx, ctrl_cnv.idx]
    case_ci <- pois.ci(case_cnv, n_case)
    ctrl_ci <- pois.ci(ctrl_cnv, n_ctrl)
    data.frame("cohort" = rep(cohort, 2),
               "cnv" = rep(cnv, 2),
               "pheno" = c("ctrl", "case"),
               "nsamp" = c(n_ctrl, n_case),
               "mean" = c(ctrl_cnv, case_cnv),
               "ci_lower" = c(ctrl_ci[1], case_ci[1]),
               "ci_upper" = c(ctrl_ci[2], case_ci[2]))
  }))
  rownames(res) <- NULL
  if(log==T){
    log.idx <- which(colnames(res) %in% c("mean", "ci_lower", "ci_upper"))
    res[, log.idx] <- apply(res[, log.idx], 2, log10)
    res[, log.idx] <- apply(res[, log.idx], 2, function(vals){replace(vals, is.infinite(vals), NA)})
  }
  return(res)
}

# Collect points & segments for per-sample plots
get.persample.plot.data <- function(del, dup, cohorts, rows.plotted, y.buffer=0.15){
  points <- do.call("rbind", lapply(1:length(cohorts), function(i){
    cohort <- cohorts[i]
    cdat <- as.numeric(sapply(c("ctrl", "case"), function(pheno){
      as.numeric(sapply(list(del, dup), function(df){
        df$mean[which(df$cohort==cohort & df$pheno==pheno)]
      }))
    }))
    cdat <- data.frame("x" = cdat, "color" = c(control.cnv.colors, cnv.colors))
    cdat <- cdat[which(!is.na(cdat$x)), ]
    y.all <- seq(i-1+rows.plotted+y.buffer, i+rows.plotted-y.buffer, length.out=2+nrow(cdat))
    cdat <- data.frame("x" = as.numeric(cdat$x),
                       "y" = as.numeric(y.all)[-c(1, length(y.all))],
                       "color" = as.character(cdat$color))
    return(cdat)
  }))
  segs <- do.call("rbind", lapply(1:length(cohorts), function(i){
    cohort <- cohorts[i]
    cdat <- do.call("rbind", lapply(c("ctrl", "case"), function(pheno){
      t(sapply(list(del, dup), function(df){
        df[which(df$cohort==cohort & df$pheno==pheno),
           which(colnames(df) %in% c("ci_lower", "ci_upper"))]
      }))
    }))
    cdat <- data.frame("x0" = as.numeric(cdat[, 1]),
                       "x1" = as.numeric(cdat[, 2]),
                       "color" = c(control.cnv.colors, cnv.colors))
    cdat <- cdat[which(!is.na(cdat$x0)), ]
    y.all <- seq(i-1+rows.plotted+y.buffer, i+rows.plotted-y.buffer, length.out=2+nrow(cdat))
    cdat <- data.frame("x0" = as.numeric(cdat$x0),
                       "x1" = as.numeric(cdat$x1),
                       "y0" = as.numeric(y.all)[-c(1, length(y.all))],
                       "y1" = as.numeric(y.all)[-c(1, length(y.all))],
                       "color" = as.character(cdat$color))
    return(cdat)
  }))
  rects <- data.frame("xleft" = par("usr")[1], "xright" = par("usr")[2],
                      "ybottom" = (1:length(cohorts))-1+rows.plotted+y.buffer,
                      "ytop" = (1:length(cohorts))+rows.plotted-y.buffer)
  return(list("points" = points, "segs" = segs, "rects" = rects))
}

# Collect plot data for 2d scatterplot
get.2d.data <- function(stats, metacohorts, log=F){
  res <- as.data.frame(do.call("rbind", lapply(c("ctrl", "case"), function(pheno){
    if(pheno=="ctrl"){
      color.vector <- control.cnv.colors
    }else{
      color.vector <- cnv.colors
    }
    do.call("rbind", lapply(c("DEL", "DUP"), function(cnv){
      persamp.idx <- which(colnames(stats)==paste(cnv, "per", pheno, sep="_"))
      size.idx <- which(colnames(stats)==paste("med", pheno, cnv, "size", sep="_"))
      n.idx <- which(colnames(stats)==paste("n", pheno, sep="_"))
      data.frame("x" = as.numeric(stats[, size.idx]),
                 "y" = as.numeric(stats[, persamp.idx]),
                 "color" = rep(color.vector[which(names(color.vector) == cnv)], nrow(stats)),
                 "cnv" = rep(cnv, nrow(stats)))

    }))
  })))
  res <- res[-which(apply(apply(res, 1, is.na), 2, any)), ]
  if(log==T){
    res[, 1:2] <- log10(res[, 1:2])
  }
  return(res)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot horizontal grouped dotplots of CNVs per sample
plot.persample <- function(stats, metacohorts, xlims=NULL, title=NULL,
                           y.axis="left", blue.bg=TRUE){
  # Collect plot data
  del <- get.cnv.persample(stats, unlist(metacohorts), "DEL", log=T)
  dup <- get.cnv.persample(stats, unlist(metacohorts), "DUP", log=T)
  if(is.null(xlims)){
    xlims <- range(c(del$ci_lower, dup$ci_lower, del$ci_upper, dup$ci_upper), na.rm=T)
  }
  n.meta <- length(metacohorts)
  n.cohorts <- sum(sapply(metacohorts, length))
  if(blue.bg==TRUE){
    sep.col <- "white"
    grid.col <- "white"
    plot.bg <- bluewhite
  }else{
    sep.col <- bluewhite
    grid.col <- NA
    plot.bg <- "white"
  }

  # Prep plot area
  par(mar=c(2.2, 6.8, 1.2, 0.3), bty="n")
  if(y.axis == "right"){
    par(mar=c(2.2, 0.3, 1.2, 6.8))
  }
  plot(x=NA, y=NA, xlim=xlims, ylim=c(n.cohorts + n.meta - 1, 0),
       xaxt="n", yaxt="n", xlab="", ylab="", yaxs="i")
  axis(1, at=log10(logscale.minor), labels=NA, tck=-0.015, col=blueblack, lwd=0.75)
  axis(1, at=log10(logscale.demi), labels=NA, tck=-0.03, col=blueblack, lwd=1.25)
  sapply(logscale.demi, function(x){
    axis(1, at=log10(x), labels=prettyNum(x, big.mark=","), tick=F, line=-0.8)
  })
  mtext(1, line=1.15, text="CNVs per Sample")
  mtext(3, line=0.1, text=title, font=2)
  # axis(2, at=-1, las=2, line=-0.8, font=3, col.axis=blueblack, labels="Source", xpd=T, tick=F)

  # Iterate over metacohorts to plot background shading, points, 95% CIs, and Y labels
  rows.plotted <- 0
  for(m in 1:length(metacohorts)){
    cohorts <- sort(metacohorts[[m]])
    pdat <- get.persample.plot.data(del, dup, cohorts, rows.plotted)
    # rect(xleft=par("usr")[1], xright=par("usr")[2],
    #      ybottom=rows.plotted, ytop=rows.plotted+length(cohorts),
    #      col=sep.col, bty="n", border=NA)
    if(blue.bg==TRUE){
      rect(xleft=pdat$rects$xleft, xright=pdat$rects$xright,
           ybottom=pdat$rects$ybottom, ytop=pdat$rects$ytop,
           bty="n", border=NA, col=plot.bg)
    }else{
      segments(x0=pdat$rects$xleft, x1=pdat$rects$xright,
               y0=pdat$rects$ybottom-0.15, y1=pdat$rects$ybottom-0.15,
               col=sep.col)
    }
    segments(x0=log10(logscale.major), x1=log10(logscale.major),
             y0=rep(rows.plotted, length(logscale.major)),
             y1=rep(rows.plotted+length(cohorts), length(logscale.major)),
             col=grid.col, lend="butt")
    segments(x0=pdat$segs$x0, x1=pdat$segs$x1,
             y0=pdat$segs$y0, y1=pdat$segs$y1,
             col=pdat$segs$color)
    points(x=pdat$points$x, y=pdat$points$y, pch=18, col=pdat$points$color)
    rect(xleft=par("usr")[1], xright=par("usr")[2],
         ybottom=rows.plotted, ytop=rows.plotted+length(cohorts),
         col=NA, bty="o", border=blueblack, xpd=T)
    # Add Y-axis labels & metacohort grouping
    if(y.axis=="right"){
      y.axis.side <- 4
    }else{
      y.axis.side <- 2
    }
    sapply(1:length(cohorts), function(i){
      axis(y.axis.side, at=mean(as.numeric(pdat$rects[i, 3:4]), na.rm=T),
           las=2, tick=F, line=-0.8, labels=cohorts[i], col.axis=blueblack,
           cex.axis=0.85)
    })
    left.mar.grouping.at <- c(rows.plotted-0.025, rows.plotted+length(cohorts)+0.025)
    axis(y.axis.side, at=left.mar.grouping.at,
         labels=NA, tck=0.015, col=blueblack, line=3.2, xpd=T)
    axis(y.axis.side, at=mean(left.mar.grouping.at), tick=F, line=2.4, las=2,
         labels=cohort.abbrevs[which(names(cohort.abbrevs) == names(metacohorts)[m])])
    rows.plotted <- rows.plotted + length(cohorts) + 1
  }
}

# Plot horizontal boxplots for CNV sizes for a single cohort
plot.boxes <- function(cohort.sizes, y, y.buffer=0.05, row.width=1){
  flat.list <- list(cohort.sizes[["DEL"]][["ctrl"]],
                    cohort.sizes[["DUP"]][["ctrl"]],
                    cohort.sizes[["DEL"]][["case"]],
                    cohort.sizes[["DUP"]][["case"]])
  all.colors <- c(control.cnv.colors, cnv.colors)
  keep.idx <- which(sapply(flat.list, length) > 0)
  n.boxes <- length(keep.idx)
  plot.list <- flat.list[keep.idx]
  plot.colors <- all.colors[keep.idx]
  box.centers <- seq(y-(row.width/2)+y.buffer, y+(row.width/2)-y.buffer, length.out=n.boxes+2)[-c(1, n.boxes+2)]
  box.width <- (row.width - (2*y.buffer))/4
  sapply(1:length(plot.list), function(i){
    boxplot(log10(plot.list[[i]]), horizontal=T, at=box.centers[i], boxwex=box.width,
            col=plot.colors[i], add=T, lty=1, staplewex=0, outline=F,
            border=plot.colors[i], xaxt="n", xlab="")
    segments(x0=median(log10(plot.list[[i]])), x1=median(log10(plot.list[[i]])),
             y0=box.centers[i]-(box.width/2), y1=box.centers[i]+(box.width/2),
             lwd=2, lend="butt", col=blueblack)
  })
}

# Plot horizontal boxplots of CNV size per cohort
plot.size <- function(sizes, metacohorts, xlims=NULL, title=NULL, y.axis="left",
                      blue.bg=TRUE){
  if(is.null(xlims)){
    xlims <- range(log10(as.numeric(unlist(unlist(sizes)))), na.rm=T)
  }
  n.meta <- length(metacohorts)
  n.cohorts <- sum(sapply(metacohorts, length))
  if(blue.bg==TRUE){
    sep.col <- "white"
    grid.col <- "white"
    plot.bg <- bluewhite
  }else{
    sep.col <- bluewhite
    grid.col <- NA
    plot.bg <- "white"
  }

  # Prep plot area
  par(mar=c(2.2, 6.8, 1.2, 0.6), bty="n")
  if(y.axis=="right"){
    par(mar=c(2.2, 0.6, 1.2, 6.8))
  }
  plot(x=NA, y=NA, xlim=xlims, ylim=c(n.cohorts + n.meta - 1, 0),
       xaxt="n", yaxt="n", xlab="", ylab="", yaxs="i")
  axis(1, at=log10(logscale.minor), labels=NA, tck=-0.015, col=blueblack, lwd=0.75)
  axis(1, at=log10(logscale.major.bp), labels=NA, tck=-0.03, col=blueblack, lwd=1.25)
  sapply(1:length(logscale.major.bp), function(i){
    axis(1, at=log10(logscale.major.bp[i]), labels=logscale.major.bp.labels[i], tick=F, line=-0.8)
  })
  mtext(1, line=1.15, text="CNV Size")
  mtext(3, line=0.1, text=title, font=2)

  # Iterate over metacohorts to plot background shading, points, 95% CIs, and Y labels
  rows.plotted <- 0
  for(m in 1:length(metacohorts)){
    cohorts <- sort(metacohorts[[m]])
    sdat <- sizes[[m]]
    if(blue.bg==TRUE){
      rect(xleft=par("usr")[1], xright=par("usr")[2],
           ybottom=(1:length(cohorts))-1+rows.plotted+0.15,
           ytop=(1:length(cohorts))+rows.plotted-0.15,
           bty="n", border=NA, col=plot.bg)
    }else{
      segments(x0=par("usr")[1], x1=par("usr")[2],
               y0=(1:length(cohorts))-1+rows.plotted,
               y1=(1:length(cohorts))-1+rows.plotted,
               col=sep.col)
    }
    segments(x0=log10(logscale.major), x1=log10(logscale.major),
             y0=rep(rows.plotted, length(logscale.major)),
             y1=rep(rows.plotted+length(cohorts), length(logscale.major)),
             col=grid.col, lend="butt")
    sapply(1:length(cohorts), function(i){
      plot.boxes(sizes[[which(names(sizes)==cohorts[i])]], y=rows.plotted+i-0.5)
    })
    rect(xleft=par("usr")[1], xright=par("usr")[2],
         ybottom=rows.plotted, ytop=rows.plotted+length(cohorts),
         col=NA, bty="o", border=blueblack, xpd=T)
    # Add Y-axis labels & metacohort grouping
    if(y.axis=="right"){
      y.axis.side <- 4
    }else{
      y.axis.side <- 2
    }
    sapply(1:length(cohorts), function(i){
      axis(y.axis.side, at=i+rows.plotted-0.5, las=2, tick=F, line=-0.8,
           labels=cohorts[i], col.axis=blueblack, cex.axis=0.85)
    })
    left.mar.grouping.at <- c(rows.plotted-0.025, rows.plotted+length(cohorts)+0.025)
    axis(y.axis.side, at=left.mar.grouping.at,
         labels=NA, tck=0.015, col=blueblack, line=3.2, xpd=T)
    axis(y.axis.side, at=mean(left.mar.grouping.at), tick=F, line=2.4, las=2,
         labels=cohort.abbrevs[which(names(cohort.abbrevs) == names(metacohorts)[m])])
    rows.plotted <- rows.plotted + length(cohorts) + 1
  }
}

# Plot 2d scatterplot of CNVs per sample vs median CNV size
plot.2dscatter <- function(stats, metacohorts, persample.lims=NULL, size.lims=NULL,
                           background=TRUE, title=NULL, y.axis="left", blue.bg=TRUE){
  # Collect plot data
  pdat <- get.2d.data(stats, metacohorts, log=T)
  if(is.null(size.lims)){
    xlims <- range(pdat$x, na.rm=T)
  }else{
    xlims <- size.lims
  }
  if(is.null(persample.lims)){
    ylims <- range(pdat$y, na.rm=T)
  }else{
    ylims <- persample.lims
  }
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
  par(mar=c(2.5, 2.5, 1.2, 1.2), bty="n")
  if(y.axis=="right"){
    par(mar=c(2.5, 1.2, 1.2, 2.5))
  }
  plot(x=NA, y=NA, xlim=xlims, ylim=ylims,
       xaxt="n", yaxt="n", xlab="", ylab="")
  mtext(3, line=0.1, text=title, font=2)

  # Add background shading
  if(background==T){
    rect(xleft=par("usr")[1], xright=par("usr")[2],
         ybottom=par("usr")[3], ytop=par("usr")[4],
         bty=plot.bty, border=plot.border, col=plot.bg)
    # abline(h=log10(logscale.minor), v=log10(logscale.minor), lwd=0.75, col="white")
    abline(h=log10(logscale.major), v=log10(logscale.major), lwd=1, col=grid.col)
  }

  # Add X axis
  axis(1, at=log10(logscale.minor), labels=NA, tck=-0.015, col=blueblack, lwd=0.75)
  axis(1, at=log10(logscale.major.bp), labels=NA, tck=-0.03, col=blueblack, lwd=1.25)
  sapply(1:length(logscale.major.bp), function(i){
    axis(1, at=log10(logscale.major.bp[i]), labels=logscale.major.bp.labels[i], tick=F, line=-0.75)
  })
  mtext(1, line=1.15, text="Median CNV Size")

  # Add Y axis
  if(y.axis=="right"){
    y.axis.side <- 4
  }else{
    y.axis.side <- 2
  }
  axis(y.axis.side, at=log10(logscale.minor), labels=NA, tck=-0.015, col=blueblack, lwd=0.75)
  axis(y.axis.side, at=log10(logscale.demi), labels=NA, tck=-0.03, col=blueblack, lwd=1.25)
  sapply(logscale.demi, function(x){
    axis(y.axis.side, at=log10(x), labels=prettyNum(x, big.mark=","), tick=F, line=-0.6, las=2)
  })
  mtext(y.axis.side, line=1.6, text="CNVs per Sample")

  # Add points
  points(x=pdat$x, y=pdat$y, pch=18, col=pdat$color)

  # Cleanup box
  box(bty="o", col=blueblack)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse)

# List of command-line options
option_list <- list(
  make_option(c("--control-hpo"), help="String to consider as a control HPO code.",
              default="HEALTHY_CONTROL")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog raw.input.tsv raw.stats.tsv rcnv.input.tsv rcnv.stats.tsv metacohorts.tsv out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 6){
  stop(paste("Four positional arguments required: raw.input.tsv, raw.stats.tsv, rcnv.input.tsv, rcnv.stats.tsv, metacohorts.tsv, and out_prefix\n", sep=" "))
}

# Writes args & opts to vars
raw.paths.in <- args$args[1]
raw.stats.in <- args$args[2]
rcnv.paths.in <- args$args[3]
rcnv.stats.in <- args$args[4]
metacohorts.in <- args$args[5]
out.prefix <- args$args[6]
control.hpo <- opts$`control-hpo`

# # DEV PARAMETERS
# setwd("~/scratch/")
# raw.paths.in <- "~/scratch/raw_cnv.input.tsv"
# raw.stats.in <- "~/scratch/raw_cnv.stats.tsv"
# rcnv.paths.in <- "~/scratch/rCNV.input.tsv"
# rcnv.stats.in <- "~/scratch/rCNV.stats.tsv"
# metacohorts.in <- "~/scratch/rCNV_metacohort_list.txt"
# out.prefix <- "~/scratch/cnv_filtering_summary"
# control.hpo <- "HEALTHY_CONTROL"

# Set colors
cnv.colors <- cnv.colors[1:2]
control.cnv.colors <- control.cnv.colors[1:2]

# Read data
raw.sizes <- load.sizes(raw.paths.in)
raw <- load.cnv.stats(raw.stats.in)
rcnv.sizes <- load.sizes(rcnv.paths.in)
rcnv <- load.cnv.stats(rcnv.stats.in)
cohorts <- sort(raw$cohort)
metacohorts <- load.metacohorts(metacohorts.in)

# Report size statistics
cat("Raw CNV size quantiles:\n")
print(quantile(unlist(raw.sizes)))
cat("Harmonized CNV size quantiles:\n")
print(quantile(unlist(rcnv.sizes)))

# Plot CNVs per sample (raw & filtered)
persamp.vals <- log10(as.numeric(unlist(lapply(list(raw, rcnv), function(df){df[, grep("_per_", colnames(df), fixed=T)]}))))
persamp.xlims <- range(persamp.vals[which(!is.infinite(persamp.vals))], na.rm=T)
persamp.pdf.dims <- c(3.7, 3.2)
pdf(paste(out.prefix, "cnv_per_sample.raw.pdf", sep="."),
    height=persamp.pdf.dims[1], width=persamp.pdf.dims[2])
plot.persample(raw, metacohorts, xlims=persamp.xlims, title="Raw Data", blue.bg=FALSE)
dev.off()
pdf(paste(out.prefix, "cnv_per_sample.filtered.pdf", sep="."),
    height=persamp.pdf.dims[1], width=persamp.pdf.dims[2])
plot.persample(rcnv, metacohorts, xlims=persamp.xlims, title="Harmonized Data",
               y.axis="right", blue.bg=FALSE)
dev.off()

# Plot CNV sizes
# size.xlims <- quantile(log10(as.numeric(unlist(raw.sizes))), probs=c(0.025, 0.975), na.rm=T)
size.xlims <- log10(c(1000, 1000000))
size.pdf.dims <- c(3.7, 3.1)
pdf(paste(out.prefix, "sizes.raw.pdf", sep="."),
    height=size.pdf.dims[1], width=size.pdf.dims[2])
plot.size(raw.sizes, metacohorts, xlims=size.xlims, title="Raw Data", blue.bg=FALSE)
dev.off()
pdf(paste(out.prefix, "sizes.filtered.pdf", sep="."),
    height=size.pdf.dims[1], width=size.pdf.dims[2])
plot.size(rcnv.sizes, metacohorts, xlims=size.xlims, title="Harmonized Data",
          y.axis="right", blue.bg=FALSE)
dev.off()

# Plot 2d scatter of CNVs per sample vs median CNV size
medsize.xlims <- log10(range(sapply(list(raw, rcnv),
                                 function(df){
                                   range(df[, grep("_size", colnames(df), fixed=T)], na.rm=T)
                                 }),
                          na.rm=T))
scatter.pdf.dim <- 2.5
pdf(paste(out.prefix, "size_vs_per_sample.raw.pdf", sep="."),
    height=scatter.pdf.dim, width=scatter.pdf.dim)
plot.2dscatter(raw, size.lims=medsize.xlims, metacohorts, persample.lims=persamp.xlims,
               title="Raw Data", blue.bg=FALSE)
dev.off()
pdf(paste(out.prefix, "size_vs_per_sample.filtered.pdf", sep="."),
    height=scatter.pdf.dim, width=scatter.pdf.dim)
plot.2dscatter(rcnv, size.lims=medsize.xlims, metacohorts, persample.lims=persamp.xlims,
               title="Harmonized Data", y.axis="right", blue.bg=FALSE)
dev.off()
