#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot summary of HPO hierarchy & sample sizes


options(stringsAsFactors=F, scipen=1000)
require(rCNV2, quietly=T)


######################
### DATA FUNCTIONS ###
######################
# Load sample overlap matrix
load.ovr <- function(ovr.in, hpos){
  ovr <- read.table(ovr.in, sep="\t", comment.char="", header=T, check.names=F)
  rownames(ovr) <- ovr$HPO
  ovr <- ovr[, -1]
  ovr <- ovr[match(hpos, rownames(ovr)), match(hpos, colnames(ovr))]

  # Insert null rows separating neuro & somatic phenotypes
  neuro.idx <- which(hpos %in% neuro.hpos)
  soma.idx <- which(hpos %in% somatic.hpos)
  ovr.wnulls <- as.data.frame(rbind(ovr[1, ],
                                    NA,
                                    ovr[neuro.idx, ],
                                    NA,
                                    ovr[soma.idx, ],
                                    NA,
                                    ovr[(max(soma.idx) + 1):nrow(ovr), ]))
  return(ovr.wnulls)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot descriptive text columns
plot.text.columns <- function(hpos, meta, x.at=c(0, 2, 8), y.lab.buffer=0.05,
                              background=T, y.buffer=0.1, y.top=-0.6){
  # Get basic plotting info
  nhpos <- nrow(meta)
  y.at <- (1:nhpos)-0.5
  hpo.labels <- meta$HPO
  hpo.labels[which(hpo.labels == "UNKNOWN")] <- "N/A"
  abbrevs <- hpo.abbrevs[match(meta$HPO, names(hpo.abbrevs), nomatch="")]

  # Prep plot area & add y-axis titles
  plot(x=NA, y=NA, xlim=c(0, 10), ylim=c(nhpos, y.top),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  sapply(1:length(x.at), function(i){
    axis(3, at=c(x.at[i]+0.125, (c(x.at[-1], par("usr")[2]) - 0.25)[i]),
         tck=0, labels=NA, xpd=T, col=blueblack, line=0.2+y.top)
  })
  height <- par("usr")[4]-par("usr")[3]
  text(x=x.at, y=par("usr")[4]+(y.lab.buffer*height)+(y.top*(sum(par("mar")[c(1,3)])/height)),
       labels=c("HPO", "Description", "Abbreviation"),
       xpd=T, font=2, pos=4)

  # Plot text
  if(background==T){
    rect(xleft=par("usr")[1], xright=par("usr")[2],
         ybottom=y.at[which(!is.na(meta$HPO))]-0.5+y.buffer,
         ytop=y.at[which(!is.na(meta$HPO))]+0.5-y.buffer,
         border=NA, bty="n", col=bluewhite)
  }
  text(x=sapply(x.at, function(x){rep(x, nhpos)}), y=rep(y.at, 3),
       labels=c(hpo.labels, meta$description, abbrevs), pos=4, xpd=T)
}

# Plot bars of samples per HPO, colored by cohort
plot.count.bars <- function(hpos, meta, scaled=F, title=NULL,
                            x.scalar=1.01, y.buffer=0.1, end.labels=F,
                            y.lab.buffer=0.05, background=T, x.axis=T,
                            y.top=-0.6){
  # Get basic plotting info
  counts <- meta[, grep("meta", colnames(meta), fixed=T)]
  counts$total <- apply(counts, 1, sum, na.rm=T)
  counts$total[which(apply(counts, 1, function(vals){any(is.na(vals))}))] <- NA
  if(scaled==TRUE){
    counts <- as.data.frame(t(apply(counts, 1, function(vals){vals/vals[length(vals)]})))
  }
  cumsums <- cbind(rep(0, nrow(counts)), t(apply(counts[, -ncol(counts)], 1, cumsum)))
  x.max <- x.scalar * max(counts, na.rm=T)
  nhpos <- nrow(meta)
  y.at <- (1:nhpos)-0.5
  ybottoms <- sapply(y.at, function(y){rep(y-0.5+y.buffer, ncol(counts)-1)})
  ytops <- sapply(y.at, function(y){rep(y+0.5-y.buffer, ncol(counts)-1)})

  # Prep plot area
  plot(x=NA, y=NA, xlim=c(0, x.max), ylim=c(nhpos, y.top),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  if(background==T){
    rect(xleft=par("usr")[1], xright=par("usr")[2],
         ybottom=y.at[which(!is.na(meta$HPO))]-0.5+y.buffer,
         ytop=y.at[which(!is.na(meta$HPO))]+0.5-y.buffer,
         border=NA, bty="n", col=bluewhite)
  }

  # Add rectangles
  rect(xleft=as.vector(t(cumsums[, -ncol(cumsums)])),
       xright=as.vector(t(cumsums[, -1])),
       ybottom=ybottoms, ytop=ytops,
       col=cohort.colors, border=blueblack, lwd=0.15)
  rect(xleft=rep(0, nrow(counts)), xright=counts$total,
       ybottom=y.at-0.5+y.buffer, ytop=y.at+0.5-y.buffer,
       col=NA, border=blueblack, xpd=T)

  # Add title, axis, and labels, if optioned
  if(x.axis==T){
    if(scaled==T){
      ax.at <- seq(0, 1, 0.25)
      # ax.lab <- paste(seq(0, 100, 25), "%", sep="")
      ax.lab <- ax.at
    }else{
      ax.at <- seq(0, round(x.max, -4), length.out=5)
      ax.lab <- paste(ax.at / 1000, "k", sep="")
      ax.lab[1] <- 0
    }
    axis(3, at=ax.at, labels=NA, line=0.2+y.top, col=blueblack, tck=-0.035)
    sapply(1:length(ax.at), function(x){
      axis(3, at=ax.at[x], tick=F, line=-0.4+y.top, labels=ax.lab[x], cex=0.8)
    })
    title.line <- 0.025
  }else{
    title.line <- 0
  }
  height <- par("usr")[4]-par("usr")[3]
  text(x=mean(par("usr")[1:2]), y=par("usr")[4]+((y.lab.buffer+title.line)*height)+(y.top*(sum(par("mar")[c(1,3)])/height)),
       labels=title, xpd=T, font=2)
  if(end.labels==T){
    text(x=counts$total, y=y.at, pos=4, xpd=T,
         labels=prettyNum(counts$total, big.mark=","))
  }
}

# Heatmap of sample overlap
plot.overlaps <- function(hpos, ovr, max.ovr=0.2, title.y.scalar=0.2, y.top=-0.6){
  # Get basic plotting info
  plot.ovr <- as.data.frame(t(apply(ovr, 1, function(vals){
    vals[which(vals>=max.ovr)] <- max.ovr
    vals / max.ovr
  })))
  nhpos <- ncol(ovr)
  nhpos.y <- nrow(ovr)

  # Prep plot area
  plot(x=NA, y=NA, xlim=c(0, nhpos), ylim=c(nhpos.y, y.top),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  axis(3, at=(1:nhpos)-0.5, tick=F, line=-0.5, las=2,
       labels=hpo.abbrevs[match(colnames(ovr), names(hpo.abbrevs), nomatch="")])

  # Add top X-axis legend rectangles
  sapply(1:nhpos, function(i){
    # Above diagonal colored by neuro
    # Below diagonal colored by severity
    hpo.color.top <- get.hpo.color(colnames(ovr)[i], color.by="neuro")
    hpo.color.bottom <- get.hpo.color(colnames(ovr)[i], color.by="severity")
      polygon(x=i-0.5+c(-0.4, -0.4, 0.4, -0.4),
              y=y.top+0.35+c(-0.7, 0, -0.7, -0.7),
              xpd=T, border=blueblack, col=hpo.color.top)
      polygon(x=i-0.5+c(-0.4, 0.4, 0.4, -0.4),
              y=y.top+0.35+c(0, -0.7, 0, 0),
              xpd=T, border=blueblack, col=hpo.color.bottom)
  })

  # Plot rectangles
  sapply(1:nhpos.y, function(row){
    row.vals <- plot.ovr[row, ]
    row.colors <- percentile.palette[as.numeric(round(100*row.vals, 0) + 1)]
    rect(xleft=(1:nhpos)-1, xright=1:nhpos,
         ybottom=row, ytop=row-1,
         col=row.colors, border=row.colors, xpd=T)
  })
  abline(v=1:nhpos, h=1:nhpos.y, col="white")

  # Plot outside borders
  na.rows <- which(!(rownames(ovr) %in% hpos))
  sapply(1:length(na.rows), function(i){
    rect(xleft=0, xright=nhpos, ybottom=na.rows[i]-1, ytop=na.rows[i],
         col="white", border="white", bty="n", xpd=T)
  })
  sapply(0:length(na.rows), function(i){
    rect(xleft=0, xright=nhpos, ybottom=c(0, na.rows)[i+1], ytop=c(na.rows, nhpos.y+1)[i+1]-1,
         col=NA, xpd=T, border=blueblack)
  })

  # Add title
  text(x=mean(par("usr")[1:2]), y=par("usr")[4]+(title.y.scalar*(par("usr")[4]-par("usr")[3])),
       labels="Fraction of Cases Overlapping Other Phenotypes", xpd=T, font=2)
  segments(x0=0.15, x1=nhpos-0.15,
           y0=par("usr")[4]+((title.y.scalar-0.012)*(par("usr")[4]-par("usr")[3])),
           y1=par("usr")[4]+((title.y.scalar-0.012)*(par("usr")[4]-par("usr")[3])),
           xpd=T, col=blueblack)
  segments(x0=c(0.15, nhpos-0.15), x1=c(0.15, nhpos-0.15),
           y0=rep(par("usr")[4]+((title.y.scalar-0.012)*(par("usr")[4]-par("usr")[3])), 2),
           y1=rep(par("usr")[4]+((title.y.scalar-(1.5*0.012))*(par("usr")[4]-par("usr")[3])), 2),
           xpd=T, col=blueblack)
}

# Master figure plot wrapper
plot.hpos <- function(hpos, meta, ovr, panel.widths=c(2, 8, 3, 4, 10),
                      parmar=c(0.3, 0.1, 5, 0.1), max.ovr=0.2,
                      y.lab.buffer.text=0.01, y.lab.buffer.bars=0.01){
  # Prep layout
  layout(matrix(c(1:5), nrow=1, byrow=T), widths=panel.widths)
  par(bty="n", mar=parmar)

  # Plot dendrogram in left margin
  plot.dendro(meta, color=blueblack, lwd=2)

  # Add text labels
  plot.text.columns(hpos, meta, x.at=c(0, 1.5, 7.5), y.lab.buffer.text)

  # Plot bars of pct per cohort and total count
  parmar.leftbars <- parmar
  parmar.leftbars[4] <- parmar[4] + 0.25
  par(mar=parmar.leftbars)
  plot.count.bars(hpos, meta, scaled=T, title="Proportion of Cases", y.buffer=0.05,
                  y.lab.buffer=y.lab.buffer.bars)
  parmar.rightbars <- parmar
  parmar.rightbars[2] <- parmar[2] + 0.25
  par(mar=parmar.rightbars)
  plot.count.bars(hpos, meta, scaled=F, title="Total Cases",
                  x.scalar=1.5, y.buffer=0.05, end.labels=T,
                  y.lab.buffer=y.lab.buffer.bars)

  # Heatmap of sample overlaps
  par(mar=parmar)
  plot.overlaps(hpos, ovr, max.ovr=max.ovr, title.y.scalar=0.18)
}

# Plot metacohort legend
plot.cohort.legend <- function(){
  par(mar=c(0.1, 0.1, 0.1, 0.1), bty="n")
  plot(x=NA, y=NA, xlim=c(0, 3), ylim=c(5, 0),
       xaxt="n", yaxt="n", xlab="", ylab="")
  x.at <- rep(1, 4)
  y.at <- 1:4
  rect(xleft=x.at-0.3, xright=x.at+0.3,
       ybottom=y.at-0.3, ytop=y.at+0.3,
       col=cohort.colors, border=blueblack)
  text(x=x.at+0.2, y=y.at, pos=4, labels=cohort.abbrevs)
}

# Plot phenotype legend
plot.pheno.legend <- function(){
  par(mar=c(0.1, 0.1, 0.1, 0.1), bty="n")
  plot(x=NA, y=NA, xlim=c(0, 3), ylim=c(4, 0),
       xaxt="n", yaxt="n", xlab="", ylab="")
  x.at <- rep(1, 3)
  y.at <- 1:3
  rect(xleft=x.at-0.3, xright=x.at+0.3,
       ybottom=y.at-0.3, ytop=y.at+0.3,
       col=pheno.colors, border=blueblack)
  text(x=x.at+0.2, y=y.at, pos=4, labels=pheno.abbrevs)
}

# Plot overlap legend
plot.overlap.legend <- function(max.ovr, title.line=3){
  par(mar=c(0.1, 0.5, 2, 0.5), bty="n")
  plot(x=NA, y=NA, xlim=c(0, 101), ylim=c(0, 1),
       xaxt="n", yaxt="n", xlab="", ylab="", yaxs="i")
  rect(xleft=0:100, xright=1:101, ybottom=0, ytop=1,
       col=percentile.palette, border=percentile.palette)
  rect(xleft=0, xright=101, ybottom=0, ytop=1, col=NA, border=blueblack, xpd=T)
  x.at <- seq(0, 100, 25)
  x.labels <- seq(0, max.ovr, length.out=length(x.at))
  x.labels[length(x.labels)] <- paste(">", x.labels[length(x.labels)], sep="")
  # axis(3, at=x.at+0.5, labels=NA, col=blueblack, line=0)
  sapply(1:length(x.at), function(x){
    axis(3, at=(x.at+0.5)[x], tick=F, labels=x.labels[x], line=-0.9)
  })
  text(x=mean(par("usr")[1:2]), y=title.line, xpd=T, labels="Fraction Overlapping")
}


#####################
### RSCRIPT BLOCK ###
#####################
require(optparse)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog hpos.txt metadata.tsv counts.tsv overlap.tsv out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 5){
  stop(paste("Five positional arguments required: hpos.txt, metadata.tsv, counts.tsv, overlap.tsv, out_prefix\n", sep=" "))
}

# Writes args & opts to vars
hpos.in <- args$args[1]
meta.in <- args$args[2]
counts.in <- args$args[3]
ovr.in <- args$args[4]
out.prefix <- args$args[5]

# # DEV PARAMETERS
# hpos.in <- "~/scratch/rCNV2_analysis_d2.reordered_hpos.txt"
# meta.in <- "~/scratch/phenotype_groups.HPO_metadata.txt"
# counts.in <- "~/scratch/HPOs_by_metacohort.table.tsv"
# ovr.in <- "~/scratch/rCNV2_analysis_d2.hpo_sample_overlap_fraction_matrix.tsv"
# out.prefix <- "~/scratch/test"

# Read data
hpos <- as.character(read.table(hpos.in, header=F, sep="\t")[, 1])
meta <- load.hpo.metadata(meta.in, counts.in, hpos)
ovr <- load.ovr(ovr.in, hpos)

# Plot figure
max.ovr <- 0.2
vertical.scale <- 0.75
panel.widths <- c(0.8, 6.75, 1.6, 2, 8)
pdf(paste(out.prefix, "pdf", sep="."),
    height=vertical.scale*11,
    width=vertical.scale*sum(panel.widths))
plot.hpos(hpos, meta, ovr, panel.widths,
          parmar=c(0.2, 0.1, 10, 0.1),
          y.lab.buffer.text=0.0075,
          y.lab.buffer.bars=0.0075,
          max.ovr=max.ovr)
dev.off()

# Plot legends
# pdf(paste(out.prefix, "cohort_legend.pdf", sep="."),
#     height=1.75, width=1.5)
# plot.cohort.legend()
# dev.off()
pdf(paste(out.prefix, "overlap_legend.pdf", sep="."),
    height=0.6, width=2.1)
plot.overlap.legend(max.ovr, title.line=2.7)
dev.off()
# pdf(paste(out.prefix, "pheno_legend.pdf", sep="."),
#     height=1.75, width=1.7)
# plot.pheno.legend()
# dev.off()

