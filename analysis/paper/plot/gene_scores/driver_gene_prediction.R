#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Map rCNV gene scores into genomic disorders to prioritize driver genes


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Map gene scores into large rCNV segments
score.segs <- function(segs, scores){
  seg.scores <- lapply(1:nrow(segs), function(i){
    if(segs$n_genes[i] > 0){
      genes <- unlist(segs$genes[i])
      cnvtype <- segs$cnv[i]
      if(cnvtype=="DUP"){
        seg.score.df <- scores[which(scores$gene %in% genes), c("gene", "pTriplo")]
      }else{
        seg.score.df <- scores[which(scores$gene %in% genes), c("gene", "pHaplo")]
      }
      colnames(seg.score.df) <- c("gene", "score")
      seg.score.df[order(-seg.score.df$score), ]
    }else{
      return(data.frame("gene"=character(), "score"=numeric()))
    }
  })
  names(seg.scores) <- segs$region_id
  return(seg.scores)
}

# Make table of predicted drivers for each segments
make.driver.table <- function(segs, seg.scores){
  cols.to.keep <- c("chr", "start", "end", "region_id", "cnv", "n_genes")
  driver.table <- as.data.frame(do.call("rbind", lapply(1:nrow(segs), function(i){
    seg.id <- segs$region_id[i]
    ss.i <- seg.scores[[which(names(seg.scores) == seg.id)]]
    hc.genes <- sort(ss.i$gene[which(ss.i$score >= 0.9)])
    if(length(hc.genes) > 0){
      hc.genes.str <- paste(hc.genes, collapse=",")
    }else{
      hc.genes.str <- "."
    }
    lc.genes <- sort(ss.i$gene[which(ss.i$score < 0.9 & ss.i$score >= 0.5)])
    if(length(lc.genes) > 0){
      lc.genes.str <- paste(lc.genes, collapse=",")
    }else{
      lc.genes.str <- "."
    }
    ns.genes <- sort(ss.i$gene[which(ss.i$score < 0.5)])
    if(length(ns.genes) > 0){
      ns.genes.str <- paste(ns.genes, collapse=",")
    }else{
      ns.genes.str <- "."
    }
    as.character(c(unlist(as.vector(segs[i, cols.to.keep])), hc.genes.str, lc.genes.str, ns.genes.str))
  })))
  colnames(driver.table) <- c(cols.to.keep, "high_conf_drivers", "low_conf_drivers", "other_genes")
  return(driver.table)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot a single horizontal swarmplot of gene scores for a single segment
horiz.swarm <- function(x, at, hc.color, lc.color, ns.color, x.add=0, end.buffer=0.05,
                        hc.cutoff=0.9, lc.cutoff=0.5, height=0.8, pt.cex=1){
  pt.colors <- sapply(x, function(x){
    if(x>=hc.cutoff){
      hc.color
    }else if(x>=lc.cutoff){
      lc.color
    }else{
      ns.color
    }
  })
  pt.pch <- sapply(x, function(x){
    if(x>=lc.cutoff){
      16
    }else{
      20
    }
  })
  beeswarm(x+x.add, horiz=T, at=at, add=T, corral="wrap", corralWidth=height,
           pwpch=pt.pch, lwd=pt.lwd, cex=pt.cex, pwcol=pt.colors)
}

# Format segment label to describe region in X-axis of segment gene-score plot
format.seg.label <- function(segs, idx){
  id.parts <- unlist(strsplit(segs$region_id[idx], split="_"))
  n.parts <- length(id.parts)
  if(id.parts[n.parts] %in% LETTERS){
    cyto <- paste(id.parts[(n.parts-1):n.parts], collapse=" ")
  }else{
    cyto <- id.parts[n.parts]
  }
  chrom <- segs$chr[idx]
  start.fmt <- format(round(segs$start[idx] / 1000000, 1), nsmall=1)
  end.fmt <- format(round(segs$end[idx] / 1000000, 1), nsmall=1)
  coords <- paste("chr", chrom, ":", start.fmt, "-", end.fmt, "Mb", sep="")
  paste(cyto, "\n(", coords, ")", sep="")
}

# Plot scores for all segments from a given CNV type
plot.all.seg.scores <- function(segs, seg.scores, cnv.type, hc.cutoff=0.9, lc.cutoff=0.5,
                                pt.cex=1, n.cols=4, parmar=c(0.2, 0.2, 3, 0.2)){
  # Get plot data
  plot.dat <- seg.scores[which(names(seg.scores) %in% segs$region_id[which(segs$cnv==cnv.type)])]
  if(cnv.type=="DEL"){
    hc.color <- cnv.colors[1]
    lc.color <- control.cnv.colors[1]
    ax.title <- "pHaplo"
  }else{
    hc.color <- cnv.colors[2]
    lc.color <- control.cnv.colors[2]
    ax.title <- "pTriplo"
  }
  n.rows <- ceiling(length(plot.dat) / n.cols)
  x.buffer <- 1
  y.buffer <- 0.1
  x.left <- seq(0, 2*(n.cols-1), 1+x.buffer)
  y.bottom <- (0:(n.rows-1))+y.buffer
  y.top <- (1:n.rows)-y.buffer

  # Prep plot area
  par(bty="n", mar=parmar)
  plot(NA, xlim=c(-x.buffer, max(x.left)+1), ylim=c(n.rows, 0),
       xaxt="n", yaxt="n", xlab="", ylab="", yaxs="i")

  # Add plot areas, gridlines, and shading
  sapply(1:n.rows, function(ridx){
    sapply(1:n.cols, function(cidx){
      seg.idx <- ((ridx-1)*n.cols) + cidx
      if(seg.idx <= length(plot.dat)){
        rect(xleft=x.left[cidx], xright=x.left[cidx]+1,
             ybottom=y.bottom[ridx], ytop=y.top[ridx],
             col=bluewhite, border=NA)
        segments(x0=x.left[cidx]+lc.cutoff, x1=x.left[cidx]+lc.cutoff,
                 y0=y.bottom[ridx], y1=y.top[ridx], col=lc.color)
        segments(x0=x.left[cidx]+hc.cutoff, x1=x.left[cidx]+hc.cutoff,
                 y0=y.bottom[ridx], y1=y.top[ridx], col=hc.color)
      }
    })
  })

  # Add swarmplots & segment annotations
  sapply(1:n.rows, function(ridx){
    sapply(1:n.cols, function(cidx){
      seg.idx <- ((ridx-1)*n.cols) + cidx
      if(seg.idx <= length(plot.dat)){
        seg.name <- names(plot.dat)[seg.idx]
        orig.vals <- plot.dat[[seg.idx]]$score
        horiz.swarm(orig.vals, at=ridx-0.5, height=1-(2*(y.buffer+0.025)),
                    hc.color=hc.color, lc.color=lc.color, ns.color=ns.color,
                    x.add=x.left[cidx], hc.cutoff=hc.cutoff, lc.cutoff=lc.cutoff,
                    pt.cex=pt.cex)
        hc.idxs <- which(orig.vals>=hc.cutoff)
        if(length(hc.idxs) > 0 & length(hc.idxs) < 3){
          text(x=orig.vals[hc.idxs[1]]+x.left[cidx]+0.04, y=ridx-0.5-0.2,
               pos=2, cex=0.85, font=3, col="black",
               labels=plot.dat[[seg.idx]]$gene[hc.idxs[1]])
          if(length(hc.idxs) == 2){
            text(x=orig.vals[hc.idxs[2]]+x.left[cidx]+0.04, y=ridx-0.5+0.2,
                 pos=2, cex=0.85, font=3, col="black",
                 labels=plot.dat[[seg.idx]]$gene[hc.idxs[2]])
          }
        }
        text(x=x.left[cidx], y=ridx-0.5, pos=2, xpd=T,
             labels=format.seg.label(segs, idx=which(segs$region_id==seg.name)))
      }
    })
  })

  # Add cleanup outer boxes
  sapply(1:n.rows, function(ridx){
    sapply(1:n.cols, function(cidx){
      seg.idx <- ((ridx-1)*n.cols) + cidx
      if(seg.idx <= length(plot.dat)){
        rect(xleft=x.left[cidx], xright=x.left[cidx]+1,
             ybottom=y.bottom[ridx], ytop=y.top[ridx],
             col=NA, border=blueblack, xpd=T)
      }
    })
  })

  # Add top X axes
  x.ax.labs <- seq(0, 1, 0.2)
  sapply(1:n.cols, function(cidx){
    x.ax.at <- seq(x.left[cidx], x.left[cidx]+1, 0.2)
    axis(3, at=x.ax.at, labels=NA, col=blueblack)
    sapply(1:length(x.ax.at), function(x){
      axis(3, at=x.ax.at[x], tick=F, line=-0.5, labels=x.ax.labs[x])
    })
    axis(3, at=x.left[cidx]+0.5, tick=F, line=0.5, labels=ax.title)
  })
}


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(optparse, quietly=T)
require(beeswarm, quietly=T)

# List of command-line options
option_list <- list(
  make_option("--score-cutoffs", metavar="tsv", default=NULL,
              help=".tsv of gene score cutoffs generated by empirical_score_cutoffs.R")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog scores.tsv segs.tsv out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 3){
  stop(paste("Three positional arguments required: scores.tsv, segs.tsv, output_prefix\n", sep=" "))
}

# Writes args & opts to vars
scores.in <- args$args[1]
segs.in <- args$args[2]
out.prefix <- args$args[3]
score.cutoffs.in <- opts$`score-cutoffs`

# # DEV PARAMETERS
# scores.in <- "~/scratch/rCNV.gene_scores.tsv.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d2.master_segments.bed.gz"
# out.prefix <- "~/scratch/driver_gene_test"
# score.cutoffs.in <- "~/scratch/rCNV2.gene_score_cutoff.test.tsv"

# Load scores & cutoffs
scores <- load.scores(scores.in)
load.score.cutoffs(score.cutoffs.in)

# Load segment table
segs <- load.segment.table(segs.in)

# Subset to segs at GW/FDR significance or Bonf significance among GDs
segs <- segs[which(segs$any_sig | segs$bonf_sig_gd), ]

# Map scores into segs while matching on CNV type
seg.scores <- score.segs(segs, scores)

# Compute fractions of segments with at least one high-confidence gene
n.segs.w.driver <- length(which(sapply(seg.scores, function(df){any(df$score >= 0.9)})))
cat(paste("Proportion of segments with at least one DS gene:",
          round(n.segs.w.driver / nrow(segs), 3), "\n"))
n.dup.segs.w.driver <- length(which(sapply(seg.scores[which(names(seg.scores) %in% segs$region_id[which(segs$cnv=="DUP")])],
                                           function(df){any(df$score >= 0.9)})))
cat(paste("Proportion of DUP segments with at least one TS gene:",
          round(n.dup.segs.w.driver / length(which(segs$cnv=="DUP")), 3), "\n"))

# Compute fractions of segments with two or fewer high-confidence drivers
n.segs.w.lt2.drivers <- length(which(sapply(seg.scores, function(df){
  length(which(df$score >= 0.9)) %in% c(1, 2)
  })))
cat(paste("Proportion of segments with no more than two DS genes:",
          round(n.segs.w.lt2.drivers / nrow(segs), 3), "\n"))
n.segs.w.one.driver <- length(which(sapply(seg.scores, function(df){
  length(which(df$score >= 0.9)) == 1
})))
cat(paste("Proportion of segments with exactly one DS gene:",
          round(n.segs.w.one.driver / nrow(segs), 3), "\n"))

# Format segments with candidate drivers as table
driver.table <- make.driver.table(segs, seg.scores)
colnames(driver.table)[1] <- paste("#", colnames(driver.table)[1], sep="")
driver.table.outpath <- paste(out.prefix, "gd_driver_genes.table.bed", sep=".")
write.table(driver.table, driver.table.outpath,
            col.names=T, row.names=F, sep="\t", quote=F)
system(paste("bgzip -f", driver.table.outpath), wait=T)

# Plot distributions of scores for all rCNV segments
n.cols <- 3
# Deletions, NAHR
n.segs <- length(which(segs$cnv=="DEL" & segs$nahr==T))
pdf(paste(out.prefix, "gd_driver_genes.scores_per_segment.DEL.NAHR.pdf", sep="."),
    height=2*(((n.segs/n.cols)*(11-1.5)/34) + 1), width=2*(8.5-1.5))
plot.all.seg.scores(segs[which(segs$nahr==T), ], seg.scores, "DEL", n.cols=n.cols,
                    hc.cutoff=hi.hc, lc.cutoff=hi.lc)
dev.off()
# Deletions, nonrecurrent
n.segs <- length(which(segs$cnv=="DEL" & segs$nahr==F))
pdf(paste(out.prefix, "gd_driver_genes.scores_per_segment.DEL.nonrecurrent.pdf", sep="."),
    height=2*(((n.segs/n.cols)*(11-1.5)/34) + 1), width=2*(8.5-1.5))
plot.all.seg.scores(segs[which(segs$nahr==F), ], seg.scores, "DEL", n.cols=n.cols,
                    hc.cutoff=hi.hc, lc.cutoff=hi.lc)
dev.off()
# Duplications, NAHR
n.segs <- length(which(segs$cnv=="DUP" & segs$nahr==T))
pdf(paste(out.prefix, "gd_driver_genes.scores_per_segment.DUP.NAHR.pdf", sep="."),
    height=2*(((n.segs/n.cols)*(11-1.5)/34) + 1), width=2*(8.5-1.5))
plot.all.seg.scores(segs[which(segs$nahr==T), ], seg.scores, "DUP", n.cols=n.cols,
                    hc.cutoff=ts.hc, lc.cutoff=ts.lc)
dev.off()
# Duplications, nonrecurrent
n.segs <- length(which(segs$cnv=="DUP" & segs$nahr==F))
pdf(paste(out.prefix, "gd_driver_genes.scores_per_segment.DUP.nonrecurrent.pdf", sep="."),
    height=2*(((n.segs/n.cols)*(11-1.5)/34) + 1), width=2*(8.5-1.5))
plot.all.seg.scores(segs[which(segs$nahr==F), ], seg.scores, "DUP", n.cols=n.cols,
                    hc.cutoff=ts.hc, lc.cutoff=ts.lc)
dev.off()
