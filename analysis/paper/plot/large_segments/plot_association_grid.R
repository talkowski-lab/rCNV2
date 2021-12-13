#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot grid summarizing large segment association across phenotypes for rCNV2 manuscript


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Annotate HPOs associated with each segment split by significance level
annotate.sig.hpos <- function(segs, assocs.in){
  # Load association data
  assocs <- read.table(assocs.in, check.names=F, sep="\t", header=T, comment.char="")
  assocs <- assocs[, c("credible_set_id", "hpo", "sig_level")]

  # Annotate GW-sig HPOs
  segs$gw_sig_hpos <- sapply(segs$constituent_assocs, function(csids){
    unique(assocs$hpo[which(assocs$credible_set_id %in% csids
                            & assocs$sig_level == "genome_wide")])
  })

  # Annotate FDR-sig HPOs
  segs$fdr_sig_hpos <- sapply(segs$constituent_assocs, function(csids){
    unique(assocs$hpo[which(assocs$credible_set_id %in% csids
                            & assocs$sig_level == "FDR")])
  })

  # Return annotated dataframe
  return(segs)
}

# Load flat table of effect sizes and P-values per region per phenotype per CNV
load.sumstats <- function(sumstats.in){
  sumstats <- read.table(sumstats.in, header=T, sep="\t")
  sumstats[, 4:8] <- apply(sumstats[, 4:8], 2, as.numeric)
  return(sumstats)
}

# Load list of clusters to plot on a single line
load.clusters <- function(clusters.in, loci){
  if(!is.null(clusters.in)){
    clusters <- strsplit(read.table(clusters.in, header=F, sep="\t")[, 1], split=",")
  }else{
    clusters <- lapply(loci$region_id, function(x){x})
  }

  # Restrict clusters to regions found in loci
  clusters <- lapply(clusters, function(rids){rids[which(rids %in% loci$region_id)]})
  clusters <- clusters[which(lapply(clusters, length) > 0)]

  # Rename clusters
  names(clusters) <- unlist(lapply(clusters, function(x){
    sort(sapply(strsplit(x, split="_"), function(l){
      l <- l[grep('^[A-Z]$', l, invert=T)]
      l[length(l)]
      }))[1]
    }))
  return(clusters)
}

# Load HPO sample sizes
load.hpo.n <- function(hpo.n.in){
  hdf <- read.table(hpo.n.in, header=T, sep="\t", comment.char="", check.names=F)
  hpo.n <- as.numeric(hdf$Total)
  names(hpo.n) <- as.character(hdf[, 1])
  return(hpo.n)
}

# Average & reformat size for a cluster of loci
format.size <- function(segs, cluster){
  size <- mean(segs$size[which(segs$region_id %in% cluster)])
  paste(prettyNum(round(size / 1000, 0), big.mark=","), "kb")
}

# Summarize & format genes for a cluster of loci
format.genes <- function(segs, cluster, column.name="genes", max.genes=6, max.characters=20, newline.at=3){
  genes <- unique(sort(unlist(segs[which(segs$region_id %in% cluster), which(colnames(segs)==column.name)])))
  ngenes <- length(genes)
  if(ngenes == 0){
    out.fmt <- "No genes"
  }else if(ngenes > max.genes){
    out.fmt <- paste(prettyNum(ngenes, big.mark=","), "genes")
  }else{
    if(ngenes > newline.at){
      genes[newline.at] <- paste(genes[newline.at], ",\n", sep="")
      genes[-c(newline.at, length(genes))] <- sapply(genes[-c(newline.at, length(genes))], function(g){paste(g, ", ", sep="")})
      out.fmt <- paste(genes, collapse="")
    }else{
      out.fmt <- paste(genes, collapse=", ")
    }
    if(nchar(out.fmt) > max.characters){
      out.fmt <- paste(prettyNum(ngenes, big.mark=","), "genes")
    }
  }
  return(out.fmt)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Prep plotting area
prep.plot <- function(sumstats, clusters, hpos,
                      ncols.prestats, prestat.widths, prestat.colnames,
                      ncols.poststats, poststat.widths, poststat.colnames,
                      shading.color=NA, shading.buffer=0.15, top.bracket.tck=0.005,
                      hpo.bracket.line=9, fdr=FALSE, parmar=c(0.5, 0.5, 11, 0.5)){
  # Get plot values
  ncols.hpo <- length(hpos)
  nrows <- length(clusters)
  if(fdr==TRUE){
    left.bracket.title <- "Large segments significant at FDR < 0.01"
  }else{
    left.bracket.title <- "Genome-wide significant large segments"
  }

  # Prep plotting area
  par(bty="n", mar=parmar)
  plot(NA, xlim=c(-sum(prestat.widths), ncols.hpo + sum(poststat.widths)), ylim=c(nrows, -0.5),
       xlab="", xaxt="n",  xaxs="i", ylab="", yaxt="n", yaxs="i")

  # Add background shading
  rect(xleft=(1:ncols.hpo)[which(!is.na(hpos))]-1+shading.buffer,
       xright=(1:ncols.hpo)[which(!is.na(hpos))]-shading.buffer,
       ybottom=nrows-shading.buffer, ytop=par("usr")[4],
       bty="n", border=NA, col=shading.color)
  rect(xleft=par("usr")[1], xright=par("usr")[2],
       ybottom=(1:nrows)-1+shading.buffer,
       ytop=(1:nrows)-shading.buffer,
       bty="n", border=NA, col=shading.color)
  sapply((1:ncols.hpo)[which(!is.na(hpos))], function(i){
    polygon(x=c(i-1+shading.buffer, i-1+shading.buffer,
                i-shading.buffer, i-1+shading.buffer),
            y=c(-0.25, -0.825, -0.825, -0.25), xpd=T,
            border=blueblack, col=get.hpo.color(hpos[i], color.by="neuro"))
    polygon(x=c(i-1+shading.buffer, i-shading.buffer,
                i-shading.buffer, i-1+shading.buffer),
            y=c(-0.25, -0.825, -0.25, -0.25), xpd=T,
            border=blueblack, col=get.hpo.color(hpos[i], color.by="severity"))
  })
  abline(v=c(-cumsum(prestat.widths), 0,
             ncols.hpo+c(0, cumsum(poststat.widths))),
         lwd=3, col="white")

  # Add axes
  prestat.hr.at <- rev(c(0, -cumsum(prestat.widths)))
  sapply(1:ncols.prestats, function(i){
    axis(3, at=c(prestat.hr.at[i]+0.5, prestat.hr.at[i+1]-0.5),
         tck=0, labels=NA, col=blueblack, line=-0.25)
    axis(3, at=mean(c(prestat.hr.at[i]+0.5, prestat.hr.at[i+1]-0.5)),
         tick=F, labels=prestat.colnames[i], line=-1.15, font=2)
  })
  axis(3, at=c(prestat.hr.at[1], prestat.hr.at[ncols.prestats+1]),
       tck=top.bracket.tck, labels=NA, col=blueblack, line=2.35, lend="round")
  axis(3, at=mean(c(prestat.hr.at[1], prestat.hr.at[ncols.prestats+1])),
       tick=F, labels=left.bracket.title, line=1.45)
  axis(3, at=(1:ncols.hpo)[which(!is.na(hpos))]-0.5, tick=F, line=-0.475, las=2,
       labels=hpo.abbrevs[hpos[which(!is.na(hpos))]])
  hpo.bracket.idxs <- list(which(hpos %in% neuro.hpos & hpos %in% adult.hpos),
                           which(hpos %in% neuro.hpos & hpos %in% developmental.hpos),
                           which(hpos %in% somatic.hpos & hpos %in% developmental.hpos),
                           which(hpos %in% somatic.hpos & hpos %in% adult.hpos))
  sapply(1:length(hpo.bracket.idxs), function(k){
    axis(3, at=range(hpo.bracket.idxs[[k]])-c(1, 0), tck=top.bracket.tck,
         labels=NA, col=blueblack, line=hpo.bracket.line, lend="round")
    axis(3, at=mean(range(hpo.bracket.idxs[[k]])-c(1, 0)),
         tick=F, line=hpo.bracket.line-0.9,
         labels=c("Weaker-effect\nneurological", "Strong-effect\nneurological",
                  "Strong-effect\nnon-neuro.", "Weaker-effect\nnon-neuro")[k])
  })
  poststat.hr.at <- c(ncols.hpo, ncols.hpo+cumsum(poststat.widths))
  sapply(1:ncols.poststats, function(i){
    axis(3, at=c(poststat.hr.at[i]+0.5, poststat.hr.at[i+1]-0.5),
         tck=0, labels=NA, col=blueblack, line=-0.25)
    axis(3, at=mean(c(poststat.hr.at[i]+0.5, poststat.hr.at[i+1]-0.5)),
         tick=F, labels=poststat.colnames[i], line=-1.15, font=2)
  })
}

# Plot a single semicircle
plot.semicircle <- function(x, y, r, side="top", color="black", border=NA, lwd=1, xpd=F){
  # Adopted from https://stackoverflow.com/questions/31538534/plotting-half-circles-in-r
  rs <- seq(0, pi, len=100)
  xm <- r*cos(rs)
  ym <- r*sin(rs)
  if(side == "bottom"){
    sign <- -1
  }else{
    sign <- 1
  }
  xc <- x + (sign * xm)
  yc <- y + (sign * ym)
  polygon(xc, yc, col=color, border=border, lwd=lwd, xpd=xpd)
}

# Plot data for a single locus-phenotype pair
plot.locus.single_hpo <- function(segs, sumstats, region_ids, hpo, x, y, max.radius=0.5,
                                  stat.size="pvalue", max.stat.size=8,
                                  stat.color="lnor", max.stat.color=4, fdr=FALSE){
  # Scale radii
  stat.size.idx <- which(colnames(sumstats) == stat.size)
  stat.color.idx <- which(colnames(sumstats) == stat.color)
  radii <- as.numeric(sapply(c("DUP", "DEL"), function(cnv){
    radius.x <- max(sumstats[which(sumstats$region_id %in% region_ids & sumstats$hpo==hpo & sumstats$cnv==cnv), stat.size.idx], na.rm=T)
    max.radius * (min(c(max(c(0, radius.x), na.rm=T), max.stat.size)) / max.stat.size)
  }))

  # Determine color based on stat.color
  colors <- sapply(c("DUP", "DEL"), function(cnv){
    color.x <- max(sumstats[which(sumstats$region_id %in% region_ids & sumstats$hpo==hpo & sumstats$cnv==cnv), stat.color.idx], na.rm=T)
    cnv.color.palettes[[cnv]][round(min(c(max(c(0, color.x), na.rm=T), max.stat.color)) * 100 / max.stat.color) + 1]
  })

  # Determine border based on significance
  sig <- sapply(c("DUP", "DEL"), function(cnv){
    if(fdr){
      sig.hpos <- unlist(segs$fdr_sig_hpos[which(segs$region_id %in% region_ids
                                                & segs$cnv==cnv)])
    }else{
      sig.hpos <- unlist(segs$gw_sig_hpos[which(segs$region_id %in% region_ids
                                                & segs$cnv==cnv)])
    }
    if(any(hpo %in% sig.hpos)){
      cnv.blacks[which(names(cnv.blacks) == cnv)]
    }else{
      NA
    }
  })

  # Plot semicircle shading
  plot.semicircle(x, y, r=radii[1], side="bottom", color=colors[1], border=NA, lwd=2)
  plot.semicircle(x, y, r=radii[2], side="top", color=colors[2], border=NA, lwd=2)

  # Plot semicircle borders
  plot.semicircle(x, y, r=radii[1], side="bottom", color=NA, border=sig[1], lwd=1.5)
  plot.semicircle(x, y, r=radii[2], side="top", color=NA, border=sig[2], lwd=1.5)
}

# Plot a row of sumstats for all hpos per region
plot.locus.row.stats <- function(loci, sumstats, hpos, region_ids, y,
                                 stat.size="pvalue", max.stat.size=8,
                                 stat.color="lnor", max.stat.color=4,
                                 outline=TRUE, outline.color="gray80",
                                 buffer=0.15, fdr=FALSE){
  x.at <- (1:length(hpos))[which(!is.na(hpos))]-0.5
  rmax <- 0.5-(buffer/2)
  if(outline==TRUE){
    sapply(x.at, draw.circle, y, radius=rmax, border=outline.color, col="white")
  }
  sapply((1:length(hpos))[which(!is.na(hpos))], function(x){
    if(!is.na(hpos[x])){
      plot.locus.single_hpo(loci, sumstats, region_ids, hpos[x], (1:length(hpos))[x]-0.5,
                            y, max.radius=rmax, stat.size, max.stat.size, stat.color,
                            max.stat.color, fdr=fdr)
    }
  })
}

# Plot a grid of all regions
plot.all.loci <- function(segs, clusters, sumstats, hpos,
                          stat.size="pvalue", max.stat.size=8,
                          stat.color="lnor", max.stat.color=4,
                          outline.color="gray60", background.color=bluewhite,
                          cex.table.text=0.9, buffer=0.1, fdr=FALSE,
                          parmar=c(0.25, 0.25, 11, 0.25)){
  # Prep plot area
  prestat.widths <- rev(c(8, 3, 4, 3, 3))
  poststat.widths <- c(4, 11, 11)
  if(fdr==TRUE){
    sig.col.title <- "FDR<0.01\nrCNVs"
  }else{
    sig.col.title <- "G-W Sig.\nrCNVs"
  }
  prestat.colnames <- c("Cytoband", sig.col.title, "Size", "NAHR", "Known\nGDs")
  poststat.colnames <- c("# Genes", "Genes Constrained\nAgainst Truncating SNVs",
                         "Phenotype-Matched\nDisease Genes")
  prep.plot(sumstats, clusters, hpos,
            ncols.prestats=5, prestat.widths=prestat.widths,
            prestat.colnames=prestat.colnames,
            ncols.poststats=3, poststat.widths=poststat.widths,
            poststat.colnames=poststat.colnames,
            shading.color=background.color, shading.buffer=buffer,
            fdr=fdr, parmar=parmar)

  # Get plot data
  nhpos <- length(hpos)
  nrows <- length(clusters)
  y.at <- (1:nrows) - 0.5
  x.at.left <- c(rev(cumsum(-prestat.widths)), 0)
  x.at.right <- nhpos + c(0, cumsum(poststat.widths))

  # Add locus names
  text(x=mean(par("usr")[1]), y=y.at, pos=4,
       labels=names(clusters), cex=cex.table.text)

  # Add significant CNV association
  sapply(1:nrows, function(i){
    cnv <- unique(segs$cnv[which(segs$region_id %in% clusters[[i]])])
    if(length(cnv) == 2){
      text(x=mean(x.at.left[2:3]), y=y.at[i], cex=cex.table.text,
           labels="Recip.", col=cnv.colors[3])
    }else{
      text(x=mean(x.at.left[2:3]), y=y.at[i], cex=cex.table.text,
           labels=cnv, col=cnv.colors[which(names(cnv.colors) == cnv)])
    }
  })

  # Add locus sizes
  text(x=x.at.left[4], y=y.at, pos=2, cex=cex.table.text,
       labels=sapply(clusters, format.size, segs=segs))

  # Add NAHR label
  lapply(1:nrows, function(i){
    rids <- clusters[[i]]
    if(any(segs$nahr[which(segs$region_id %in% rids)])){
      text(x=mean(x.at.left[4:5]), y=y.at[i], cex=cex.table.text,
           labels="Yes")
    }else{
      text(x=mean(x.at.left[4:5]), y=y.at[i], cex=cex.table.text,
           labels="No", col=control.cnv.colors[2])
    }
  })

  # Add GD overlap
  gd.hits <- sapply(clusters, function(rids){
    chrom <- unique(segs$chr[which(segs$region_id %in% rids)])
    start <- min(segs$start[which(segs$region_id %in% rids)])
    end <- max(segs$end[which(segs$region_id %in% rids)])
    get.gd.overlap(chrom, start, end, segs)
  })
  sapply(1:nrows, function(i){
    hits <- gd.hits[[i]]
    nhits <- length(hits)
    if(nhits == 0){
      text(x=mean(x.at.left[5:6]), y=y.at[i], cex=cex.table.text,
           labels="-", col=control.cnv.colors[2])
    }else if(nhits == 1){
      text(x=mean(x.at.left[5:6]), y=y.at[i], cex=cex.table.text,
           labels=hits, col=cnv.colors[which(names(cnv.colors) == hits)])
    }else{
      text(x=mean(x.at.left[5:6]), y=y.at[i], cex=cex.table.text,
           labels="Recip.", col=cnv.colors[3])
    }
  })

  # Plot association sumstats each region
  sapply(1:length(clusters), function(ridx){
    region_ids <- clusters[[ridx]]
    plot.locus.row.stats(segs, sumstats, hpos, region_ids, ridx-0.5,
                         stat.size, max.stat.size, stat.color, max.stat.color,
                         outline=T, outline.color, buffer, fdr=fdr)
  })

  # Add column for number of genes
  sapply(1:length(clusters), function(i){
    n.genes <- length(unique(unlist(segs[which(segs$region_id %in% clusters[[i]]), "genes"])))
    if(n.genes == 0){
      text.color <- control.cnv.colors[2]
    }else{
      text.color <- "black"
    }
    text(x=mean(x.at.right[1:2]), y=y.at[i], font=1, col=text.color,
         cex=0.95*cex.table.text, labels=n.genes)
  })

  # Add text for constrained genes and HPO-matched genes
  glists <- c("gnomAD_constrained_genes", "HPOmatched_genes")
  sapply(1:length(glists), function(gi){
    genes <- sapply(clusters, format.genes, segs=segs,
                    column.name=glists[gi],
                    max.genes=5, newline.at=100, max.characters=25)
    sapply(1:length(clusters), function(i){
      if(genes[[i]] == "No genes"){
        font <- 1
        text.color <- control.cnv.colors[2]
      }else{
        text.color <- "black"
        if(length(grep(" genes", genes[[i]], fixed=T)) == 0){
          font <- 3
        }else{
          font <- 1
        }
      }
      text(x=mean(x.at.right[gi + 1:2]), y=y.at[i], font=font, col=text.color,
           cex=0.95*cex.table.text, labels=genes[[i]])
    })
  })
}

# Plot P-value legend for locus grid
plot.p.legend <- function(max.p=6, max.is.gw=TRUE){
  pvals <- unique(round(seq(0, max.p, length.out=4)))[-1]
  # Prep plot area
  par(mar=c(rep(0.1, 3), 4), bty="n")
  plot(NA, xlim=c(0, 1), ylim=c(0, 3), asp=1,
       xaxt="n", xlab="", xaxs="i", yaxt="n", ylab="", yaxs="i")
  sapply(1:length(pvals), function(i){
    if(i<length(pvals)){
      sapply(c("top", "bottom"), function(side){
        plot.semicircle(x=0.5, y=i-0.5, r=0.5*((pvals[i])/max.p), side=side,
                        color=control.cnv.colors[2], border=control.cnv.colors[2], xpd=T)
      })
      axis(4, at=i-0.5, las=2, line=-0.8, tick=F,
           labels=bquote(italic(P) == 10^-.(pvals[i])))
    }else{
      if(max.is.gw==TRUE){
        border <- blueblack
      }else{
        border <- control.cnv.colors[2]
      }
      sapply(c("top", "bottom"), function(side){
        plot.semicircle(x=0.5, y=i-0.5, r=0.5*((pvals[i])/max.p), side=side,
                        color=control.cnv.colors[2], border=border, lwd=2, xpd=T)
      })
      axis(4, at=i-0.5, las=2, line=-0.8, tick=F,
           labels=bquote(italic(P) <= 10^-.(pvals[i])))
    }
  })
}

# Plot odds ratio legend for locus grid
plot.lnor.legend <- function(max.lnor=log(8)){
  # Prep plot area
  par(mar=c(0.1, 2, 1.5, 0.8), bty="n")
  plot(NA, xlim=c(0, 101), ylim=c(0, 2),
       xaxt="n", xlab="", xaxs="i", yaxt="n", ylab="", yaxs="i")

  # Add gradients
  sapply(1:2, function(i){
    rect(xleft=0:100, xright=1:101, ybottom=i-1, ytop=i,
         col=cnv.color.palettes[[i]], border=cnv.color.palettes[[i]])
    rect(xleft=0, xright=101, ybottom=i-1, ytop=i,
         col=NA, border=blueblack, xpd=T)
  })

  # Add labels
  axis(2, at=(1:2)-0.5, tick=F, las=2, line=-0.9, labels=c("DEL", "DUP"))
  top.at <- log(c(1, 2, 4, 8)) * (100/max.lnor)
  axis(3, at=top.at, tck=-0.15, labels=NA, col=blueblack)
  axis(3, at=top.at[1], tick=F, line=-0.75, labels=expression(""<=1))
  axis(3, at=top.at[2:3], tick=F, line=-0.75, labels=c(2, 4))
  axis(3, at=top.at[4], tick=F, line=-0.75, labels=expression("">=8))
}


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(optparse, quietly=T)
require(plotrix, quietly=T)
require(funr, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--gw-clusters"), help="Comma-delimited list of GW-sig regions to collapse."),
  make_option(c("--fdr-clusters"), help="Comma-delimited list of FDR-sig regions to collapse."),
  make_option(c("--hpo-jaccard"), help="Tsv matrix of Jaccard indexes for all HPO pairs"),
  make_option(c("--hpo-sample-sizes"), help="Tsv of HPO information with sample sizes")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog loci.bed associations.bed segs.tsv sumstats.tsv hpos.tsv out.prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 6){
  stop(paste("Six positional arguments required: loci.bed, associations.bed, segs.tsv, sumstats.tsv, hpos.tsv, out.prefix\n", sep=" "))
}

# Writes args & opts to vars
loci.in <- args$args[1]
assocs.in <- args$args[2]
segs.in <- args$args[3]
sumstats.in <- args$args[4]
hpos.in <- args$args[5]
out.prefix <- args$args[6]
gw.clusters.in <- opts$`gw-clusters`
fdr.clusters.in <- opts$`fdr-clusters`
hpo.jac.in <- opts$`hpo-jaccard`
hpo.n.in <- opts$`hpo-sample-sizes`

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# assocs.in <- "~/scratch/rCNV.final_segments.associations.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d2.master_segments.bed.gz"
# sumstats.in <- "~/scratch/rCNV2_analysis_d2.all_segs.all_sumstats.tsv.gz"
# hpos.in <- "~/scratch/rCNV2_analysis_d2.reordered_hpos.txt"
# gw.clusters.in <- "~/scratch/locus_clusters.genome_wide.txt"
# fdr.clusters.in <- "~/scratch/locus_clusters.FDR.txt"
# hpo.jac.in <- "~/scratch/rCNV2_analysis_d2.hpo_jaccard_matrix.tsv"
# hpo.n.in <- "~/scratch/HPOs_by_metacohort.table.tsv"
# out.prefix <- "~/scratch/test"

# Load loci & segment table
loci <- load.loci(loci.in)
segs.raw <- load.segment.table(segs.in)
segs <- merge.loci.segs(loci, segs.raw)

# Annotate each segment with separate lists of GW-sig and FDR-sig HPOs
segs <- annotate.sig.hpos(segs, assocs.in)

# Load other data
sumstats <- load.sumstats(sumstats.in)
hpos <- as.character(read.table(hpos.in, header=F, sep="\t")[, 1])
gw.clusters <- load.clusters(gw.clusters.in, segs[which(segs$gw_sig), ])
fdr.clusters <- load.clusters(fdr.clusters.in, segs[which(segs$fdr_sig), ])

# Prune phenotypes based on Jaccard indexes for HPO pairs, if optioned
if(!is.null(hpo.jac.in)){
  hpo.jac <- load.hpo.jaccard.matrix(hpo.jac.in)
  hpo.n <- load.hpo.n(hpo.n.in)
  filtered.hpo.jac <- prune.hpo.jaccard.matrix(hpo.jac, hpo.n, max.jac=0.25)
  hpos <- hpos[which(hpos %in% colnames(filtered.hpo.jac))]
}


# Modify HPOs to have breaks after phenotype clusters
# Ordering: all cases, adult neuro, NDD, non-neuro developmental, adult non-neuro, unknown
hpos <- c(NA, hpos[which(hpos == "HP:0000118")], NA,
  hpos[which(hpos %in% adult.hpos & hpos %in% neuro.hpos)], NA,
  hpos[which(hpos %in% developmental.hpos & hpos %in% neuro.hpos)], NA,
  hpos[which(hpos %in% developmental.hpos & hpos %in% somatic.hpos)], NA,
  hpos[which(hpos %in% adult.hpos & hpos %in% somatic.hpos)], NA,
  hpos[which(hpos == "UNKNOWN")], NA)

# Plot locus grid for genome-wide significant segments
pdf(paste(out.prefix, "large_segments.association_grid.GW_sig.pdf", sep="."),
    height=17, width=16)
plot.all.loci(segs, gw.clusters, sumstats, hpos,
              stat.size="pvalue", max.stat.size=6,
              stat.color="lnor", max.stat.color=log(8),
              outline.color=control.cnv.colors[2], background.color=bluewhite,
              buffer=0.15, cex.table.text=0.85)
dev.off()

# Plot locus grid for FDR-significant segments
pdf(paste(out.prefix, "large_segments.association_grid.FDR_sig.pdf", sep="."),
    height=18, width=16)
plot.all.loci(segs, fdr.clusters, sumstats, hpos,
              stat.size="pvalue", max.stat.size=5,
              stat.color="lnor", max.stat.color=log(8),
              outline.color=control.cnv.colors[2], background.color=bluewhite,
              buffer=0.15, cex.table.text=0.85,
              fdr=TRUE)
dev.off()

# Plot legends
pdf(paste(out.prefix, "pvalue_legend.gw_sig.pdf", sep="."),
    height=1, width=1.04)
plot.p.legend(max.p=6)
dev.off()
pdf(paste(out.prefix, "pvalue_legend.fdr.pdf", sep="."),
    height=1, width=1.04)
plot.p.legend(max.p=5, max.is.gw=FALSE)
dev.off()
pdf(paste(out.prefix, "lnor_legend.pdf", sep="."),
    height=0.7, width=1.4)
plot.lnor.legend()
dev.off()
