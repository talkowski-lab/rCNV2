#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot grid summarizing large segment association across phenotypes for rCNV2 manuscript


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
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
  names(clusters) <- unlist(lapply(clusters, function(x){sort(sapply(strsplit(x, split="_"), function(l){l[length(l)]}))[1]}))
  return(clusters)
}

# Average & reformat size for a cluster of loci
format.size <- function(gw, cluster){
  size <- mean(gw$size[which(gw$region_id %in% cluster)])
  paste(prettyNum(round(size / 1000, 0), big.mark=","), "kb")
}

# Summarize & format genes for a cluster of loci
format.genes <- function(gw, cluster, column.name="genes", max.genes=6, max.characters=20, newline.at=3){
  genes <- unique(sort(unlist(gw[which(gw$region_id %in% cluster), which(colnames(gw)==column.name)])))
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
                      shading.color=NA, shading.buffer=0.15, nominal=FALSE,
                      parmar=c(0.5, 0.5, 8, 0.5)){
  # Get plot values
  ncols.hpo <- length(hpos)
  nrows <- length(clusters)
  if(nominal==TRUE){
    left.bracket.title <- "Nominally significant large segments from literature"
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
  rect(xleft=(1:ncols.hpo)[which(!is.na(hpos))]-1+shading.buffer,
       xright=(1:ncols.hpo)[which(!is.na(hpos))]-shading.buffer,
       ybottom=-0.25, ytop=-0.5, border=blueblack,
       col=sapply(hpos[which(!is.na(hpos))], get.hpo.color), xpd=T)
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
       tck=0.01, labels=NA, col=blueblack, line=2.35, lend="round")
  axis(3, at=mean(c(prestat.hr.at[1], prestat.hr.at[ncols.prestats+1])),
       tick=F, labels=left.bracket.title, line=1.45)
  
  # sapply(1:ncols.hpo, function(i){
  #   axis(3, at=c(i-1+shading.buffer, i-shading.buffer),
  #        tck=0, labels=NA, col=blueblack)
  # })
  axis(3, at=(1:ncols.hpo)[which(!is.na(hpos))]-0.5, tick=F, line=-0.9, las=2,
       labels=hpo.abbrevs[hpos[which(!is.na(hpos))]])
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
plot.locus.single_hpo <- function(loci, sumstats, region_ids, hpo, x, y, max.radius=0.5,
                                  stat.size="pvalue", max.stat.size=8,
                                  stat.color="lnor", max.stat.color=4){
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
    if(any(hpo %in% unlist(loci$hpos[which(loci$region_id %in% region_ids & loci$cnv==cnv)]))){
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
                                 outline=TRUE, outline.color="gray80", buffer=0.15){
  x.at <- (1:length(hpos))[which(!is.na(hpos))]-0.5
  rmax <- 0.5-(buffer/2)
  if(outline==TRUE){
    sapply(x.at, draw.circle, y, radius=rmax, border=outline.color, col="white")
  }
  sapply((1:length(hpos))[which(!is.na(hpos))], function(x){
    if(!is.na(hpos[x])){
      plot.locus.single_hpo(loci, sumstats, region_ids, hpos[x], (1:length(hpos))[x]-0.5, y, max.radius=rmax,
                            stat.size, max.stat.size, stat.color, max.stat.color)
    }
  })
}

# Plot a grid of all regions
plot.all.loci <- function(gw, clusters, sumstats, hpos, 
                          stat.size="pvalue", max.stat.size=8, 
                          stat.color="lnor", max.stat.color=4,
                          outline.color="gray60", background.color=bluewhite,
                          cex.table.text=0.9, buffer=0.1, nominal=FALSE, 
                          parmar=c(0.25, 0.25, 8, 0.25)){
  # Prep plot area
  prestat.widths <- rev(c(7, 3, 4, 3, 3))
  poststat.widths <- rep(10, 3)
  if(nominal==TRUE){
    prestat.colnames <- c("Cytoband", "Nominal\nrCNVs", "Size", "NAHR", "Known\nGDs")
    poststat.colnames <- c("All Genes", "Genes Constrained\nAgainst Truncating SNVs", "All OMIM\nDisease Genes")
  }else{
    prestat.colnames <- c("Cytoband", "Signif.\nrCNVs", "Size", "NAHR", "Known\nGDs")
    poststat.colnames <- c("All Genes", "Genes Constrained\nAgainst Truncating SNVs", "Phenotype-Matched\nDisease Genes")
  }
  prep.plot(sumstats, clusters, hpos,
            ncols.prestats=5, prestat.widths=prestat.widths, 
            prestat.colnames=prestat.colnames,
            ncols.poststats=3, poststat.widths=poststat.widths,
            poststat.colnames=poststat.colnames,
            shading.color=background.color, shading.buffer=buffer,
            nominal=nominal, parmar=parmar)
  
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
    cnv <- unique(gw$cnv[which(gw$region_id %in% clusters[[i]])])
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
       labels=sapply(clusters, format.size, gw=gw))
  
  # Add NAHR label
  lapply(1:nrows, function(i){
    rids <- clusters[[i]]
    if(any(gw$nahr[which(gw$region_id %in% rids)])){
      text(x=mean(x.at.left[4:5]), y=y.at[i], cex=cex.table.text,
           labels="Yes")
    }else{
      text(x=mean(x.at.left[4:5]), y=y.at[i], cex=cex.table.text,
           labels="No", col=control.cnv.colors[2])
    }
  })
  
  # Add GD overlap
  gd.hits <- sapply(clusters, function(rids){
    chrom <- unique(gw$chr[which(gw$region_id %in% rids)])
    start <- min(gw$start[which(gw$region_id %in% rids)])
    end <- max(gw$end[which(gw$region_id %in% rids)])
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
    plot.locus.row.stats(loci, sumstats, hpos, region_ids, ridx-0.5,
                         stat.size, max.stat.size, stat.color, max.stat.color,
                         outline=T, outline.color, buffer)
  })
  
  # Add genes, constrained genes, and HPO-matched genes
  if(nominal==TRUE){
    glists <- c("genes", "gnomAD_constrained_genes", "OMIM_genes")
  }else{
    glists <- c("genes", "gnomAD_constrained_genes", "HPOmatched_genes")
  }
  sapply(1:length(glists), function(gi){
    genes <- sapply(clusters, format.genes, gw=gw, 
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
      text(x=mean(x.at.right[gi:(gi+1)]), y=y.at[i], font=font, col=text.color,
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
require(optparse, quietly=T)
require(plotrix, quietly=T)
require(funr, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--gw-clusters"), help="Comma-delimited list of gw-sig regions to collapse."),
  make_option(c("--lit-clusters"), help="Comma-delimited list of literature regions to collapse."),
  make_option(c("--rcnv-config"), help="rCNV2 config file to be sourced.")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog loci.bed segs.tsv sumstats.tsv hpos.tsv out.prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 5){
  stop(paste("Five positional arguments required: loci.bed, segs.tsv, sumstats.tsv, hpos.tsv, out.prefix\n", sep=" "))
}

# Writes args & opts to vars
loci.in <- args$args[1]
segs.in <- args$args[2]
sumstats.in <- args$args[3]
hpos.in <- args$args[4]
out.prefix <- args$args[5]
gw.clusters.in <- opts$`gw-clusters`
lit.clusters.in <- opts$`lit-clusters`
rcnv.config <- opts$`rcnv-config`

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d1.master_segments.bed.gz"
# sumstats.in <- "~/scratch/rCNV2_analysis_d1.all_segs.all_sumstats.tsv.gz"
# hpos.in <- "~/scratch/rCNV2_analysis_d1.reordered_hpos.txt"
# gw.clusters.in <- "~/scratch/locus_clusters.gw_sig.txt"
# lit.clusters.in <- "~/scratch/locus_clusters.all_gds.txt"
# rcnv.config <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/config/rCNV2_rscript_config.R"
# out.prefix <- "~/scratch/test"
# script.dir <- "~/Desktop/Collins/Talkowski/CNV_DB/rCNV_map/rCNV2/analysis/paper/plot/large_segments/"

# Source rCNV2 config, if optioned
if(!is.null(rcnv.config)){
  source(rcnv.config)
}

# Source common functions
script.dir <- funr::get_script_path()
source(paste(script.dir, "common_functions.R", sep="/"))

# Load loci & segment table
loci <- load.loci(loci.in)
segs <- load.segment.table(segs.in)

# Merge loci & segment data for genome-wide significant sites only
gw <- merge.loci.segs(loci, segs)

# Subset segment data for nominal (but not genome-wide significant) sites
nom <- segs[which(segs$nom_sig & !segs$gw_sig), ]

# Print size distributions for various subsets of segments:
cat("Size stats for genome-wide significant segments:")
quantile(gw$size, probs=seq(0, 1, 0.25))

# Load other data
sumstats <- load.sumstats(sumstats.in)
hpos <- as.character(read.table(hpos.in, header=F, sep="\t")[, 1])
gw.clusters <- load.clusters(gw.clusters.in, loci)
lit.clusters <- load.clusters(lit.clusters.in, nom)

# Modify HPOs to have breaks after phenotype clusters
hpos <- c(NA, hpos[1], NA,
          hpos[2:max(which(hpos %in% neuro.hpos))], NA,
          hpos[min(which(hpos %in% somatic.hpos)):max(which(hpos %in% somatic.hpos))], NA,
          hpos[length(hpos)], NA)

# Plot locus grid for genome-wide significant segments
pdf(paste(out.prefix, "large_segments.association_grid.gw_sig.pdf", sep="."), 
    height=10.25, width=16)
plot.all.loci(gw, gw.clusters, sumstats, hpos, 
              stat.size="pvalue", max.stat.size=6, 
              stat.color="lnor", max.stat.color=log(8),
              outline.color=control.cnv.colors[2], background.color=bluewhite, 
              buffer=0.15, cex.table.text=0.85)
dev.off()

# Plot locus grid for nominal (not genome-wide significant) segments
pdf(paste(out.prefix, "large_segments.association_grid.nominal.pdf", sep="."), 
    height=11, width=16)
plot.all.loci(nom, lit.clusters, sumstats, hpos, 
              stat.size="pvalue", max.stat.size=5, 
              stat.color="lnor", max.stat.color=log(8),
              outline.color=control.cnv.colors[2], background.color=bluewhite, 
              buffer=0.15, cex.table.text=0.85,
              nominal=TRUE)
dev.off()

# Plot legends
pdf(paste(out.prefix, "pvalue_legend.gw_sig.pdf", sep="."),
    height=1, width=1.04)
plot.p.legend(max.p=6)
dev.off()
pdf(paste(out.prefix, "pvalue_legend.nominal.pdf", sep="."),
    height=1, width=1.04)
plot.p.legend(max.p=5, max.is.gw=FALSE)
dev.off()
pdf(paste(out.prefix, "lnor_legend.pdf", sep="."),
    height=0.7, width=1.4)
plot.lnor.legend()
dev.off()
