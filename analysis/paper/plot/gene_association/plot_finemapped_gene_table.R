#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Create annotated 2x2 table-like figure of fine-mapped genes for rCNV2 paper


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Load a list of OMIM gene lists per HPO
load.omim.lists <- function(omimlists.in){
  # Expects a two-column tsv of (hpo, path) pairs
  olist <- read.table(omimlists.in, sep="\t", header=F)
  hpos <- as.character(olist[, 1])
  glists <- lapply(1:nrow(olist), function(i){
    unique(sort(as.character(read.table(olist[i, 2], header=F)[, 1])))
  })
  names(glists) <- hpos
  return(glists)
}

# Annotate a single gene based on OMIM, HPO, and constraint
annotate.gene <- function(gene, hpos, omim.genes, lof.genes, mis.genes){
  # OMIM comparisons
  omim.matches <- sapply(omim.genes, function(glist){gene %in% glist})
  if(any(omim.matches[hpos])){
    omim.label <- "D+"
  }else if(any(omim.matches)){
    omim.label <- "D"
  }else{
    omim.label <- NA
  }

  # Constraint comparisons
  if(gene %in% lof.genes){
    lof <- TRUE
  }else{
    lof <- FALSE
  }
  if(gene %in% mis.genes){
    mis <- TRUE
  }else{
    mis <- FALSE
  }

  return(list("omim"=omim.label, "lof"=lof, "mis"=mis))
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Format gene symbol for plotting
format.gene.label <- function(gene, assocs, omim.genes, lof.genes, mis.genes){
  # Get data
  hpos <- assocs$hpo[which(assocs$gene == gene)]
  annos <- annotate.gene(gene, hpos, omim.genes, lof.genes, mis.genes)

  # Set base
  base <- paste("italic('", toupper(gene), "')", sep="")
  if(annos$lof==TRUE){
    base <- paste("bold", base, sep="")
  }
  if(annos$mis==TRUE){
    base <- paste("underline(", base, ")", sep="")
  }

  # Set superscript
  super.vals <- c(annos$omim)
  if(any(!is.na(super.vals))){
    super <- paste(super.vals[which(!is.na(super.vals))], sep=", ")
    gene.fmt <- paste(base, "^'", super, "'", sep="")
  }else{
    gene.fmt <- base
  }

  return(gene.fmt)
}

# Get gene color by quadrant
get.gene.color.byquadrant <- function(gene, gene.groups, top, conf){
  quad.name <- paste(top, conf, sep=".")
  is.del <- gene %in% gene.groups$DEL[[quad.name]]
  is.dup <- gene %in% gene.groups$DUP[[quad.name]]
  if(is.del==TRUE){
    if(is.dup==TRUE){
      return(cnv.colors[3])
    }else{
      return(cnv.colors[1])
    }
  }else if(is.dup==TRUE){
    return(cnv.colors[2])
  }else{
    return(NA)
  }
}

# Print genes as text to a single quadrant
print.quadrant <- function(genes, gene.groups, credsets, assocs,
                           omim.genes, lof.genes, mis.genes,
                           x.start, y.start, top, conf,
                           text.cex=5/6, break.width=5,
                           max.n.genes=50, y.min=NULL){
  y.add <- 0
  x.add <- 0
  k <- 0
  if(length(genes) > 0){
    if(length(genes) < max.n.genes){
      for(i in 1:length(genes)){
        gene <- genes[i]
        gene.fmt <- format.gene.label(gene, assocs, omim.genes, lof.genes, mis.genes)
        gene.col <- get.gene.color.byquadrant(gene, gene.groups, top, conf)
        text(x=x.start + x.add, y=y.start + y.add + 0.5,
             labels=parse(text=gene.fmt), pos=4, cex=text.cex, col=gene.col)
        x.add <- x.add + 1
        k <- k + 1
        if(k >= break.width){
          k <- 0
          x.add <- 0
          y.add <- y.add + 1
        }
      }
    }else{
      gene.counts <- table(sapply(genes, get.gene.color.byquadrant, gene.groups=gene.groups, top=top, conf=conf))
      sum.labels <- paste(gene.counts[cnv.colors], "genes")
      text(x=x.start + (1:3) - 1, y=y.start + 0.5,
           labels=sum.labels, pos=4, cex=text.cex, col=cnv.colors)
    }
  }else{
    text(x=x.start + (break.width/2), y=(y.start + y.min)/2,
         labels="(No genes)", cex=text.cex, col=ns.color)
  }
}

# Master function to plot 2x2 gene grid
plot.gene.grid <- function(gene.groups, credsets, assocs,
                           omim.genes, lof.genes, mis.genes,
                           quadrant.width=5, max.per.quadrant=50,
                           parmar=c(0.1, 6, 1.5, 0.1)){
  # Get genes by quadrant
  genes.topleft <- get.quadrant.genes(gene.groups, "top", "vconf")
  genes.topright <- get.quadrant.genes(gene.groups, "nottop", "vconf")
  genes.bottomleft <- get.quadrant.genes(gene.groups, "top", "conf")
  genes.bottomright <- get.quadrant.genes(gene.groups, "nottop", "conf")
  genes.otherleft <- get.quadrant.genes(gene.groups, "top", "notconf")
  genes.otherright <- get.quadrant.genes(gene.groups, "nottop", "notconf")

  # Gather plot dimensions
  n.cols <- 2*quadrant.width
  n.rows.top <- max(ceiling(sapply(list(genes.topleft, genes.topright), length) / quadrant.width))
  n.rows.bottom <- max(ceiling(sapply(list(genes.bottomleft, genes.bottomright), length) / quadrant.width))
  n.rows.total <- n.rows.top + n.rows.bottom + 1

  # Prepare plot area
  par(mar=parmar, bty="n")
  plot(NA, xlim=c(0, n.cols), ylim=c(n.rows.total, 0),
       xaxt="n", xlab="", xaxs="i", yaxt="n", ylab="", yaxs="i")

  # Add grid
  rect(xleft=0, xright=quadrant.width, ytop=0, ybottom=n.rows.top+n.rows.bottom, xpd=T,
       col="white", border=NA, bty="n", lwd=1.5)
  rect(xleft=0, xright=quadrant.width, ytop=0, ybottom=n.rows.top+n.rows.bottom, xpd=T,
       col=adjustcolor(highlight.color, alpha=0.3), border=NA, bty="n", lwd=1.5)
  rect(xleft=0, xright=n.cols, ytop=0, ybottom=n.rows.total, xpd=T,
       col=NA, border=blueblack, lwd=1.5)
  segments(x0=0, x1=n.cols,
           y0=c(n.rows.top, n.rows.top+n.rows.bottom),
           y1=c(n.rows.top, n.rows.top+n.rows.bottom),
           col=blueblack)
  segments(x0=quadrant.width, x1=quadrant.width, y0=0, y1=n.rows.total, col=blueblack)

  # Add axes
  axis(3, at=quadrant.width/2, tick=F, line=-0.85, labels="Top gene in at least one credible set")
  axis(3, at=(3/2)*quadrant.width, tick=F, line=-0.85, labels="Not top gene in any credible set")
  axis(2, at=(n.rows.top/2)-0.9, tick=F, line=-0.8, las=2, labels="Highly")
  axis(2, at=(n.rows.top/2), tick=F, line=-0.8, las=2, labels="confident")
  axis(2, at=(n.rows.top/2)+0.9, tick=F, line=-0.8, las=2,
       labels=expression(("PIP" >= 0.85)), cex.axis=(5.5/6))
  axis(2, at=n.rows.top+(n.rows.bottom/2)-0.2, tick=F, line=-0.8, las=2, padj=0,
       labels="Confident")
  axis(2, at=n.rows.top+(n.rows.bottom/2)+0.2, tick=F, line=-0.8, las=2, padj=1,
       labels=expression(("PIP" >= 0.15)), cex.axis=(5.5/6))
  axis(2, at=n.rows.top+n.rows.bottom+0.5-0.2, tick=F, line=-0.8, las=2, padj=0,
       labels="Unlikely")
  axis(2, at=n.rows.top+n.rows.bottom+0.5+0.2, tick=F, line=-0.8, las=2, padj=1,
       labels=expression(("PIP" < 0.15)), cex.axis=(5.5/6))

  # Add genes to each quadrant
  print.quadrant(genes.topleft, gene.groups, credsets, assocs, omim.genes, lof.genes, mis.genes,
                 x.start=0.08, y.start=0, "top", "vconf", break.width=quadrant.width,
                 y.min=n.rows.top, max.n.genes=max.per.quadrant)
  print.quadrant(genes.topright, gene.groups, credsets, assocs, omim.genes, lof.genes, mis.genes,
                 x.start=quadrant.width+0.08, y.start=0, "nottop", "vconf",
                 break.width=quadrant.width, y.min=n.rows.top, max.n.genes=max.per.quadrant)
  print.quadrant(genes.bottomleft, gene.groups, credsets, assocs, omim.genes, lof.genes, mis.genes,
                 x.start=0.08, y.start=n.rows.top, "top", "conf",
                 break.width=quadrant.width, y.min=n.rows.top+n.rows.bottom,
                 max.n.genes=max.per.quadrant)
  print.quadrant(genes.bottomright, gene.groups, credsets, assocs, omim.genes, lof.genes, mis.genes,
                 x.start=quadrant.width+0.08, y.start=n.rows.top, "nottop", "conf",
                 break.width=quadrant.width, y.min=n.rows.top+n.rows.bottom,
                 max.n.genes=max.per.quadrant)
  print.quadrant(genes.otherleft, gene.groups, credsets, assocs, omim.genes, lof.genes, mis.genes,
                 x.start=0.08, y.start=n.rows.top+n.rows.bottom, "top", "notconf",
                 break.width=quadrant.width, y.min=n.rows.total)
  print.quadrant(genes.otherright, gene.groups, credsets, assocs, omim.genes, lof.genes, mis.genes,
                 x.start=quadrant.width+0.08, y.start=n.rows.top+n.rows.bottom, "nottop", "notconf",
                 break.width=quadrant.width, y.min=n.rows.total)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(optparse, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog credsets.bed assocs.bed",
                                            "lof.constr.list mis.constr.list",
                                            "omim.lists out.prefix"),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 6){
  stop(paste("Six positional arguments required: credsets.bed, assocs.bed,",
             "lof.constrained.list, mis.constrained.list, omim.lists, and out.prefix\n"))
}

# Writes args & opts to vars
credsets.in <- args$args[1]
assocs.in <- args$args[2]
lof.in <- args$args[3]
mis.in <- args$args[4]
omimlists.in <- args$args[5]
out.prefix <- args$args[6]

# # DEV PARAMETERS
# credsets.in <- "~/scratch/rCNV.final_genes.credible_sets.bed.gz"
# assocs.in <- "~/scratch/rCNV.final_genes.associations.bed.gz"
# lof.in <- "~/scratch/gene_lists/gnomad.v2.1.1.lof_constrained.genes.list"
# mis.in <- "~/scratch/gene_lists/gnomad.v2.1.1.mis_constrained.genes.list"
# omimlists.in <- "~/scratch/omim.gene_lists.tsv"
# out.prefix <- "~/scratch/test_finemapped_genes_table"

# Load credible sets and associations
credsets <- load.credsets(credsets.in)
assocs <- load.gene.associations(assocs.in)

# Load gene lists for annotations
lof.genes <- read.table(lof.in, header=F)[, 1]
mis.genes <- read.table(mis.in, header=F)[, 1]
omim.genes <- load.omim.lists(omimlists.in)

# Split fine-mapped genes into categories based on top/not top status and conf/vconf
gene.groups <- categorize.genes(credsets)

# Plot gene grid for all genes
pdf(paste(out.prefix, "finemapped_genes_grid.all.pdf", sep="."), height=2*6, width=2*6.5)
plot.gene.grid(gene.groups, credsets, assocs, omim.genes, lof.genes, mis.genes,
               quadrant.width=5, max.per.quadrant=10e10, parmar=c(0.5, 4.5, 1.2, 0.2))
dev.off()

## TODO: ADD SPLITS BY SIGNIFICANCE & ADULT/DEV
