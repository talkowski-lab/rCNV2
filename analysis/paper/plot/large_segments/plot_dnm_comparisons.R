#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Evaluation of de novo mutations within large segments for rCNV paper


options(stringsAsFactors=F, scipen=1000)
csqs <- c("lof", "mis", "syn")


######################
### DATA FUNCTIONS ###
######################
# Load table of mutation rates
load.mutrates <- function(mutrates.in){
  mutrates <- read.table(mutrates.in, header=T, sep="\t", check.names=F, comment.char="")
  colnames(mutrates)[1] <- "gene"
  return(mutrates)
}

# Load count of de novo mutations for a single study and residualize vs. expected
load.dnms <- function(dnms.in, mutrates){
  dnms <- read.table(dnms.in, header=T, sep="\t", check.names=F, comment.char="")
  colnames(dnms)[1] <- "gene"
  dnms <- merge(dnms, mutrates, all.x=T, all.y=F, sort=F, by="gene")
  for(csq in csqs){
    counts <- dnms[, which(colnames(dnms) == csq)]
    mus <- dnms[, which(colnames(dnms) == paste("mu", csq, sep="_"))]
    n.dnms <- sum(counts, na.rm=T)
    scaled.mus <- (mus / sum(mus, na.rm=T)) * n.dnms
    residuals <- counts - scaled.mus
    dnms[paste("excess", csq, sep="_")] <- residuals
  }
  return(dnms)
}

# Combine DNM data from two cohorts
combine.dnms <- function(d1, d2, mutrates){
  d12 <- aggregate(. ~ gene, rbind(d1, d2), sum)
  d12[, grep("^mu_", colnames(d12))] <- NULL
  d12 <- merge(d12, mutrates, by="gene", all.x=T, all.y=F, sort=F)
  d12[, colnames(d1)]
}

# Calculate oligogenicity index (O_i)
# Defined as the minimum number of genes that captures X% of the total DNM excess in a segment
calc.oligo.index <- function(segs, dnms, csq, cutoff=0.90){
  sapply(segs$genes, function(genes){
    x <- dnms[which(dnms$gene %in% genes), which(colnames(dnms) == paste("excess", csq, sep="_"))]
    x <- x[which(!is.na(x) & x>0)]
    if(length(x) > 0){
      min(which(cumsum(x[order(-x)])/sum(x) >= cutoff))
    }else{
      NA
    }
  })
}

# Calculate cdf matrix of de novo residuals per segment
calc.dnm.excess.cdf.matrix <- function(segs, dnms, csq, n.max.genes=10){
  res <- as.data.frame(t(sapply(segs$genes, function(genes){
    x <- dnms[which(dnms$gene %in% genes), which(colnames(dnms) == paste("excess", csq, sep="_"))]
    if(length(x) > 0){
      x[which(x<0 | is.na(x))] <- 0
      x <- x[order(-x)]
      cdf <- sapply(1:n.max.genes, function(i){
        sum(head(x, i))
      })
      cdf <- c(cdf, sum(x))
    }else{
      cdf <- rep(NA, times=n.max.genes+1)
    }
    return(cdf)
  })))
  rownames(res) <- segs$region_id
  colnames(res) <- c(paste(1:n.max.genes, "genes", sep="_"), "all_genes")
  res[order(-res[, ncol(res)]), ]
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Plot matrix of excess de novo cdfs for all loci
dnm.excess.cdf.barplots <- function(segs, dnms, csq, n.max.genes=5, norm=F, min.excess=1,
                                    subset_to_regions=NULL,
                                    sort.firstwo.only=FALSE,
                                    xtitle.suffix="Segments", ytitle=NULL,
                                    legend=F, cnv.marker.wex=0.035,
                                    parmar=c(1.1, 2.6, 0.4, 0.1)){
  # Subset segments if optioned
  if(!is.null(subset_to_regions)){
    keep.idx <- which(segs$region_id %in% subset_to_regions)
    segs <- segs[keep.idx, ]
  }
  # Helper function for ordering segments by decile
  order.by.decile <- function(res){
    deciles <- seq(1, 0, -0.1)
    ranks <- sapply(deciles, function(d){
      apply(res, 1, function(vals){min(which(vals>=d))})
    })
    order(ranks[, 1], ranks[, 2], ranks[, 3], ranks[, 4], ranks[, 5], ranks[, 6],
          ranks[, 7], ranks[, 8], ranks[, 9], ranks[, 10], ranks[, 11])
  }
  # # Helper function for ordering segments by quartile
  # order.by.quartile <- function(res){
  #   quartiles <- seq(1, 0, -0.2)
  #   ranks <- sapply(quartiles[-1], function(d){
  #     apply(res, 1, function(vals){min(which(vals>=d))})
  #   })
  #   order(ranks[, 1], ranks[, 2], ranks[, 3], ranks[, 4])
  # }

  # Get plot data
  res <- calc.dnm.excess.cdf.matrix(segs, dnms, csq, n.max.genes)
  if(norm==T){
    res <- res[which(apply(res, 1, sum) > min.excess), ]
    res <- t(apply(res, 1, function(vals){vals/max(vals, na.rm=T)}))
    if(sort.firstwo.only==FALSE){
      res <- res[order(res[, 5], res[, 4], res[, 3], res[, 2], res[, 1], decreasing=T), ]
      res <- res[order.by.decile(res), ]
    }else{
      res <- res[order(res[, 2], res[, 1], decreasing=T), ]
      res <- res[order.by.decile(res[, c(1, 2, ncol(res))]), ]
    }
  }
  res <- res[which(!is.na(res[, 1])), ]
  colors <- rev(viridis(n.max.genes + 1))
  ymax <- max(res, na.rm=T)

  # Prep plot area
  par(bty="n", mar=parmar)
  plot(NA, xlim=c(-1, nrow(res)+1), ylim=c(-cnv.marker.wex*ymax, ymax),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")
  # rect(xleft=par("usr")[1], xright=par("usr")[2],
  #      ybottom=par("usr")[3], ytop=par("usr")[4],
  #      border=NA, bty="n", col=bluewhite)
  # abline(h=axTicks(2), col="white")

  # Add rectangles
  sapply(1:nrow(res), function(i){
    rect(xleft=i-1, xright=i,
         ybottom=rev(c(0, res[i, -ncol(res)])), ytop=rev(res[i, ]),
         col=rev(colors), border=rev(colors), xpd=T, lwd=0.25, lend="butt")
    rect(xleft=i-1, xright=i, ybottom=-cnv.marker.wex*ymax, ytop=0,
         col=cnv.colors[segs$cnv[which(segs$region_id==rownames(res)[i])]],
         border=cnv.colors[segs$cnv[which(segs$region_id==rownames(res)[i])]],
         lwd=0.1, xpd=T)
  })

  # Add axes
  abline(h=0, col=blueblack)
  axis(1, at=c(-10e10, 10e10), col=blueblack, labels=NA)
  mtext(1, line=0.15, text=parse(text=paste(prettyNum(nrow(res)), "~", xtitle.suffix)))
  axis(2, at=c(0, 10e10), col=blueblack, labels=NA, tck=0)
  axis(2, at=axTicks(2), labels=NA, col=blueblack, tck=-0.025)
  axis(2, at=axTicks(2), tick=F, las=2, line=-0.65)
  mtext(2, line=1.55, text=ytitle)

  # Add legend, if optioned
  if(legend==T){
    legend.labs <- c("Top Gene", expression(2^"nd" ~ "Gene"), expression(3^"rd" ~ "Gene"),
                     sapply(4:100, function(x){bquote(.(x)^"th" ~ "Gene")}))
    legend.labs <- c(legend.labs[1:n.max.genes], "All Others")
    legend("topright", fill=colors, legend=legend.labs, bty="n", border=NA)
  }
}


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(optparse, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog loci.bed segs.tsv ddd.tsv asc.tsv asc.control.tsv mutrates.tsv out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 7){
  stop(paste("Seven positional arguments required: loci.bed, segs.tsv, ddd.tsv, asc.tsv, asc.control.tsv, mutrates.tsv, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
loci.in <- args$args[1]
segs.in <- args$args[2]
ddd.dnms.in <- args$args[3]
asc.dnms.in <- args$args[4]
asc.control.dnms.in <- args$args[5]
gene.mutrates.in <- args$args[6]
out.prefix <- args$args[7]

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d2.master_segments.bed.gz"
# ddd.dnms.in <- "~/scratch/ddd_dnm_counts.tsv.gz"
# asc.dnms.in <- "~/scratch/asc_dnm_counts.tsv.gz"
# asc.control.dnms.in <- "~/scratch/asc_dnm_counts.unaffecteds.tsv.gz"
# gene.mutrates.in <- "~/scratch/gene_mutation_rates.tsv.gz"
# out.prefix <- "~/scratch/test_oligogenicity"

# Load dnm counts and residualize by mutation rates
mutrates <- load.mutrates(gene.mutrates.in)
ddd <- load.dnms(ddd.dnms.in, mutrates)
asc <- load.dnms(asc.dnms.in, mutrates)
asc.control <- load.dnms(asc.control.dnms.in, mutrates)
ddd.plus.asc <- combine.dnms(ddd, asc, mutrates)
dnms <- list("DDD"=ddd, "ASC"=asc, "ASC_unaffected"=asc.control,
             "DDD_plus_ASC"=ddd.plus.asc)

# Load loci & segment table
loci <- load.loci(loci.in)
segs.all <- load.segment.table(segs.in)

# Create subset of all GD segs & those with nominal significance
segs.all <- segs.all[which(segs.all$any_gd | segs.all$any_sig), ]

# Make analysis subset of only discovery segments at GW or FDR, or lit GDs at Bonferroni
segs <- segs.all[which(segs.all$any_sig | segs.all$bonf_sig_gd), ]

# For DNM analyses, subset strictly to NDD segments (to match DDD & ASC phenotypes)
dev.seg.ids <- get.developmental.region_ids(loci, segs)
neuro.seg.ids <- get.neuro.region_ids(loci, segs)
ndd.seg.ids <- intersect(dev.seg.ids, neuro.seg.ids)
ndd.segs <- segs[which(segs$region_id %in% ndd.seg.ids), ]

# Annotate segs with oligogenicity indexes
for(csq in csqs){
  for(cohort in names(dnms)){
    ndd.segs[paste(cohort, csq, "oligo", sep=".")] <- calc.oligo.index(ndd.segs, dnms[[cohort]], csq=csq)
  }
}

# Swarmplots of average excess DNMs per gene split by CNV type
sapply(names(dnms), function(cohort){
  sapply(c("lof", "mis"), function(csq){
    if(csq=="lof"){
      csq.name <- "PTV"
    }else if(csq=="mis"){
      csq.name <- "Mis."
    }else{
      csq.name <- "Syn."
    }
    pdf(paste(out.prefix, "sig_plus_litGDs", cohort, csq, "excess_per_gene.NDD_only.pdf", sep="."),
        height=2, width=1.6)
    segs.simple.vioswarm(ndd.segs,
                         y=ndd.segs[, paste(cohort, "dnm", csq, "norm_excess_per_gene", sep="_")],
                         add.pvalue=T,
                         ytitle=paste("Excess dn", csq.name, " per Gene", sep=""),
                         pt.cex=0.4,
                         parmar=c(1.2, 2.65, 1.5, 0))
    dev.off()
  })
})

# Swarmplots of oligogenicity index for all deletions & duplications
sapply(names(dnms), function(cohort){
  sapply(csqs, function(csq){
    pdf(paste(out.prefix, "sig_plus_litGDs", cohort, csq, "oligogenicity_index.NDD_only.pdf", sep="."),
        height=2, width=1.6)
    segs.simple.vioswarm(ndd.segs, y=ndd.segs[, paste(cohort, csq, "oligo", sep=".")], add.pvalue=T,
                         ytitle="Oligogenicity Index",
                         pt.cex=0.4, parmar=c(1.2, 2.65, 1.5, 0))
    dev.off()
  })
})

# Scatterplots of oligogenicity index vs number of genes in segment for all deletions & duplications
abline.slopes <- c(1/2, 1/4, 1/10)
abline.labels <- c("50%", "25%", "10%")
sapply(c("ASC", "DDD"), function(cohort){
  sapply(csqs, function(csq){
    pdf(paste(out.prefix, "sig_plus_litGDs", cohort, csq, "oligogenicity_index_vs_ngenes.NDD_only.pdf", sep="."),
        height=2, width=2.2)
    segs.scatter(ndd.segs, x=ndd.segs$n_genes, add.cor=F,
                 y=ndd.segs[, paste(cohort, csq, "oligo", sep=".")],
                 ylims=c(0, ceiling(0.5*max(ndd.segs$n_genes, na.rm=T))),
                 xtitle="Genes in Segment", x.title.line=1.3,
                 ytitle="Oligogenicity Index", y.title.line=1.5,
                 abline.a=rep(1, length(abline.slopes)),
                 abline.b=abline.slopes,
                 abline.lty=rep(2, length(abline.slopes)),
                 blue.bg=FALSE, add.lm=T, pt.cex=0.4,
                 parmar=c(2.3, 2.5, 0.6, 1.8))
    sapply(1:length(abline.slopes), function(i){
      axis(4, at=(par("usr")[2]*abline.slopes[i])+1, tck=-0.02,
           labels=NA, col=blueblack, xpd=T)
      axis(4, at=(par("usr")[2]*abline.slopes[i])+1, tick=F, labels=abline.labels[i],
           line=-0.7, las=2, xpd=T, col.axis=blueblack, cex.axis=0.8)
    })
    dev.off()
  })
})

# Swarmplot of DDD oligogenicity indexes vs. mechanism
pdf(paste(out.prefix, "sig_plus_litGDs.segs_by_mechanism.ddd_lof_oligogenicity_index.NDD_only.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(ndd.segs, x.bool=ndd.segs$nahr, y=ndd.segs$ddd.lof.oligo,
           x.labs=c("Non-NAHR", "NAHR"), violin=T, add.pvalue=T,
           add.y.axis=T, ytitle="Oligogenicity Index", pt.cex=0.4,
           parmar=c(1.2, 3, 2.5, 0))
dev.off()
pdf(paste(out.prefix, "sig_plus_litGDs.segs_by_mechanism.ddd_mis_oligogenicity_index.NDD_only.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(ndd.segs, x.bool=ndd.segs$nahr, y=ndd.segs$ddd.mis.oligo,
           x.labs=c("Non-NAHR", "NAHR"), violin=T, add.pvalue=T,
           add.y.axis=T, ytitle="Oligogenicity Index", pt.cex=0.4,
           parmar=c(1.2, 3, 2.5, 0))
dev.off()

# Swarmplot of ASC oligogenicity indexes vs. mechanism
pdf(paste(out.prefix, "sig_plus_litGDs.segs_by_mechanism.asc_lof_oligogenicity_index.NDD_only.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(ndd.segs, x.bool=ndd.segs$nahr, y=ndd.segs$asc.lof.oligo,
           x.labs=c("Non-NAHR", "NAHR"), violin=T, add.pvalue=T,
           add.y.axis=T, ytitle="Oligogenicity Index", pt.cex=0.4,
           parmar=c(1.2, 3, 2.5, 0))
dev.off()
pdf(paste(out.prefix, "sig_plus_litGDs.segs_by_mechanism.asc_mis_oligogenicity_index.NDD_only.pdf", sep="."),
    height=2.25, width=2.6)
segs.swarm(ndd.segs, x.bool=ndd.segs$nahr, y=ndd.segs$asc.mis.oligo,
           x.labs=c("Non-NAHR", "NAHR"), violin=T, add.pvalue=T,
           add.y.axis=T, ytitle="Oligogenicity Index", pt.cex=0.4,
           parmar=c(1.2, 3, 2.5, 0))
dev.off()

# Distribution of excess DNMs across neuro segments
sapply(names(dnms), function(cohort){
  if(cohort %in% c("DDD", "DDD_plus_ASC")){
    min.excess <- 3
  }else{
    min.excess <- 1
  }
  sapply(csqs, function(csq){
    if(csq=="lof"){
      ytitle.1 <- bquote("Excess" ~ italic("dn") * "PTVs")
      ytitle.2 <- bquote("Prop. of Excess" ~ italic("dn") * "PTVs")
      xtitle.suffix.2 <- paste('"Segs. with " >= "', min.excess, ' Excess" ~ italic("dn") * "PTV"', sep="")
    }else if(csq=="mis"){
      ytitle.1 <- bquote("Excess" ~ italic("dn") * "Mis.")
      ytitle.2 <- bquote("Prop. of Excess" ~ italic("dn") * "Mis.")
      xtitle.suffix.2 <- paste('"Segs. with " >= "', min.excess, ' Excess" ~ italic("dn") * "Mis."', sep="")
    }else{
      ytitle.1 <- bquote("Excess" ~ italic("dn") * "Syn.")
      ytitle.2 <- bquote("Prop. of Excess" ~ italic("dn") * "Syn.")
      xtitle.suffix.2 <- paste('"Segs. with " >= "', min.excess, ' Excess" ~ italic("dn") * "Syn."', sep="")
    }
    pdf(paste(out.prefix, "sig_plus_litGDs", cohort, csq, "excess_dnm_distrib_bygene.NDD_only.pdf", sep="."),
        height=2, width=3)
    dnm.excess.cdf.barplots(ndd.segs, dnms[[cohort]], csq=csq, norm=F, legend=T,
                            xtitle.suffix='"NDD rCNV Segments"',
                            ytitle=ytitle.1)
    dev.off()
    pdf(paste(out.prefix, "sig_plus_litGDs", cohort, csq, "excess_dnm_distrib_bygene.norm.NDD_only.pdf", sep="."),
        height=2, width=3)
    dnm.excess.cdf.barplots(ndd.segs, dnms[[cohort]], csq=csq, norm=T, legend=F,
                            min.excess=min.excess, xtitle.suffix=xtitle.suffix.2,
                            ytitle=ytitle.2)
    dev.off()
  })
})

