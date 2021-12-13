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
# Compute excess DNM density of all possible N-gene segments
enumerate.all.excess.dnms <- function(segs, dnms, csq, gbed){
  max.n.genes <- ceiling(1.25 * max(segs$n_genes, na.rm=T))

  # Reorder genes by linear genomic position and split by chromosome
  dnms.ordered <- merge(gbed, dnms, by="gene", sort=F, all=F)
  dnms.chromsplit <- lapply(sort(unique(dnms.ordered$chrom)), function(contig){
    dnms.ordered[which(dnms.ordered$chrom == contig), ]
    })

  # Compute rolling sum of excess DNMs for all possible k-gene blocks
  lapply(1:max.n.genes, function(k){
    as.numeric(unlist(lapply(dnms.chromsplit, function(df){
      vals <- rollapply(df[, paste("excess", csq, sep="_")], k, sum, fill=NA, partial=F)
      return(vals[which(!is.na(vals))])
    })))
  })
}

# Assign DNM excess percentile to each segment
assign.excess.pct <- function(segs, enumerated.excesses){
  for(cohort in names(enumerated.excesses)){
    ex.allcsq <- enumerated.excesses[[cohort]]
    for(csq in names(ex.allcsq)){
      ex <- ex.allcsq[[csq]]
      pcts <- sapply(1:nrow(segs), function(i){
        n.genes <- segs$n_genes[i]
        if(n.genes == 0){
          return(NA)
        }else{
          excess <- segs[i, paste(cohort, "dnm", csq, "obs_wMu", sep="_")] - segs[i, paste(cohort, "dnm", csq, "exp_wMu", sep="_")]
          null.dist <- ex[[n.genes]]
          return(length(which(excess >= null.dist)) / length(null.dist))
        }
      })
      segs[, paste(cohort, "dnm", csq, "excess_pct", sep="_")] <- pcts
    }
  }
  return(segs)
}

# Quantify enrichment of segments in high percentiles of DNM excess
quantify.excesses <- function(segs, outfile, quantiles=c(0.5, 0.75, 0.9, 0.95)){
  res <- do.call("rbind", lapply(colnames(segs)[grep("_excess_pct", colnames(segs), fixed=T)],
                          function(cname){
    do.call("rbind", lapply(list(c("DEL"), c("DUP"), c("DEL", "DUP")), function(cnv.types){
      n.segs <- length(which(segs$cnv %in% cnv.types))

      # Compute enrichments per quantile
      res <- do.call("rbind", lapply(quantiles, function(q){
        n.pass <- length(which(segs[, cname] >= q & segs$cnv %in% cnv.types))
        stats <- binom.test(n.pass, n.segs, p=(1-q), alternative="greater")
        c(cname, paste(cnv.types, collapse="_"), "binomial", q, n.segs, (1-q)*n.segs, n.pass,
          stats$p.value, stats$estimate / (1-q))
      }))

      # Perform K-S test to evaluate overall deviation from expectation (uniform)
      stats <- ks.test(segs[which(segs$cnv %in% cnv.types), cname], "punif")
      rbind(res, c(cname, paste(cnv.types, collapse="_"), "K-S", NA, n.segs,
                   NA, NA, stats$p.value, NA))
    }))
  }))
  res <- as.data.frame(res)
  colnames(res) <- c("feature", "cnv", "test", "quantile", "n_segs", "expected",
                     "observed", "p_value", "fold_enrichment")
  write.table(res, outfile, sep="\t", col.names=T, row.names=F, quote=F)
}

# Gather ordered vector of DNM excesses for a single segment
get.dnm.excess.byGene <- function(segs, dnms, csq, region.id){
  genes <- segs$genes[[which(segs$region_id == region.id)]]
  x <- dnms[which(dnms$gene %in% genes), c("gene", paste("excess", csq, sep="_"))]
  vals <- as.numeric(x[, 2])
  names(vals) <- as.character(x[, 1])
  vals[order(-vals)]
}

# Calculate cdf matrix of de novo residuals per segment
calc.dnm.excess.cdf.matrix <- function(segs, dnms, csq, n.max.genes=10){
  res <- as.data.frame(t(sapply(segs$region_id, function(region.id){
    x <- get.dnm.excess.byGene(segs, dnms, csq, region.id)
    if(length(x) > 0){
      x[which(x<0 | is.na(x))] <- 0
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
# Scatterplot of excess DNMs per segment vs. expected
plot.excess.dnms.vs.expected <- function(segs, enumerated.excesses, cohort, csq,
                                         quantiles=c(0.5, 0.75, 0.9, 0.95),
                                         label.top.n=12, min.label.quantile=0.9,
                                         pt.cex=0.85, y.lims=NULL,
                                         parmar=c(2.2, 2.7, 0.25, 2.45)){
  # Get core data
  segs <- segs[which(segs$n_genes > 0), ]
  quantiles <- sort(quantiles, decreasing=T)
  en.ex <- enumerated.excesses[[cohort]][[csq]]
  max.n.genes <- max(segs$n_genes, na.rm=T)
  segs$plot.vals <- segs[, paste(cohort, "dnm", csq, "obs_wMu", sep="_")] - segs[, paste(cohort, "dnm", csq, "exp_wMu", sep="_")]
  expected.curves <- lapply(quantiles, function(q){
    data.frame("n_genes"=1:length(en.ex),
               "excess"=sapply(en.ex, quantile, na.rm=T, probs=q))
  })
  names(expected.curves) <- quantiles

  # Get plot parameters
  if(is.null(y.lims)){
    y.lims <- range(c(segs$plot.vals, unlist(sapply(expected.curves, function(df){df[, 2]}))), na.rm=T)
    y.lims[2] <- 1.05 * y.lims[2]
  }
  csq.abbrev <- c("lof" = "PTVs", "mis" = "Missense", "syn" = "Syn.")[csq]
  ytitle <- bquote("Excess" ~ italic("De Novo") ~ .(csq.abbrev))
  pct.grad <- colorRampPalette(c("white", ns.color))(length(quantiles) + 2)
  pct.lab.grad <- colorRampPalette(c("gray30", ns.color))(length(quantiles))
  pct.cex <- seq(1.25, 0.6, length.out=length(quantiles) + 1) * pt.cex
  seg.pcts <- segs[, paste(cohort, "dnm", csq, "excess_pct", sep="_")]
  pt.cex <- pct.cex[sapply(seg.pcts, function(p){min(which(p - c(quantiles, 0) >= 0))})]
  segs$pt.cex <- pt.cex

  # Prep plot area
  par(bty="n", mar=parmar)
  plot(NA, xlim=c(1, max.n.genes + 2), ylim=y.lims,
       xaxt="n", xlab="", xaxs="i", yaxt="n", ylab="")

  # Add expected quantiles
  lapply(1:length(quantiles), function(i){
    q <- quantiles[i]
    qpct <- round(100*q, 1)
    qlab <- bquote(.(qpct) ^ "th")
    ex <- expected.curves[[which(names(expected.curves) == q)]]
    polygon(x=c(ex[, 1], 10e10, -10e10), y=c(ex[, 2], -10e10, -10e10),
            border=NA, col=pct.grad[i])
    points(ex, type="l", col=ns.color)
    axis(4, at=ex[max.n.genes + 2, 2] + 0.02*diff(par("usr")[3:4]),
         labels=qlab, tick=F, las=2, cex.axis=5/6, line=-0.8,
         col.axis=pct.lab.grad[i])
  })

  # Add axes
  x.ax.at <- sort(unique(c(1, axTicks(1))))
  axis(1, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(1, at=x.ax.at, labels=NA, tck=-0.025, col=blueblack)
  axis(1, at=x.ax.at, tick=F, line=-0.7)
  y.ax.at <- axTicks(2)
  axis(2, at=c(-10e10, 10e10), tck=0, labels=NA, col=blueblack)
  axis(2, at=y.ax.at, labels=NA, tck=-0.025, col=blueblack)
  axis(2, at=y.ax.at, tick=F, line=-0.65, las=2)

  # Add axis titles
  mtext(1, text="Genes per Segment", line=1.1)
  mtext(2, text=ytitle, line=1.7)
  mtext(4, text="Percentile (Genome-Wide)", line=1.45, col=ns.color.dark)

  # Label top N points
  lab.segs <- segs[which(seg.pcts >= min.label.quantile), ]
  sapply(head(order(lab.segs$plot.vals, decreasing=T), label.top.n), function(i){
    x <- lab.segs$n_genes[i]; y <- lab.segs$plot.vals[i]
    if(x > mean(par("usr")[1:2])){lab.pos <- 2}else{lab.pos <- 4}
    lab <- unlist(strsplit(lab.segs$cytoband[i], split="-"))[1]
    text(x=x, y=y, labels=lab, cex=4.5/6, col=lab.segs$pt.border[i], pos=lab.pos, xpd=T)
  })

  # Add points in order by percentile and significance
  sapply((length(quantiles) + 1):1, function(i){
    q.idxs <- which(seg.pcts > c(quantiles, 0)[i] & seg.pcts <= c(1, quantiles)[i])
    plot.idxs <- list(intersect(q.idxs, which(segs$any_gd & !segs$any_sig)),
                      intersect(q.idxs, which(segs$fdr_sig)),
                      intersect(q.idxs, which(segs$gw_sig)))
    lapply(plot.idxs, function(p.idxs){
      segs.tmp <- segs[p.idxs, ]
      points(segs.tmp$n_genes, segs.tmp$plot.vals, xpd=T, pch=segs.tmp$pt.pch,
             bg=segs.tmp$pt.bg, col=segs.tmp$pt.border, cex=segs.tmp$pt.cex)
    })
  })
}

# Plot detailed highlight barplots of DNM excess for top segments
dnm.excess.cdf.highlights <- function(segs, dnms, cohort, csq,
                                      top.n.segs=5, top.n.genes=4,
                                      annotate.weaker.genes=FALSE,
                                      rect.buffer=0.1, parmar=c(1, 7, 2.35, 1)){
  # Subset segments to top.n with strongest DNM excesses
  pct.colname <- paste(cohort, "dnm", csq, "raw_excess_per_gene", sep="_")
  keep.idx <- head(order(-segs[, pct.colname] * segs$n_genes), top.n.segs)
  segs <- segs[keep.idx, ]

  # Get plot data
  res <- lapply(segs$region_id, get.dnm.excess.byGene, segs=segs, dnms=dnms, csq=csq)
  names(res) <- segs$region_id
  bar.pal <- rev(viridis(top.n.genes + 1))
  csq.lab <- c("lof" = "PTVs", "mis" = "Mis.", "syn" = "Syn.")[csq]
  cnv.lab <- c("lof" = "DEL", "mis" = "DUP", "syn" = "rCNV")[csq]

  # Prep plot area
  par(bty="n", mar=parmar)
  plot(NA, xlim=c(0, 1), ylim=c(top.n.segs+rect.buffer, -rect.buffer),
       xaxt="n", xlab="", xaxs="i", yaxt="n", ylab="", yaxs="i")

  # Plot each segment as a single strip
  sapply(1:nrow(segs), function(i){
    # Clean up DNM values
    vals <- res[[i]]
    vals[which(vals < 0 | is.na(vals))] <- 0
    named <- head(vals[which(vals > 0)], top.n.genes)
    nonnamed <- vals[which(!names(vals) %in% names(named))]
    n.nonnamed <- length(nonnamed)
    if(length(nonnamed) > 1){
      nonnamed <- sum(nonnamed)
      names(nonnamed) <- paste(n.nonnamed, "Others")
    }
    vals <- c(named, nonnamed)
    vals.scaled.marginal <- vals / sum(vals)
    vals.scaled <- cumsum(vals) / sum(vals)

    # Plot rectangle
    if(annotate.weaker.genes){
      rect.y.coords <- c((i-1)+rect.buffer, i-0.5)
    }else{
      rect.y.coords <- c((i-1)+rect.buffer, i-rect.buffer)
    }
    rect(xleft=c(0, vals.scaled[-length(vals.scaled)]), xright=vals.scaled,
         ybottom=rect.y.coords[1], ytop=rect.y.coords[2], col=bar.pal, xpd=T)
    rect.mid <- mean(rect.y.coords)

    # Add gene labels
    direct.labs <- which(vals.scaled.marginal >= 0.2)
    sapply(direct.labs, function(j){
      if(j<3){text.col <- "black"}else{text.col <- ns.color.light}
      if(length(grep("Others", names(vals)[j])) == 1){
        lab <- names(vals)[j]
      }else if(vals.scaled.marginal[j] > 0.4){
        lab <- bquote(bolditalic(.(names(vals)[j])) ~ (.(paste("+", round(vals[j], 0), " ", csq.lab, sep=""))))
      }else if(vals.scaled.marginal[j] > 0.3){
        lab <- bquote(bolditalic(.(names(vals)[j])) ~ (.(paste("+", round(vals[j], 0), sep=""))))
      }else{
        lab <- bquote(bolditalic(.(names(vals)[j])))
      }
      text(x=mean(c(c(0, vals.scaled)[j], vals.scaled[j])), y=rect.mid,
           labels=lab, cex=5/6, col=text.col)
    })

    # Add minor gene annotations in margins, if optioned
    indirect.labs <- which(vals.scaled.marginal < 0.2)
    if(annotate.weaker.genes & length(indirect.labs) > 0){
      indirect.xpos <- seq(0.9, 0.1, length.out=top.n.genes)[length(indirect.labs):1]
      sapply(1:length(indirect.labs), function(j){
        lab.idx <- indirect.labs[j]
        if(vals.scaled.marginal[lab.idx] > 0.1 & vals[lab.idx] > 3
           & length(grep("Others", names(vals)[lab.idx])) == 0){
          lab <- bquote(italic(.(names(vals)[lab.idx])) ~ (.(paste("+", round(vals[lab.idx], 0),  sep=""))))
        }else{
          lab <- bquote(italic(.(names(vals)[lab.idx])))
        }
        text(x=indirect.xpos[j], y=i-0.15, labels=lab, cex=5/6, col=ns.color.dark, xpd=T)
        segments(x0=mean(c(0, vals.scaled)[lab.idx + 0:1]),
                 x1=mean(c(0, vals.scaled)[lab.idx + 0:1]),
                 y0=i-0.45, y1=i-0.5, col="black", xpd=T)
        segments(x0=indirect.xpos[j], x1=mean(c(0, vals.scaled)[lab.idx + 0:1]),
                 y0=i-0.35, y1=i-0.45, col=ns.color, xpd=T)
        segments(x0=indirect.xpos[j], x1=indirect.xpos[j],
                 y0=i-0.35, y1=i-0.3, col=ns.color, xpd=T)
      })
    }

    # Add region label
    axis(2, at=rect.mid, tick=F, line=-0.85, las=2, labels=segs$cytoband[i])
  })

  # Add top X axis
  x.ax.at <- axTicks(3)
  axis(3, at=c(-10e10, 10e10), col=blueblack, tck=0, labels=NA)
  axis(3, at=x.ax.at, col=blueblack, tck=-0.025, labels=NA)
  sapply(x.ax.at, function(x){
    axis(3, at=x, tick=F, line=-0.7, cex.axis=5.5/6,
         labels=paste(round(100 * x, 0), "%", sep=""))
  })
  x.ax.title <- bquote("Prop. of Excess" ~ italic("dn") * .(csq.lab) ~ "per" ~ .(cnv.lab) ~ "Segment")
  mtext(3, line=1.25, text=x.ax.title)

  # Add bottom X axis
  x.ax.at <- seq(0, 1, 0.25)
  axis(1, tck=-0.015, col=ns.color.light, labels=NA, xpd=T, line=0.35)
  axis(1, tck=0.015, col=ns.color.light, labels=NA, xpd=T, line=0.35)
  mtext(1, text="Genes per Segment Ranked by DNM Excess",
        cex=5/6, font=3, col=ns.color.dark, line=-0.25)
}

# Plot matrix of excess de novo cdfs for all loci
dnm.excess.cdf.barplots <- function(segs, dnms, cohort, csq, n.max.genes=5, norm=F,
                                    pct.cutoff=0.5, sort.firstwo.only=FALSE,
                                    order.by="decile", fancy.sort=FALSE,
                                    xtitle.suffix="Segments", ytitle=NULL,
                                    legend=F, cnv.marker.wex=0.035,
                                    parmar=c(1.1, 2.6, 0.4, 0.1)){
  # Subset segments based on pct.cutoff
  pct.colname <- paste(cohort, "dnm", csq, "excess_pct", sep="_")
  keep.idx <- which(segs[, pct.colname] >= pct.cutoff & segs$n_genes > 0)
  segs <- segs[keep.idx, ]

  # Helper function for ordering segments by decile
  order.by.decile <- function(res){
    deciles <- seq(1, 0, -0.1)
    ranks <- sapply(deciles, function(d){
      apply(res, 1, function(vals){min(which(vals>=d))})
    })
    order(ranks[, 1], ranks[, 2], ranks[, 3], ranks[, 4], ranks[, 5], ranks[, 6],
          ranks[, 7], ranks[, 8], ranks[, 9], ranks[, 10], ranks[, 11])
  }
  # Helper function for ordering segments by quintile
  order.by.quintile <- function(res){
    quintiles <- seq(1, 0, -0.2)
    ranks <- sapply(quintiles, function(d){
      apply(res, 1, function(vals){min(which(vals>=d))})
    })
    order(ranks[, 1], ranks[, 2], ranks[, 3], ranks[, 4], ranks[, 5])
  }
  # Helper function for fancy sorting
  fancy.order <- function(res){
    n <- ncol(res)
    order(apply(res, 1, function(vals){
      vals <- rev(rev(vals) - c(rev(vals)[-1], 0))
      vals <- vals / sum(vals)
      sum(vals * (n:1)^2)
    }), decreasing=T)
  }

  # Get plot data
  res <- calc.dnm.excess.cdf.matrix(segs, dnms, csq, n.max.genes)
  if(norm==T){
    res <- t(apply(res, 1, function(vals){vals/max(vals, na.rm=T)}))
    if(sort.firstwo.only==FALSE){
      res <- res[order(res[, 5], res[, 4], res[, 3], res[, 2], res[, 1], decreasing=T), ]
      if(fancy.sort){
        res <- res[fancy.order(res), ]
      }else{
        if(order.by=="decile"){
          res <- res[order.by.decile(res), ]
        }else if(order.by=="quintile"){
          res <- res[order.by.quintile(res), ]
        }
      }
    }else{
      res <- res[order(res[, 2], res[, 1], decreasing=T), ]
      if(fancy.sort){
        res <- res[fancy.order(res), ]
      }else{
        if(order.by=="decile"){
          res <- res[order.by.decile(res[, c(1, 2, ncol(res))]), ]
        }else if(order.by=="quintile"){
          res <- res[order.by.quintile(res[, c(1, 2, ncol(res))]), ]
        }
      }
    }
  }
  res <- res[which(!is.na(res[, 1])), ]
  colors <- rev(viridis(n.max.genes + 1))
  ymax <- max(res, na.rm=T)

  # Prep plot area
  par(bty="n", mar=parmar)
  plot(NA, xlim=c(-1, nrow(res)+1), ylim=c(-cnv.marker.wex*ymax, ymax),
       xaxt="n", yaxt="n", xlab="", ylab="", xaxs="i", yaxs="i")

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
    res.sums <- apply(res, 2, sum, na.rm=T)
    res.pct <- round(100 * rev(rev(res.sums) - c(rev(res.sums)[-1], 0)) / max(res.sums))
    legend.labs <- c(as.expression(bquote(.(paste("Top Gene (", res.pct[1], "%)", sep="")))),
                     as.expression(bquote(2^"nd" ~ "Gene" ~ .(paste("(", res.pct[2], "%)", sep="")))),
                     as.expression(bquote(3^"rd" ~ "Gene" ~ .(paste("(", res.pct[3], "%)", sep="")))),
                     sapply(4:n.max.genes, function(x){
                       as.expression(bquote(.(x)^"th" ~ "Gene" ~ .(paste("(", res.pct[x], "%)", sep=""))))
                       }),
                     as.expression(bquote("All Others" ~ .(paste("(", res.pct[n.max.genes+1], "%)", sep="")))))
    legend("topright", fill=colors, legend=legend.labs, bty="n", border=NA)
  }
}


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(zoo, quietly=T)
require(optparse, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog loci.bed segs.tsv ddd.tsv asc.tsv asc.control.tsv mutrates.tsv genes.bed out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 8){
  stop(paste("Eight positional arguments required: loci.bed, segs.tsv, ddd.tsv, asc.tsv, asc.control.tsv, mutrates.tsv, genes.bed, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
loci.in <- args$args[1]
segs.in <- args$args[2]
ddd.dnms.in <- args$args[3]
asc.dnms.in <- args$args[4]
asc.control.dnms.in <- args$args[5]
gene.mutrates.in <- args$args[6]
basic.gene.features.in <- args$args[7]
out.prefix <- args$args[8]

# # DEV PARAMETERS
# loci.in <- "~/scratch/rCNV.final_segments.loci.bed.gz"
# segs.in <- "~/scratch/rCNV2_analysis_d2.master_segments.bed.gz"
# ddd.dnms.in <- "~/scratch/ddd_dnm_counts.tsv.gz"
# asc.dnms.in <- "~/scratch/asc_dnm_counts.tsv.gz"
# asc.control.dnms.in <- "~/scratch/asc_dnm_counts.unaffecteds.tsv.gz"
# gene.mutrates.in <- "~/scratch/gene_mutation_rates.tsv.gz"
# basic.gene.features.in <- "~/scratch/gencode.v19.canonical.pext_filtered.genomic_features.bed.gz"
# out.prefix <- "~/scratch/test_oligogenicity"

# Load dnm counts and residualize by mutation rates
mutrates <- load.mutrates(gene.mutrates.in)
ddd <- load.dnms(ddd.dnms.in, mutrates)
asc <- load.dnms(asc.dnms.in, mutrates)
asc.control <- load.dnms(asc.control.dnms.in, mutrates)
ddd.plus.asc <- combine.dnms(ddd, asc, mutrates)
dnms <- list("DDD"=ddd, "ASC"=asc, "ASC_unaffected"=asc.control,
             "DDD_plus_ASC"=ddd.plus.asc)

# Load basic gene features for ordering by chromosome
gbed <- read.table(basic.gene.features.in, header=T, sep="\t", comment.char="")[, c(1, 4)]
colnames(gbed) <- c("chrom", "gene")

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

# Enumerate distribution of expected DNM excess by segment size
enumerated.excesses <- lapply(names(dnms), function(cohort){
  cres <- lapply(csqs, function(csq){
    enumerate.all.excess.dnms(ndd.segs, dnms[[cohort]], csq, gbed)
  })
  names(cres) <- csqs
  return(cres)
})
names(enumerated.excesses) <- names(dnms)

# Annotate each segment with its relative excess percentile
ndd.segs <- assign.excess.pct(segs, enumerated.excesses)

# Quantify deletion and duplication excesses and write to file
quantify.excesses(ndd.segs, paste(out.prefix, "quantified_dnm_excesses.tsv", sep="."))

# Scatterplots of observed vs. expected excess DNMs per segment
sapply(names(dnms), function(cohort){
  sapply(c("lof", "mis"), function(csq){
    cnvtype <- c("lof" = "DEL", "mis" = "DUP")[csq]
    pdf(paste(out.prefix, cohort, csq, "observed_vs_expected_excess.pdf", sep="."),
        height=2.75, width=3.3)
    plot.excess.dnms.vs.expected(ndd.segs[which(ndd.segs$cnv == cnvtype), ], label.top.n=12,
                                 enumerated.excesses, cohort, csq, min.label.quantile=0.5)
    dev.off()
  })
})

# Detailed distributions of proportional DNM excess for best segments
sapply(names(dnms), function(cohort){
  sapply(csqs, function(csq){
    if(csq=="lof"){
      cnv.types <- c("DEL")
    }else if(csq=="mis"){
      cnv.types <- c("DUP")
    }else{
      cnv.types <- c("DEL", "DUP")
    }
    pdf(paste(out.prefix, cohort, csq, "excess_dnm_distrib_bygene.best_loci.pdf", sep="."),
        height=3.75, width=5)
    dnm.excess.cdf.highlights(ndd.segs[which(ndd.segs$cnv %in% cnv.types), ],
                              dnms[[cohort]], cohort, csq, top.n.segs=12, top.n.genes=5)
    dev.off()
  })
})

# Collapsed distribution of excess DNMs across all NDD segments
pct.cutoff <- 0
sapply(names(dnms), function(cohort){
  sapply(csqs, function(csq){
    cnv.marker.wex <- 0
    if(csq=="lof"){
      cnv.types <- c("DEL")
      cnv.lab <- "DEL"
      ytitle.1 <- bquote("Excess" ~ italic("dn") * "PTVs")
      ytitle.2 <- bquote("Prop. of Excess" ~ italic("dn") * "PTVs")
      xtitle.suffix.2 <- paste('"DEL Segs. with Excess" ~ italic("dn") * "PTV"', sep="")
    }else if(csq=="mis"){
      cnv.types <- c("DUP")
      cnv.lab <- "DUP"
      ytitle.1 <- bquote("Excess" ~ italic("dn") * "Mis.")
      ytitle.2 <- bquote("Prop. of Excess" ~ italic("dn") * "Mis.")
      xtitle.suffix.2 <- paste('"DUP Segs. with Excess" ~ italic("dn") * "Mis."', sep="")
    }else{
      cnv.types <- c("DEL", "DUP")
      cnv.lab <- "rCNV"
      cnv.marker.wex <- 0.035
      ytitle.1 <- bquote("Excess" ~ italic("dn") * "Syn.")
      ytitle.2 <- bquote("Prop. of Excess" ~ italic("dn") * "Syn.")
      xtitle.suffix.2 <- paste('"Segs. with Excess" ~ italic("dn") * "Syn."', sep="")
    }
    pdf(paste(out.prefix,  cohort, csq, "excess_dnm_distrib_bygene.pdf", sep="."),
        height=2, width=3)
    dnm.excess.cdf.barplots(ndd.segs[which(ndd.segs$cnv %in% cnv.types), ],
                            dnms[[cohort]], cohort=cohort, csq=csq,
                            norm=F, legend=T, pct.cutoff=pct.cutoff,
                            cnv.marker.wex=cnv.marker.wex,
                            xtitle.suffix=paste('"NDD-Associated', cnv.lab, 'Segs."'),
                            ytitle=ytitle.1)
    dev.off()
    pdf(paste(out.prefix, cohort, csq, "excess_dnm_distrib_bygene.norm.pdf", sep="."),
        height=2, width=3)
    dnm.excess.cdf.barplots(ndd.segs[which(ndd.segs$cnv %in% cnv.types), ],
                            dnms[[cohort]], cohort=cohort, csq=csq,
                            norm=T, legend=F, pct.cutoff=pct.cutoff,
                            cnv.marker.wex=cnv.marker.wex, fancy.sort=TRUE,
                            xtitle.suffix=xtitle.suffix.2,
                            ytitle=ytitle.2)
    dev.off()
  })
})

