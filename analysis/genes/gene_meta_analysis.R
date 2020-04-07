#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2019 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Perform meta-analysis of gene-based CNV associations


options(scipen=1000, stringsAsFactors=F)


############
###FUNCTIONS
############

# Compute basic OR with adjustment for zero-inflation
calc.or <- function(control_ref, control_alt, case_ref, case_alt, adj=0.5){
  case.odds <- (case_alt + adj) / (case_ref + adj)
  control.odds <- (control_alt + adj) / (control_ref + adj)
  case.odds / control.odds
}


# Read an input file of association statistics
read.stats <- function(stats.in, prefix, p.is.phred){
  # Read data & subset to necessary columns
  stats <- read.table(stats.in, header=T, sep="\t", comment.char="")
  colnames(stats)[1] <- "chr"
  cols.to.keep <- c("chr", "start", "end", "gene", "case_alt", "case_ref",
                    "control_alt", "control_ref", "fisher_phred_p")
  stats <- stats[, which(colnames(stats) %in% cols.to.keep)]
  colnames(stats)[which(colnames(stats)=="fisher_phred_p")] <- "p_value"
  if(p.is.phred==T){
    stats$p_value <- 10^-stats$p_value
  }
  stats$odds_ratio <- calc.or(stats$control_ref, stats$control_alt,
                              stats$case_ref, stats$case_alt)
  colnames(stats)[-(1:4)] <- paste(prefix, colnames(stats)[-(1:4)], sep=".")
  return(stats)
}


# Single density-colored scatterplot
dens.scatter <- function(x, y, pt.cex=1, parmar=rep(0.5, 4), add.cor=T){
  # Based on heatscatter.R from Colby Chiang
  # (https://github.com/cc2qe/voir/blob/master/bin/heatscatter.R)
  par(mar=parmar, bty="o")
  plot.df <- data.frame("x"=x, "y"=y)
  plot.df <- plot.df[which(!is.infinite(plot.df$x) & !is.infinite(plot.df$y)
                           & !is.na(plot.df$x) & !is.na(plot.df$y)), ]
  dens <- densCols(plot.df$x, plot.df$y, colramp=colorRampPalette(c("black", "white")))
  plot.df$dens <- col2rgb(dens)[1, ] + 1L
  palette <- colorRampPalette(c("#440154", "#33638d", "#218f8d",
                                "#56c667", "#fde725"))(256)
  plot.df$col <- palette[plot.df$dens]
  plot.df <- plot.df[order(plot.df$dens), ]
  plot(plot.df$x, plot.df$y, type="n",
       xlab="", xaxt="n", ylab="", yaxt="n")
  abline(lm(y ~ x, data=plot.df), lwd=2)
  points(plot.df$x, plot.df$y, col=plot.df$col, pch=19, cex=pt.cex)
  if(add.cor==T){
    r <- cor(plot.df$x, plot.df$y)
    text(x=par("usr")[2], y=par("usr")[3]+(0.075*(par("usr")[4]-par("usr")[3])),
         labels=bquote(italic(R) == .(formatC(round(r, 3), digits=3))), pos=2)
  }
}


# Scatterplot grid of log odds ratios between cohorts
or.corplot.grid <- function(stats.list, pt.cex=1){
  ncohorts <- length(stats.list)
  cohorts <- names(stats.list)
  par(mfrow=c(ncohorts, ncohorts))
  ymar <- 3
  xmar <- 3.7
  sapply(1:ncohorts, function(r){
    sapply(1:ncohorts, function(c){
      parmar <- c(ymar/2, xmar/2, ymar/2, xmar/2)
      # Set margins
      if(c==1){
        parmar[2] <- xmar-0.2; parmar[4] <- 0.2
      }
      if(c==ncohorts){
        parmar[2] <- 0.2; parmar[4] <- xmar-0.2
      }
      if(r==1){
        parmar[1] <- 0.2; parmar[3] <- ymar-0.2
      }
      if(r==ncohorts){
        parmar[1] <- ymar-0.2; parmar[3] <- 0.2
      }
      # Don't plot self-self correlations
      if(c==r){
        par(mar=parmar)
        plot(x=c(0, 1), y=c(0, 1), type="n",
             xaxt="n", xlab="", yaxt="n", ylab="")
        rect(xleft=par("usr")[1], xright=par("usr")[2],
             ybottom=par("usr")[3], ytop=par("usr")[4],
             bty="n", col="gray95")
        text(x=0.5, y=0.5, labels=bquote(italic(R)==1))
        # Otherwise, plot as normal
      }else{
        dens.scatter(x=log10(stats.list[[c]][, grep("odds_ratio", colnames(stats.list[[c]]), fixed=T)]), 
                     y=log10(stats.list[[r]][, grep("odds_ratio", colnames(stats.list[[r]]), fixed=T)]),
                     parmar=parmar, pt.cex=pt.cex)
      }
      # Add headers & axes
      if(c==1){
        mtext(2, line=0.2, text=cohorts[r], font=2)
      }
      if(r==1){
        mtext(3, line=0.2, text=cohorts[c], font=2)
      }
      if(c==ncohorts){
        axis(4, labels=NA)
        axis(4, tick=F, line=-0.5, las=2)
        mtext(4, line=2.25, text=bquote(log[10](OR)), cex=0.8)
      }
      if(r==ncohorts){
        axis(1, labels=NA)
        axis(1, tick=F, line=-0.5)
        mtext(1, line=2, text=bquote(log[10](OR)) ,cex=0.8)
      }
    })
  })
}


# Merge a list of association statistics
combine.stats <- function(stats.list){
  # Merge all cohorts
  merged <- stats.list[[1]]
  for(i in 2:length(stats.list)){
    merged <- merge(merged, stats.list[[i]], 
                    by=c("chr", "start", "end", "gene"),
                    all=F, sort=F)
  }
  merged[, -c(1:4)] <- apply(merged[, -c(1:4)], 2, as.numeric)
  
  # Count number of nominally significant individual cohorts
  n_nom_sig <- apply(merged[, grep(".p_value", colnames(merged), fixed=T)],
                     1, function(pvals){length(which(pvals<=0.05))})
  merged$n_nominal_cohorts <- n_nom_sig
  
  # Determine most significant cohort per row
  merged$top_cohort <- unlist(apply(merged[, grep(".p_value", colnames(merged), fixed=T)], 1, 
                                    function(pvals){
                                      head(names(stats.list)[which(pvals == min(pvals, na.rm=T))], 1)
                                    }))
  merged <- merged[, -grep(".p_value", colnames(merged), fixed=T)]
  return(merged)
}


# Apply empirical continuity correction to meta-analysis data frame
# Per Sweeting et al., Stat. Med., 2004 (section 3.3)
sweeting.correction <- function(meta.df, cc.sum=0.01){
  # Count number of carriers & non-carriers
  n.alt <- sum(meta.df[, grep("_alt", colnames(meta.df), fixed=T)])
  n.ref <- sum(meta.df[, grep("_ref", colnames(meta.df), fixed=T)])
  # Require at least one CNV to be observed
  if(n.alt>0){
    nt <- n.alt
    R <- n.ref/n.alt
    # Pooled odds ratio estimate of all non-zero studies with at least one case sample
    nonzero.studies <- intersect(which(apply(meta.df[, grep("_alt", colnames(meta.df), fixed=T)], 1, sum)>0),
                                 which(apply(meta.df[, grep("case_", colnames(meta.df), fixed=T)], 1, sum)>0))
    nonzero.case.odds <- sum(meta.df$case_alt[nonzero.studies])/sum(meta.df$case_ref[nonzero.studies])
    nonzero.control.odds <- sum(meta.df$control_alt[nonzero.studies])/sum(meta.df$control_ref[nonzero.studies])
    if(!is.nan(nonzero.case.odds) & !is.nan(nonzero.control.odds)){
      if(nonzero.control.odds>0){
        ohat <- nonzero.case.odds/nonzero.control.odds
        # Otherwise, apply standard continuity correction of 0.5 to pooled estimate if no CNVs observed in controls
      }else{
        nonzero.case.odds <- (sum(meta.df$case_alt[nonzero.studies])+0.5)/(sum(meta.df$case_ref[nonzero.studies])+0.5)
        nonzero.control.odds <- (sum(meta.df$control_alt[nonzero.studies])+0.5)/(sum(meta.df$control_ref[nonzero.studies])+0.5)
        ohat <- nonzero.case.odds/nonzero.control.odds
      }
      # Otherwise, apply standard continuity correction to pooled estimate of *all* studies
    }else{
      nonzero.case.odds <- (sum(meta.df$case_alt)+0.5)/(sum(meta.df$case_ref)+0.5)
      nonzero.control.odds <- (sum(meta.df$control_alt)+0.5)/(sum(meta.df$control_ref)+0.5)
      ohat <- nonzero.case.odds/nonzero.control.odds
    }
    
    # Solve for kc & kt
    kc <- R/(R+ohat)
    kt <- ohat/(R+ohat)
    # Compute continuity corrections
    cor.case_alt <- cc.sum * kt
    cor.case_ref <- cc.sum * kc
    cor.control_alt <- cc.sum * (nt + kt)
    cor.control_ref <- cc.sum * ((nt*R) + kc)
    # Apply continuity corrections
    meta.df$case_alt <- meta.df$case_alt + cor.case_alt
    meta.df$case_ref <- meta.df$case_ref + cor.case_ref
    meta.df$control_alt <- meta.df$control_alt + cor.control_alt
    meta.df$control_ref <- meta.df$control_ref + cor.control_ref
  }
  return(meta.df)
}


# Make meta-analysis data frame for a single gene
make.meta.df <- function(stats.merged, cohorts, row.idx, empirical.continuity=T, drop_top_cohort=F){
  if(drop_top_cohort==T){
    top_cohort <- stats.merged$top_cohort[row.idx]
    cohorts <- cohorts[which(cohorts != top_cohort)]
    stats.merged <- stats.merged[, grep(top_cohort, colnames(stats.merged), fixed=T, invert=T)]
  }
  ncohorts <- length(cohorts)
  meta.df <- data.frame("cohort"=1:ncohorts,
                        "control_ref"=as.numeric(stats.merged[row.idx, grep("control_ref", colnames(stats.merged), fixed=T)]),
                        "case_ref"=as.numeric(stats.merged[row.idx, grep("case_ref", colnames(stats.merged), fixed=T)]),
                        "control_alt"=as.numeric(stats.merged[row.idx, grep("control_alt", colnames(stats.merged), fixed=T)]),
                        "case_alt"=as.numeric(stats.merged[row.idx, grep("case_alt", colnames(stats.merged), fixed=T)]),
                        "cohort_name"=cohorts)
  if(empirical.continuity==T){
    meta.df <- sweeting.correction(meta.df)
  }
  return(meta.df)
}


# Perform meta-analysis for a single gene
meta.single <- function(stats.merged, cohorts, row.idx, model="fe", empirical.continuity=T, drop_top_cohort=F){
  # If no CNVs are observed, return all NAs
  if(sum(stats.merged[row.idx, grep("_alt", colnames(stats.merged), fixed=T)])>0){
    meta.df <- make.meta.df(stats.merged, cohorts, row.idx, empirical.continuity, drop_top_cohort)
    # If strictly zero case CNVs are observed, unable to estimate effect size
    if(all(meta.df$case_alt==0)){
      out.v <- c(rep(NA, 4), 0)
    }else{
      # Meta-analysis
      if(model=="re"){
        meta.res <- tryCatch(rma.uni(ai=control_ref, bi=case_ref, ci=control_alt, di=case_alt,
                                     measure="OR", data=meta.df, method="REML", random = ~ 1 | cohort, slab=cohort_name,
                                     add=0, drop00=F, correct=F, digits=5, control=list(maxiter=100, stepadj=0.5)),
                             error=function(e){
                               print(paste("row", row.idx, "failed to converge. Retrying with more iterations...", sep=" "))
                               rma.uni(ai=control_ref, bi=case_ref, ci=control_alt, di=case_alt,
                                       measure="OR", data=meta.df, method="REML", random = ~ 1 | cohort, slab=cohort_name,
                                       add=0, drop00=F, correct=F, digits=5, control=list(maxiter=10000, stepadj=0.4))
                             })
        out.v <- as.numeric(c(meta.res$b[1,1], meta.res$ci.lb, meta.res$ci.ub,
                              meta.res$zval, -log10(meta.res$pval)))
      }else if(model=="mh"){
        meta.res <- rma.mh(ai=control_ref, bi=case_ref, ci=control_alt, di=case_alt,
                           measure="OR", data=meta.df, slab=cohort_name,
                           add=0, drop00=F, correct=F)
        out.v <- as.numeric(c(meta.res$b, meta.res$ci.lb, meta.res$ci.ub,
                              meta.res$zval, -log10(meta.res$MHp)))
      }else if(model=="fe"){
        meta.res <- tryCatch(rma.uni(ai=control_ref, bi=case_ref, ci=control_alt, di=case_alt,
                                     measure="OR", data=meta.df, method="FE", slab=cohort_name,
                                     add=0, drop00=F, correct=F, digits=5, control=list(maxiter=100, stepadj=0.5)),
                             error=function(e){
                               print(paste("row", row.idx, "failed to converge. Retrying with more iterations...", sep=" "))
                               rma.uni(ai=control_ref, bi=case_ref, ci=control_alt, di=case_alt,
                                       measure="OR", data=meta.df, method="FE", slab=cohort_name,
                                       add=0, drop00=F, correct=F, digits=5, control=list(maxiter=10000, stepadj=0.4))
                             })
        out.v <- as.numeric(c(meta.res$b[1,1], meta.res$ci.lb, meta.res$ci.ub,
                              meta.res$zval, -log10(meta.res$pval)))
      }
      # Force to p-values reflecting Ha : OR > 1
      if(!is.na(out.v[1]) & !is.na(out.v[5])){
        if(out.v[1] < 0){
          out.v[5] <- 0
        }
      }
    }
    return(out.v)
  }else{
    rep(NA, 5)
  }
}


# Make meta-analysis lookup table to shorten time required to run full meta-analysis
make.meta.lookup.table <- function(stats.merged, cohorts, model, empirical.continuity=T){
  unique.counts.df <- unique(stats.merged[, sort(unique(c(grep("_ref", colnames(stats.merged), fixed=T),
                                                          grep("_alt", colnames(stats.merged), fixed=T))))])
  
  unique.stats <- t(sapply(1:nrow(unique.counts.df), function(i){
    meta.single(unique.counts.df, cohorts, i, model, empirical.continuity)
  }))
  
  lookup.table <- cbind(unique.counts.df, unique.stats)
  stat.colnames <- c("meta_lnOR", "meta_lnOR_lower", "meta_lnOR_upper", "meta_z", "meta_phred_p")
  if(model=="mh"){
    stat.colnames[length(stat.colnames)-1] <- "CMH_chisq"
  }
  colnames(lookup.table)[(ncol(lookup.table)-4):ncol(lookup.table)] <- stat.colnames
  
  return(lookup.table)
}


# Apply saddlepoint approximation to vector of Z-scores to generate adjusted P-values
saddlepoint.adj <- function(zscores, phred=T){
  mu.hat <- mean(zscores, na.rm=T)
  sd.hat <- sd(zscores, na.rm=T)
  cumuls <- gaussianCumulants(mu.hat, sd.hat)
  dx <- 0.01
  x <- seq(-40, 40, dx)
  saddle.pdf <- saddlepoint(x, 1, cumuls)$approx
  saddle.cdf <- cumsum(saddle.pdf * 0.01)
  calc.saddle.p <- function(z){if(!is.na(z)){1 - tail(saddle.cdf[which(x<z)], 1)}else{NA}}
  new.pvals <- sapply(zscores, calc.saddle.p)
  if(phred==T){
    return(-log10(new.pvals))
  }else{
    return(new.pvals)
  }
}


# Wrapper function to perform a meta-analysis on all genes
meta <- function(stats.merged, cohorts, model="fe", saddle=T, secondary=T){
  # Make meta-analysis lookup table
  meta.lookup.table <- make.meta.lookup.table(stats.merged, cohorts, model, 
                                              empirical.continuity=T)
  
  # Merge stats into full list
  meta.res <- merge(stats.merged, meta.lookup.table, sort=F, all.x=T, all.y=F)
  meta.res <- meta.res[with(meta.res, order(chr, start, gene)), ]
  
  # Adjust P-values using saddlepoint approximation of null distribution, if optioned
  if(saddle==T){
    meta.res$meta_phred_p <- saddlepoint.adj(meta.res$meta_z)
  }
  
  # Compute secondary P-value
  if(secondary==T){
    meta.res.secondary <- as.data.frame(t(sapply(1:nrow(meta.res), function(i){
      meta.single(meta.res, cohorts, i, model, empirical.continuity=T, drop_top_cohort=T)
    })))
    colnames(meta.res.secondary) <- c("meta_lnOR", "meta_lnOR_lower", "meta_lnOR_upper", "meta_z", "meta_phred_p")

    # Saddlepoint on secondary, if optioned
    if(saddle==T){
      meta.res.secondary$meta_phred_p <- saddlepoint.adj(meta.res.secondary$meta_z)
    }
    
    meta.res$meta_lnOR_secondary <- meta.res.secondary$meta_lnOR
    meta.res$meta_lnOR_lower_secondary <- meta.res.secondary$meta_lnOR_lower
    meta.res$meta_lnOR_upper_secondary <- meta.res.secondary$meta_lnOR_upper
    meta.res$meta_z_secondary <- meta.res.secondary$meta_z
    meta.res$meta_phred_p_secondary <- meta.res.secondary$meta_phred_p
  }
  
  # Compute pooled carrier frequencies
  case.alt <- apply(meta.res[, grep("case_alt", colnames(meta.res), fixed=T)], 1, sum, na.rm=T)
  case.ref <- apply(meta.res[, grep("case_ref", colnames(meta.res), fixed=T)], 1, sum, na.rm=T)
  meta.res$case_freq <- case.alt / (case.alt + case.ref)
  control.alt <- apply(meta.res[, grep("control_alt", colnames(meta.res), fixed=T)], 1, sum, na.rm=T)
  control.ref <- apply(meta.res[, grep("control_ref", colnames(meta.res), fixed=T)], 1, sum, na.rm=T)
  meta.res$control_freq <- control.alt / (control.alt + control.ref)
  
  # Format output
  return(as.data.frame(cbind(meta.res[, which(colnames(meta.res) %in% c("chr", "start", "end", "gene",
                                                                        "n_nominal_cohorts", "top_cohort",
                                                                        "case_freq", "control_freq"))],
                             meta.res[, grep("meta_", colnames(meta.res), fixed=T)])))
}


################
###RSCRIPT BLOCK
################

# Load required libraries
require(optparse, quietly=T)
require(metafor, quietly=T)
require(EQL, quietly=T)

# List of command-line options
option_list <- list(
  make_option(c("--or-corplot"), type="character", default=NULL, 
              help="output .jpg file for pairwise odds ratio correlation plot [default %default]",
              metavar="path"),
  make_option(c("--model"), type="character", default="fe", 
              help="specify meta-analysis model ('re': random effects, 'fe': fixed effects, 'mh': Mantel-Haenszel) [default '%default']",
              metavar="string"),
  make_option(c("--p-is-phred"), action="store_true", default=FALSE, 
              help="provided P-values are Phred-scaled (-log10(P)) [default %default]"),
  make_option(c("--spa"), action="store_true", default=FALSE, 
              help="apply saddlepoint approximation of null distribution [default %default]"),
  make_option(c("--no-secondary"), action="store_true", default=FALSE, 
              help="do not compute secondary P-value [default %default]")
)

# Get command-line arguments & options
args <- parse_args(OptionParser(usage="%prog infile outfile",
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Writes args & opts to variable
infile <- args$args[1]
outfile <- args$args[2]
corplot.out <- opts$`or-corplot`
model <- opts$model
p.is.phred <- opts$`p-is-phred`
spa <- opts$`spa`
secondary <- !(opts$`no-secondary`)

# # Dev parameters
# infile <- "~/scratch/dec18_gene_burden/gene.meta_test.input.txt2"
# outfile <- "~/scratch/dec18_gene_burden/gene.meta_test.results.bed"
# corplot.out <- "~/scratch/dec18_gene_burden/gene_corplot.test.jpg"
# model <- "fe"
# p.is.phred <- T
# spa <- T
# secondary <- T

# Read list of cohorts to meta-analyze
cohort.info <- read.table(infile, header=F, sep="\t")
ncohorts <- nrow(cohort.info)
stats.list <- lapply(1:ncohorts, function(i){read.stats(cohort.info[i, 2], 
                                                        cohort.info[i, 1],
                                                        p.is.phred)})
names(stats.list) <- cohort.info[, 1]

# Plot correlations of odds ratios between cohorts, if optioned
if(!is.null(corplot.out)){
  jpeg(corplot.out, res=300,
       height=300*(3.5+(ncohorts/2)),
       width=300*(4+(ncohorts/2)))
  or.corplot.grid(stats.list, pt.cex=0.5)
  dev.off()
}

# Conduct meta-analysis & write to file
stats.merged <- combine.stats(stats.list)
stats.meta <- meta(stats.merged, cohort.info[, 1], model=model, 
                   saddle=spa, secondary=secondary)
colnames(stats.meta)[1] <- "#chr"
write.table(stats.meta, outfile, sep="\t",
            row.names=F, col.names=T, quote=F)