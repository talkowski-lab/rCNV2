#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020 Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Plot distributions of primary & secondary p-values from crb association meta analysis

# DEV NOTE: this code snippet is mostly saved for convenience, and references paths
#           to local files on RLC's computer that were crbrated by crb_burden_analysis.sh


options(scipen=1000, stringsAsFactors=F)


# Plot histogram of observed primary p-values from rCNV sliding crb meta-analysis

x <- read.table("~/scratch/rCNV.observed_pval_matrix.txt.gz", header=T, sep="\t")
x <- t(apply(x, 1, function(x){10^-as.numeric(x)}))

pdf("~/scratch/rCNV.observed_pval_distribution.crb_burden_meta_analysis.pdf",
    height=5, width=8)
hist(x[which(x<1)], breaks=200, col="gray20", freq=F,
     xlab="Meta-analysis P-value", 
     main=paste("Distribution of meta-analysis P-values",
                "across all phenotypes, CRBs, and CNV classes",
                sep="\n"))
mtext(3, line=-1.5, text=paste("Note: CRBs with P=1 not shown (i.e., where zero case CNVs are observed).", 
                               "These comprise the majority of all crbs.", sep="\n"),
      font=3, cex=0.85)
dev.off()


# Plot histogram of observed secondary p-values from rCNV sliding crb meta-analysis

s <- read.table("~/scratch/rCNV.observed_pval_matrix_secondary.txt.gz", header=T, sep="\t")
s <- t(apply(s, 1, function(s){10^-as.numeric(s)}))

pdf("~/scratch/rCNV.observed_secondary_pval_distribution.crb_burden_meta_analysis.pdf",
    height=5, width=8)
hist(s[which(s<1)], breaks=200, col="gray20", freq=F,
     xlab="Secondary meta-analysis P-value", 
     main=paste("Distribution of secondary meta-analysis P-values",
                "across all phenotypes, CRBs, and CNV classes",
                sep="\n"))
mtext(3, line=-1.5, text=paste("Note: CRBs with P=1 not shown (i.e., where zero case CNVs are observed).", 
                               "These comprise the majority of all crbs.", sep="\n"),
      font=3, cex=0.85)
dev.off()


# Scatterplot of all primary vs secondary P-values

pairs <- data.frame("p"=as.vector(unlist(x)),
                    "s"=as.vector(unlist(s)))
pairs$p[which(is.na(pairs$p))] <- 1
pairs$s[which(is.na(pairs$s))] <- 1
pairs <- pairs[-which(pairs$p>=0.05 & pairs$s>=0.05), ]
sig.pairs.idx <- which(-log10(pairs$p)>6 & pairs$s<0.05)

png("~/scratch/rCNV.crb_burden.primary_vs_secondary_pvalues.png", height=300*5, width=300*5, res=300)
plot(-log10(pairs$p), -log10(pairs$s), cex=0.2,
     panel.first=c(abline(0, 1, lty=2, col="gray50")),
     xlim=c(0, 20), ylim=c(0, 20),
     xlab="-log10(P) Primary", ylab="-log10(P) Secondary",
     main="Primary vs. Secondary P-Values\nAll phenotypes & all CNV types")
points(-log10(pairs$p)[sig.pairs.idx], -log10(pairs$s)[sig.pairs.idx], 
       pch=19, cex=0.2, col="red")
rect(xleft=0, xright=-log10(0.05), ybottom=0, ytop=-log10(0.05), col="black")
dev.off()

