#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Reformat de novo mutation data from Fu et al., 2021 (ASC exome manuscript)


require(rCNV2, quietly=T)
options(stringsAsFactors=FALSE)


# Read all data
args <- commandArgs(trailingOnly=TRUE)
asc.df <- read.table(args[1], header=T, sep="\t")
spark.df <- read.table(args[2], header=T, sep="\t")
elig.genes <- unique(sort(as.character(read.table(args[3], header=F)[, 1])))


# Collapse ASC & SPARK data for probands
pro.df <- aggregate(cbind(dn.ptv, dn.misb, dn.misa, dn.syn) ~ gene,
                    rbind(asc.df, spark.df), sum)
pro.df$dn.mis <- pro.df$dn.misb + pro.df$dn.misa
pro.df[, c("dn.misb", "dn.misa")] <- NULL
colnames(pro.df) <- c("#gene", "lof", "syn", "mis")


# Collapse ASC & SPARK data for probands
sib.df <- aggregate(cbind(dn.ptv.sib, dn.misb.sib, dn.misa.sib, dn.syn.sib) ~ gene,
                    rbind(asc.df, spark.df), sum)
sib.df$dn.mis <- sib.df$dn.misb + sib.df$dn.misa
sib.df[, c("dn.misb", "dn.misa")] <- NULL
colnames(sib.df) <- c("#gene", "lof", "syn", "mis")


# Write out reformatted data
write.table(pro.df[which(pro.df$`#gene` %in% elig.genes),
                   c("#gene", "lof", "mis", "syn")],
            "fu_asc_spark_dnm_counts.tsv",
            col.names=T, row.names=F, sep="\t", quote=F)
system("gzip -f fu_asc_spark_dnm_counts.tsv", wait=T, intern=F)
write.table(sib.df[which(sib.df$`#gene` %in% elig.genes),
                   c("#gene", "lof", "mis", "syn")],
            "fu_asc_spark_dnm_counts.unaffecteds.tsv",
            col.names=T, row.names=F, sep="\t", quote=F)
system("gzip -f fu_asc_spark_dnm_counts.unaffecteds.tsv", wait=T, intern=F)
