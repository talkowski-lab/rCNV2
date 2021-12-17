#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2020-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Determine gene features predictive of dosage sensitivity with regression


options(stringsAsFactors=F, scipen=1000)


######################
### DATA FUNCTIONS ###
######################
# Perform linear elastic net regression of gene features versus a continuous response
gene.glm <- function(scores, feats, response.name, glm.alpha=NULL, seed=2020){
  glm.df <- feats[which(feats$gene %in% scores$gene), ]
  if(!(response.name %in% colnames(scores))){
    stop(paste("Error: response.name \"", response.name, "\" not found in scores."))
  }
  Y <- as.numeric(sapply(as.vector(glm.df$gene), function(g){
    scores[which(scores$gene==g), which(colnames(scores)==response.name)]
  }))
  Y <- as.vector(as.numeric(scale(Y)))
  glm.df <- apply(glm.df[, -1], 2, function(vals){as.vector(as.numeric(unlist(vals)))})

  # Use caret::train to determine the optimal alpha if none is provided
  if(is.null(glm.alpha)){
    glm.df <- as.data.frame(cbind(Y, glm.df))
    set.seed(seed)
    model <- caret::train(Y ~ ., data=glm.df, method="glmnet",
                          trControl=trainControl("cv", number=10),
                          tuneLength=10)
    fit <- model$bestTune
    print(fit)
    coefs <- coef(model$finalModel, model$bestTune$lambda)
  }else{
    cv <- cv.glmnet(glm.df, Y, alpha=glm.alpha, nfolds=10)
    fit <- glmnet(glm.df, Y, alpha=glm.alpha, lambda=cv$lambda.1se)
    coefs <- coef(fit)
  }
  coefs <- coefs[order(coefs[, 1]), ]
  return(coefs)
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Barplot of top N regression coefficients
plot.top.coefs <- function(coefs, feat.meta, colors, borders,
                           lab.subscripts, N=10, parmar=c(0.1, 0.5, 3, 0.5)){
  # Get plotting data
  coefs <- coefs[which(names(coefs) != "(Intercept)")]
  top.coefs <- rev(sort(head(coefs[order(-abs(coefs))], N)))
  bar.colors <- sapply(top.coefs, function(val){
    if(val < 0){colors[1]}else if(val >= 0){colors[2]}
  })
  bar.borders <- sapply(top.coefs, function(val){
    if(val < 0){borders[1]}else if(val >= 0){borders[2]}
  })
  lab.pos <- sapply(top.coefs, function(val){
    if(val < 0){4}else if(val >= 0){2}
  })
  lab.names <- sapply(names(top.coefs), function(fname){
    feat.meta$name[which(feat.meta$feature==fname)]
  })
  xlims <- c(-1, 1)*max(abs(top.coefs), na.rm=T)

  # Prep plot area
  par(bty="n", mar=parmar)
  plot(NA, xlim=xlims, ylim=c(N, 0),
       xaxt="n", yaxt="n", xlab="", ylab="", yaxs="i")

  # Add gridlines
  x.ax.at <- axTicks(3)
  abline(v=x.ax.at, col="white")

  # Add rectangles
  rect(xleft=rep(0, N), xright=top.coefs,
       ybottom=(1:N)-0.85, ytop=(1:N)-0.15,
       bty="o", border=bar.borders, col=bar.colors)
  text(x=rep(0, N), y=(1:N)-0.5, pos=lab.pos, labels=lab.names, xpd=T, cex=5/6)

  # Add top axis
  axis(3, at=c(-100, 100), labels=NA, tck=0, col=blueblack)
  axis(3, at=x.ax.at, labels=NA, col=blueblack, tck=-0.025)

  if(length(lab.subscripts) == 2){
    axis(3, at=x.ax.at, labels=abs(x.ax.at), tick=F, line=-0.65)
    axis(3, at=mean(c(0, par("usr")[1])), tick=F, line=0.4,
         labels=bquote(hat(beta)[.(lab.subscripts[1])]))
    axis(3, at=mean(c(0, par("usr")[2])), tick=F, line=0.4,
         labels=bquote(hat(beta)[.(lab.subscripts[2])]))
  }else{
    axis(3, at=x.ax.at, labels=x.ax.at, tick=F, line=-0.65)
    axis(3, at=0, tick=F, line=0.4,
         labels=bquote(hat(beta)[.(lab.subscripts[1])]))
  }
  abline(v=0, col=blueblack, lwd=2)
}


#####################
### RSCRIPT BLOCK ###
#####################
require(rCNV2, quietly=T)
require(optparse, quietly=T)
require(glmnet, quietly=T)
require(caret, quietly=T)

# List of command-line options
option_list <- list()

# Get command-line arguments & options
args <- parse_args(OptionParser(usage=paste("%prog scores.tsv features.bed feature.meta.tsv out_prefix", sep=" "),
                                option_list=option_list),
                   positional_arguments=TRUE)
opts <- args$options

# Checks for appropriate positional arguments
if(length(args$args) != 4){
  stop(paste("Four positional arguments required: scores.tsv, features.bed, feature.meta.tsv, and output_prefix\n", sep=" "))
}

# Writes args & opts to vars
scores.in <- args$args[1]
features.in <- args$args[2]
feature.meta.in <- args$args[3]
out.prefix <- args$args[4]

# # DEV PARAMETERS
# scores.in <- "~/scratch/rCNV.gene_scores.tsv.gz"
# features.in <- "~/scratch/gencode.v19.canonical.pext_filtered.all_features.no_variation.transformed.bed.gz"
# feature.meta.in <- "~/scratch/gene_feature_metadata.tsv"
# out.prefix <- "~/scratch/test_gene_score_feature_regressions"

# Load scores & classify genes into subgroups based on scores
scores <- load.scores(scores.in)
ds.groups <- classify.genes.by.score(scores, hc.cutoff=0.9, lc.cutoff=0.5)

# Compute regression gradients
scores$ds.gradient <- apply(scores[, -1], 1, function(vals){
  min(as.numeric(vals), na.rm=T)
})
scores$hits.gradient <- scores$pTriplo - scores$pHaplo

# Load & normalize features
feats <- load.features(features.in, fill="mean", norm=T)
feats <- mod.features(feats)
feat.meta <- load.gene.feature.metadata(feature.meta.in)

# Compute regression coefficients
ds.coefs <- gene.glm(scores, feats, "ds.gradient", glm.alpha=NULL)
hits.coefs <- gene.glm(scores[which(scores$ds.gradient >= 0.5), ], feats, "hits.gradient", glm.alpha=NULL)
# hits.coefs.all <- gene.glm(scores, feats, "hits.gradient", glm.alpha=NULL)

# Plot top N features for ds and hi/ts gradient coefficients
coefs.to.plot <- 16
pdf(paste(out.prefix, "gradient_regression.top_coefs.ds.pdf", sep="."),
    height=3.7, width=4)
plot.top.coefs(ds.coefs, feat.meta, colors=c(ns.color, cnv.colors[3]),
               borders=c("gray30", cnv.blacks[3]), lab.subscripts=c("min(pHaplo, pTriplo)"),
               N=coefs.to.plot)
dev.off()
pdf(paste(out.prefix, "gradient_regression.top_coefs.hits.pdf", sep="."),
    height=3.7, width=4)
plot.top.coefs(hits.coefs, feat.meta, colors=cnv.colors[1:2],
               borders=cnv.blacks[1:2], N=coefs.to.plot, lab.subscripts=c("pTriplo-pHaplo"))
dev.off()

