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

# Compute correlation coefficient between a gene score and all features
gene.cor <- function(scores, feats, response.name, method="spearman"){
  cor.df <- merge(scores[, c("gene", response.name)], feats, by="gene",
                  sort=F, all=F)
  coefs <- cor(cor.df[, which(colnames(cor.df) != "gene")])[, response.name]
  return(coefs[which(names(coefs) != response.name)])
}


##########################
### PLOTTING FUNCTIONS ###
##########################
# Barplot of top N regression coefficients
plot.top.coefs <- function(coefs, feat.meta, colors, borders,
                           x.title, parse.x.title=T, N=10, lab.cex=5/6,
                           parmar=c(0.1, 1, 3.2, 1)){
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
    feat.meta$name[which(feat.meta$feature==gsub("`", "", fname))]
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
  text(x=rep(0, N), y=(1:N)-0.5, pos=lab.pos, labels=lab.names, xpd=T, cex=lab.cex)

  # Add top axis
  axis(3, at=c(-100, 100), labels=NA, tck=0, col=blueblack)
  axis(3, at=x.ax.at, labels=NA, col=blueblack, tck=-0.025)

  if(length(x.title) == 2){
    axis(3, at=x.ax.at, labels=abs(x.ax.at), tick=F, line=-0.65)
    axis(3, at=mean(c(0, par("usr")[1])), tick=F, line=0.4,
         labels=if(parse.x.title){parse(text=x.title[1])}else{x.title[1]})
    axis(3, at=mean(c(0, par("usr")[2])), tick=F, line=0.4,
         labels=if(parse.x.title){parse(text=x.title[2])}else{x.title[2]})
  }else{
    axis(3, at=x.ax.at, labels=x.ax.at, tick=F, line=-0.65)
    axis(3, at=0, tick=F, line=0.4,
         labels=if(parse.x.title){parse(text=x.title[1])}else{x.title[1]})
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

# Load scores
scores <- load.scores(scores.in)

# Compute regression gradients
scores$ds.gradient <- apply(scores[, -1], 1, function(vals){
  min(as.numeric(vals), na.rm=T)
})
scores$hits.gradient <- scores$pTriplo - scores$pHaplo

# Load & normalize features
feats <- load.features(features.in, fill="mean", norm=T)
feats <- mod.features(feats)
feat.meta <- load.gene.feature.metadata(feature.meta.in)

# Compute regression coefficients including constraint
cat("DS Gradient (including constraint features):\n")
ds.coefs <- gene.glm(scores, feats, "ds.gradient", glm.alpha=NULL)
ds.cor <- gene.cor(scores, feats, "ds.gradient")
cat("TS-HI Gradient (including constraint features):\n")
hits.coefs <- gene.glm(scores[which(scores$ds.gradient >= 0.5), ],
                       feats, "hits.gradient", glm.alpha=NULL)
hits.cor <- gene.cor(scores, feats, "hits.gradient")

# Compute regression coefficients *without* constraint
constr.feat.names <- c(feat.meta$feature[which(feat.meta$category == "constraint")], "episcore")
feats.noConstr <- feats[, which(!(colnames(feats) %in% constr.feat.names))]
cat("DS Gradient (excluding constraint features):\n")
ds.coefs.noConstr <- gene.glm(scores, feats.noConstr, "ds.gradient", glm.alpha=NULL)
ds.cor.noConstr <- gene.cor(scores, feats.noConstr, "ds.gradient")
cat("TS-HI Gradient (excluding constraint features):\n")
hits.coefs.noConstr <- gene.glm(scores[which(scores$ds.gradient >= 0.5), ],
                                feats.noConstr, "hits.gradient", glm.alpha=NULL)
hits.cor.noConstr <- gene.cor(scores, feats.noConstr, "hits.gradient")

# Plot top N features for ds and hi/ts gradient coefficients
coefs.to.plot <- 16
pdf.dims <- c(3.7, 4.7)
sapply(list(list("ds", c(ns.color, cnv.colors[3]), c("gray30", cnv.blacks[3]),
                 "min(pHaplo, pTriplo)", ds.coefs, ds.coefs.noConstr),
            list("hits", cnv.colors[1:2], cnv.blacks[1:2],
                 "pTriplo-pHaplo", hits.coefs, hits.coefs.noConstr)),
       function(vals){
         # Including all features
         pdf(paste(out.prefix, "gradient_regression.top_coefs", vals[[1]], "pdf", sep="."),
             height=pdf.dims[1], width=pdf.dims[2])
         plot.top.coefs(vals[[5]], feat.meta, colors=vals[[2]], borders=vals[[3]],
                        x.title=c(bquote(hat(beta)[.(vals[[4]])])), N=coefs.to.plot)
         dev.off()
         # Excluding constraint features
         pdf(paste(out.prefix, "gradient_regression.top_coefs", vals[[1]],
                   "no_constraint_features.pdf", sep="."),
             height=pdf.dims[1], width=pdf.dims[2])
         plot.top.coefs(vals[[6]], feat.meta, colors=vals[[2]], borders=vals[[3]],
                        x.title=c(bquote(hat(beta)[.(vals[[4]])])), N=coefs.to.plot)
         dev.off()
       })

# Also generate raw Spearman correlation coefficients for more interpretable feature prioritizations
# TODO: fix top x-axis labels
sapply(list(list("ds", c(ns.color, cnv.colors[3]), c("gray30", cnv.blacks[3]),
                 "min(pHaplo, pTriplo)", ds.cor, ds.cor.noConstr),
            list("hits", cnv.colors[1:2], cnv.blacks[1:2],
                 "pTriplo-pHaplo", hits.cor, hits.cor.noConstr)),
       function(vals){
         # Including all features
         pdf(paste(out.prefix, "gradient_regression.top_corr_coeffs", vals[[1]], "pdf", sep="."),
             height=pdf.dims[1], width=pdf.dims[2])
         plot.top.coefs(vals[[5]], feat.meta, colors=vals[[2]], borders=vals[[3]],
                        x.title=c(bquote("Spearman's" ~ rho ~ "vs." ~ .(vals[[4]]) )),
                        N=coefs.to.plot)
         dev.off()
         # Excluding constraint features
         pdf(paste(out.prefix, "gradient_regression.top_corr_coefs", vals[[1]],
                   "no_constraint_features.pdf", sep="."),
             height=pdf.dims[1], width=pdf.dims[2])
         plot.top.coefs(vals[[6]], feat.meta, colors=vals[[2]], borders=vals[[3]],
                        x.title=c(bquote("Spearman's" ~ rho ~ "vs." ~ .(vals[[4]]) )),
                        N=coefs.to.plot)
         dev.off()
       })


