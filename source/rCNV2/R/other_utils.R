#!/usr/bin/env Rscript

######################
#    rCNV Project    #
######################

# Copyright (c) 2021-Present Ryan L. Collins and the Talkowski Laboratory
# Distributed under terms of the MIT License (see LICENSE)
# Contact: Ryan L. Collins <rlcollins@g.harvard.edu>

# Miscellaneous project-wide utility functions


#' Calculate & format permuted P-value
#'
#' Computes a P-value by comparing an observed value to a vector of precomputed permutations
#'
#' @param perm.vals vector of results from precomputed permutations
#' @param obs.val empirically observed result
#' @param alternative specify sidedness of test (see `Details` below)
#'
#' @return list of two values:
#' 1. `$p` : numerical P-value
#' 2. `$formatted` : P-value formatted for plotting with `rCNV2::format.pval()`
#'
#' @export
calc.perm.p <- function(perm.vals, obs.val, alternative="greater"){
  if(alternative=="greater"){
    p <- length(which(perm.vals >= obs.val)) / length(perm.vals)
  }else{
    p <- length(which(perm.vals <= obs.val)) / length(perm.vals)
  }
  if(p==0){
    p.fmt <- format.pval(1/length(perm.vals), equality="<")
  }else{
    p.fmt <- format.pval(p)
  }
  return(list("p" = p, "formatted" = p.fmt))
}


#' Fit a robust linear regression with confidence interval
#'
#' Fits an outlier-robust version of univariate linear regression with confidence interval
#'
#' @param x independent variable vector
#' @param y dependent variable vector
#' @param conf width of confidence interval \[default: 0.95\]
#'
#' @details Effectively wraps `MASS:rlm()` and handles data reformatting as needed
#'
#' @return list of two objects:
#' 1. `$fit` : fitted model
#' 2. `$ci`: data frame with bounds of confidence interval
#'
#' @export
robust.lm <- function(x, y, conf=0.95){
  require(MASS, quietly=T)
  fit.df <- data.frame("Y"=y, "X"=x)
  xrange <- range(x[which(!is.infinite(x))], na.rm=T)
  xspan <- xrange[2] - xrange[1]
  fit <- MASS::rlm(Y ~ X, data=fit.df)
  pred.df <- data.frame("X"=seq(xrange[1] - 2*xspan, xrange[2] + 2*xspan, length.out=1000))
  pred <- predict(fit, pred.df, interval="confidence", level=conf)
  pred.out <- data.frame("x"=pred.df$X, "lower"=pred[, 2], "upper"=pred[, 3])
  return(list("fit" = fit, "ci" = pred.out))
}


#' Fit exponential decay
#'
#' Fits a univariate exponential decay with automatic parameter estimation
#'
#' @param x independent variable vector
#' @param y dependent variable vector
#'
#' @details Based on example from `https://rpubs.com/mengxu/exponential-model`
#'
#' @return Fitted model as `nls()` object with estimated parameters
#'
#' @export
fit.exp.decay <- function(x, y){
  # Following example on https://rpubs.com/mengxu/exponential-model

  train.df <- data.frame("x"=as.numeric(x),
                         "y"=as.numeric(y))
  train.df <- train.df[which(train.df$y > 0), ]

  # Estimate initial parameters
  theta.0 <- 0.5 * min(train.df$y)
  model.0 <- lm(log(y - theta.0) ~ x, data=train.df)
  alpha.0 <- exp(coef(model.0)[1])
  beta.0 <- coef(model.0)[2]
  start <- list(alpha = alpha.0, beta = beta.0, theta = theta.0)

  # Re-fit the model with estimated starting parameters
  return(nls(y ~ alpha * exp(beta * x) + theta, start = start, data=train.df,
             control=nls.control(maxiter=1000, warnOnly=T)))
}


#' Prune HPO Jaccard matrix
#'
#' Prune a matrix of Jaccard indexes for all HPO pairs based on maximum Jaccard index
#'
#' @param hpo.jac square matrix of Jaccard indexes for all HPOs
#' @param sample.sizes named vector of sample sizes for each HPO
#' @param max.jac maximum Jaccard index to retain \[default: 0.5\]
#' @param ignore.hpos vector of HPOs to ignore \[default: "HP:0000118" & "UNKNOWN"\]
#'
#' @return data.frame filtered on `max.jac`
#'
#' @details Pairs will be pruned while retaining the HPO with the larger sample size
#'
#' @seealso [load.hpo.jaccard.matrix]
#'
#' @export
prune.hpo.jaccard.matrix <- function(hpo.jac, sample.sizes, max.jac=0.5,
                                     ignore.hpos=c("HP:0000118", "UNKNOWN")){
  # Get lists of HPOs to process
  all.hpos <- union(colnames(hpo.jac), names(sample.sizes))
  elig.hpos <- all.hpos[which(!(all.hpos %in% ignore.hpos)
                              & all.hpos %in% colnames(hpo.jac))]

  # Melt matrix & order by decreasing Jaccard index
  jac <- as.data.frame(do.call("rbind", lapply(elig.hpos, function(ha){
    do.call("rbind", lapply(elig.hpos, function(hb){
      if(sample.sizes[ha] > sample.sizes[hb]){
        c(ha, hb, hpo.jac[ha, hb])
      }
    }))
  })))
  colnames(jac) <- c("hpo.A", "hpo.B", "jac")
  jac$jac <- as.numeric(jac$jac)
  jac <- jac[order(jac$jac, decreasing=T), ]

  # Prune HPOs sequentially by Jaccard index
  drop.hpos <- c()
  for(i in which(jac$jac > max.jac)){
    if(length(intersect(drop.hpos, as.vector(jac[i, c("hpo.A", "hpo.B")]))) == 0){
      drop.hpos <- union(drop.hpos, jac$hpo.B[i])
    }
  }

  # Return pruned Jaccard matrix
  if(length(drop.hpos) > 0){
    return(hpo.jac[-which(rownames(hpo.jac) %in% drop.hpos),
                   -which(colnames(hpo.jac) %in% drop.hpos)])
  }else{
    return(hpo.jac)
  }
}

