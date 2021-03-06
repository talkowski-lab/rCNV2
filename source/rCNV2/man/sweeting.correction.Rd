% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assoc.R
\name{sweeting.correction}
\alias{sweeting.correction}
\title{Empirical continuity correction for meta-analyses}
\usage{
sweeting.correction(meta.df, cc.sum = 0.01)
}
\arguments{
\item{meta.df}{meta-analysis data frame as generated by \code{\link[=make.meta.df]{make.meta.df()}}}

\item{cc.sum}{scaling factor for total continuity correction to be distributed
among cases and controls}
}
\value{
data frame of meta-analysis stats after applying continuity correction
}
\description{
Apply empirical continuity correction to meta-analysis data frame per Sweeting
et al., Stat. Med., 2004 (section 3.3)
}
\seealso{
\code{\link[=make.meta.df]{make.meta.df()}}
}
