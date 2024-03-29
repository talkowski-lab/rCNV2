% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/assoc.R
\name{load.otf.meta.dat}
\alias{load.otf.meta.dat}
\title{Load on-the-fly meta-analysis data}
\usage{
load.otf.meta.dat(meta.inputs.in)
}
\arguments{
\item{meta.inputs.in}{path to input .tsv (see \code{details})}
}
\value{
list with one element per row in \code{meta.inputs.in}.
Each cohort has four named features:
\itemize{
\item \code{cnvs} : deletions and duplications as loaded by \code{\link[=load.otf.cohort.cnvs]{load.otf.cohort.cnvs()}}
\item \code{n_case} : vector of number of cases per cohort
\item \code{n_control} : vector of number of controls per cohort
\item \code{cohorts} : vector of cohort names
}
}
\description{
Load all data required for on-the-fly meta-analysis
}
\details{
\code{meta.inputs.in} is expected to be a .tsv with the following five columns:
\enumerate{
\item \code{cohort} : name of cohort
\item \code{n_case} : number of case samples
\item \code{n_control} : number of control samples
\item \code{DEL_path} : path to deletion BED file
\item \code{DUP_path} : path to duplication BED file
}
}
\seealso{
\code{\link[=load.otf.cohort.cnvs]{load.otf.cohort.cnvs()}}
}
