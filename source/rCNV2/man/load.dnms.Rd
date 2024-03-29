% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_io.R
\name{load.dnms}
\alias{load.dnms}
\title{Load dnm counts}
\usage{
load.dnms(dnms.in, mutrates)
}
\arguments{
\item{dnms.in}{.tsv of dnms per gene}

\item{mutrates}{data.frame of gene specific mutation rates as loaded by \code{\link[=load.mutrates]{load.mutrates()}}}
}
\value{
data.frame of raw and normalized dnms per gene
}
\description{
Load count of de novo mutations for a single study and residualize vs. expected
}
\seealso{
\code{\link[=load.mutrates]{load.mutrates()}}
}
