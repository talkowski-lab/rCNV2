% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/highlight_utils.R
\name{pileup.cnvs.for.highlight}
\alias{pileup.cnvs.for.highlight}
\title{Format CNVs for locus highlight}
\usage{
pileup.cnvs.for.highlight(
  cnvs,
  start = NULL,
  end = NULL,
  dx = 100,
  cnv.height = 1,
  cnv.buffer = 0,
  bevel.switch.pct = 0.025,
  col = blueblack,
  highlight.hpo = NA,
  highlight.col = NULL
)
}
\arguments{
\item{cnvs}{data.frame of CNVs to be plotted}

\item{start}{left-most plotting coordinate}

\item{end}{right-most plotting coordinate}

\item{dx}{plotting resolution along X axis, specified as total number of steps [default: 100]}

\item{cnv.height}{relative height for each CNV [default: 1]}

\item{cnv.buffer}{relative buffer between adjacent CNVs [default: 0]}

\item{bevel.switch.pct}{fraction of start/end of each CNV to be beveled [default: 0.025]}

\item{col}{CNV color [default: blueblack]}

\item{highlight.hpo}{highlight CNVs from this HPO in a different color
[default: plot all CNVs in the same color]}

\item{highlight.color}{color to be used for \code{highlight.hpo}
[default: plot all CNVs in the same color]}
}
\value{
list ofwith the following two elements:
\verb{$cnvs}: plotting values for each individual CNV
\verb{$counts}: total number of CNVs overlapping each bin on the X axis
}
\description{
Transform CNV coordinates to plotting values (with colors) for locus highlights
}
\seealso{
\link{load.cnvs.from.region}
}
