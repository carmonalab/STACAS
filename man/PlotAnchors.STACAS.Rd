% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{PlotAnchors.STACAS}
\alias{PlotAnchors.STACAS}
\title{PlotAnchors.STACAS}
\usage{
PlotAnchors.STACAS(ref.anchors = NULL, obj.names = NULL, anchor.coverage = 0.5)
}
\arguments{
\item{ref.anchors}{A set of anchors calculated using \code{FindAnchors.STACAS}, containing the pairwise distances between anchors.}

\item{obj.names}{Vector of object names, one for each dataset in ref.anchors}

\item{anchor.coverage}{Quantile of rPCA distance distribution}
}
\value{
A plot of the distribution of rPCA distances
}
\description{
Plot distribution of rPCA distances between pairs of datasets
}
