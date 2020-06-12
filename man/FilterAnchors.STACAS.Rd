\name{FilterAnchors.STACAS}
\alias{FilterAnchors.STACAS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Filter Anchors with STACAS}
\description{
    Filters integration anchors based on pairwise rPCA distance, calculated using \code{FindAnchors.STACAS}. The automatic heuristic threshold
    is normally rather conservative, you may want choose to apply a stricter threshold after inspecting the distance distributions made by
    \code{PlotAnchors.STACAS}
}
\usage{
FilterAnchors.STACAS(ref.anchors = NULL,
    dist.thr = NULL,
    dist.pct = 0.8)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ref.anchors}{
   A set of anchors calculated using \code{FindAnchors.STACAS}, containing the pairwise distances between anchors.
  }
   \item{dist.thr}{
    Distance threshold for anchor filtering. Distances are calculated by reciprocal PCA projection between pairs of datasets. If not specified, \code{dist.thr} defaults to the \code{dist.pct} percentile of the distance between the two closest datasets.
  }
  \item{dist.pct}{
    If \code{dist.thr} is not specified, calculate the filtering threshold as the \code{dist.pct} percentile of the distance between the two closest datasets (default=0.8)
  }
}
\details{
}
\value{
    Returns a filtered AnchorSet object, which can be directly applied to Seurat integration using \code{Seurat::IntegrateData}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
   \code{Run.STACAS}

   \code{FindAnchors.STACAS}

   \code{PlotAnchors.STACAS}
}
\examples{
    #default, automatic threshold
    ref.anchors.filtered <- FilterAnchors.STACAS(ref.anchors)
    #user-defined threshold for filtering
    ref.anchors.filtered <- FilterAnchors.STACAS(ref.anchors, dist.thr=0.7)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
