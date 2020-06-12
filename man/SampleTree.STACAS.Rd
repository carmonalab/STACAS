\name{SampleTree.STACAS}
\alias{SampleTree.STACAS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute integration tree with STACAS anchors
}
\description{
    Determine automatically a hierarchical tree for Seurat integration using anchors calculated in STACAS
}
\usage{
SampleTree.STACAS(anchorset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
   \item{anchorset}{
    The anchor object calculated using \code{Run.STACAS} or \code{FilterAnchors.STACAS}
  }
}
\details{
    Computes the integration tree based on dataset similarity and number of shared anchors between datasets/tree leaves.
}
\value{
    Returns a \code{sampletree} that can be passed as a parameter to the Seurat function \code{IntegrateData} (see example below)
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

\seealso{
   \code{Run.STACAS}

   \code{FindAnchors.STACAS}

   \code{PlotAnchors.STACAS}

   \code{FilterAnchors.STACAS}
}
\examples{
    data(STACAS.sampledata)
    STACAS.anchors <- Run.STACAS(STACAS.sampledata)
    tree <- SampleTree.STACAS(STACAS.anchors)

    #Then perform integration with Seurat
    ref.integrated <- Seurat::IntegrateData(anchorset=STACAS.anchors, sample.tree=tree, preserve.order=T)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
