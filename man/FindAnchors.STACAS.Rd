\name{FindAnchors.STACAS}
\alias{FindAnchors.STACAS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find Anchors in STACAS}
\description{
    \code{FindAnchors.STACAS} computes anchors between datasets for dataset integration. It is based on the Seurat function \code{FindIntegrationAnchors},
    but is optimized for integration of heterogenous data sets containing only partially overlapping cells subsets. It also returns a measure of distance
    between candidate anchors, which can be used at a later stage for anchor filtering using \code{FilterAnchors.STACAS}

}
\usage{
FindAnchors.STACAS(object.list = NULL,
    assay = NULL,
    reference = NULL,
    anchor.features = 500,
    dims = 1:10,
    k.anchor = 5,
    k.score = 30,
    verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object.list}{
    A list of Seurat objects. Anchors will be determined between pairs of objects, and can subsequently be used for Seurat dataset integration.
  }
  \item{assay}{
    A vector containing the assay to use for each Seurat object in object.list. If not specified, uses the default assay.
  }
  \item{reference}{
    A vector specifying the (indices of the) objects to be used as a reference during integration. If NULL (default), all pairwise anchors are found.
  }
  \item{anchor.features}{
    Can be either:
    \itemize{
      \item{A numeric value. This will call \code{Seurat::SelectIntegrationFeatures} to select \code{anchor.feautures} genes for anchor finding}
      \item{A pre-calculated vector of feautures to be used for anchor search.}
    }
  }
  \item{dims}{
    The number of dimensions used for PCA reduction
  }
 \item{k.anchor}{
    The number of neighbors to use for identifying anchors
  }
  \item{k.score}{
    The number of neighbors to use for scoring anchors
  }
  \item{verbose}{
    Print detailed output and progress
  }
}
\details{
}
\value{
    Returns an AnchorSet object, which can be directly applied to Seurat integration using \code{Seurat::IntegrateData}, or optionally first filtered by
    anchor pairwise distanceusing \code{FilterAnchors.STACAS}
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

   \code{PlotAnchors.STACAS}

   \code{FilterAnchors.STACAS}
}
\examples{
     #Basic command
     ref.anchors <- FindAnchors.STACAS(object.list, anchor.features = 800)

     #Provide a pre-computed list of integration genes
     integr.genes <- Seurat::SelectIntegrationFeatures(object.list, nfeatures=800)
     ref.anchors <- FindAnchors.STACAS(object.list, anchor.features = integr.genes)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
