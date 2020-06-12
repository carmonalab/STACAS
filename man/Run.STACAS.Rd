\name{Run.STACAS}
\alias{Run.STACAS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Find integration anchors in STACAS}
\description{
    \code{Run.STACAS} is a wrapper for the identification, analysis and filtering of cell anchors for STACAS integration.
    For more control over paramenters and plots, you may choose to run the individual steps separately, as outlines in the Examples section below.
}
\usage{
Run.STACAS(object.list = NULL,
    assay = NULL, reference = NULL,
    anchor.features = 500,
    dims = 1:10,
    dist.thr = NULL,
    dist.pct = 0.8,
    k.anchor = 5,
    k.score = 30,
    plot.file = "anchor.dist.mean.pairwise.png",
    verbose = TRUE)
}

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
  \item{dist.thr}{
    Distance threshold for anchor filtering. Distances are calculated by reciprocal PCA projection between pairs of datasets. If not specified, \code{dist.thr} defaults to the \code{dist.pct} percentile of the distance between the two closest datasets.
  }
  \item{dist.pct}{
    If \code{dist.thr} is not specified, calculate the filtering threshold as the \code{dist.pct} percentile of the distance between the two closest datasets (default=0.8)
  }
  \item{k.anchor}{
    The number of neighbors to use for identifying anchors
  }
  \item{k.score}{
    The number of neighbors to use for scoring anchors
  }
  \item{plot.file}{
    Path to save the distribution plots
  }
  \item{verbose}{
    Print detailed output and progress
  }
}
\details{
}
\value{
    Returns an AnchorSet object, which can be directly applied to Seurat integration using \code{Seurat::IntegrateData}
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
   \code{FindAnchors.STACAS}

   \code{PlotAnchors.STACAS}

   \code{FilterAnchors.STACAS}
}
\examples{
### Run.STACAS is a wrapper for the functions FindAnchors.STACAS, PlotAnchors.STACAS and
#   FilterAnchors.STACAS as detailed below:

# 1) Determine anchors using FindAnchors.STACAS
    ref.anchors <- FindAnchors.STACAS(object.list, anchor.features = 500)

# 2) Plot distribution of anchor distances
    library(grid)
    library(gridExtra)
    plots <- PlotAnchors.STACAS(ref.anchors)
    g.cols <- 3
    g.rows <- as.integer((length(plots) + 2)/g.cols)
    g <- do.call("arrangeGrob", c(plots, ncol = g.cols, nrow = g.rows))
    ggsave("anchor.dist.mean.pairwise.png", plot = g, width = 6 *
        g.cols, height = 6 * g.rows)
    plot.new()
    grid.draw(g)

# 3) Filter anchors based on dist.thr - determined automatically or specified after
#     inspecting the distance plots from point 2

    ref.anchors.filtered <- FilterAnchors.STACAS(ref.anchors)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
