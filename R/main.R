#' Find integration anchors using STACAS
#'
#' This function computes anchors between datasets for dataset integration. It is based on the Seurat function
#' \code{FindIntegrationAnchors}, but is optimized for integration of heterogenous data sets containing only 
#' partially overlapping cells subsets. It also returns a measure of distance between candidate anchors, 
#' which can be used at a later stage for anchor filtering using \code{FilterAnchors.STACAS}
#'
#' @param object.list A list of Seurat objects. Anchors will be determined between pairs of objects, 
#' and can subsequently be used for Seurat dataset integration.
#' @param assay A vector containing the assay to use for each Seurat object in object.list.
#' If not specified, uses the default assay.
#' @param reference  A vector specifying the (indices of the) objects to be used as a reference during integration.
#' If NULL (default), all pairwise anchors are found.
#' @param anchor.features Can be either: \itemize{
#'   \item{A numeric value. This will call \code{Seurat::SelectIntegrationFeatures} to select \code{anchor.feautures}
#'       genes for anchor finding.}
#'   \item{A pre-calculated vector of features to be used for anchor search.}}
#' @param dims The number of dimensions used for PCA reduction
#' @param normalization.method Which normalization method was used to prepare the data - either LogNormalize (default) or SCT
#' @param k.anchor The number of neighbors to use for identifying anchors
#' @param k.score The number of neighbors to use for scoring anchors
#' @param cell.labels A metadata column name, storing cell type annotations. These will be taken into account
#' for semi-supervised alignment (optional). Cells annotated as NA or NULL will not be penalized in semi-supervised
#' alignment
#' @param accept_rate_ss Probability of accepting inconsistent anchors with score above \code{quantile_ss} 
#' @param quantile_ss Distribution quantile on anchor scores to determine which inconsistent anchors to retain
#' @param verbose Print all output
#' 
#' @return Returns an AnchorSet object, which can be directly applied to Seurat integration using
#'  \code{Seurat::IntegrateData}, or optionally first filtered/weighted by anchor pairwise distance
#'  using \code{FilterAnchors.STACAS}
#' @import Seurat
#' @export

FindAnchors.STACAS <- function (
  object.list = NULL,
  assay = NULL,
  reference = NULL,
  anchor.features = 500,
  dims = 1:10,
  normalization.method = c("LogNormalize", "SCT"),
  k.anchor = 5,
  k.score = 30,
  cell.labels = NULL,
  accept_rate_ss = 0.5,
  quantile_ss = 0.8,
  verbose = FALSE
) {
  
  scale = FALSE
  reduction = "rpca"
  normalization.method <- match.arg(arg = normalization.method)
  
  #default assay, or user-defined assay
  if (!is.null(assay)) {
    if (length(assay) != length(object.list)) {
      stop("If specifying the assay, please specify one assay per object in the object.list")
    }
    object.list <- sapply(
      X = 1:length(object.list),
      FUN = function(x) {
        DefaultAssay(object = object.list[[x]]) <- assay[x]
        return(object.list[[x]])
      }
    )
  } else {
    assay <- sapply(X = object.list, FUN = DefaultAssay)
  }
  
  #Calculate anchor genes
  if (is.numeric(anchor.features)) {
    
    #Genes to exclude from variable features
    genes.block <- get.blocklist(object.list[[1]])
    
    n.this <- anchor.features
    if (verbose) {
      message("Computing ", anchor.features, " integration features")
    }
    object.list <- lapply(object.list, function(x) {
      select.variable.genes(x, nfeat = n.this, min.exp=0.01, max.exp=3, blacklist=genes.block)
    })
    
    #Combine variable features from multiple samples into single list
    anchor.features <- SelectIntegrationFeatures(
      object.list = object.list,
      nfeatures = n.this,
      assay = assay
    )
  }
  
  #prepare PCA without data-rescaling
  message("Preparing PCA embeddings for objects...")
  for (i in 1:length(object.list)) {
    object.list[[i]] <- ScaleData(object.list[[i]], assay=assay[i], model.use="linear", do.center=FALSE, do.scale=FALSE,
                                  features = row.names(object.list[[i]]), verbose=FALSE)
    cat(paste0(" ",i,"/",length(object.list)))
    object.list[[i]] <- RunPCA(object.list[[i]], features = anchor.features, ndims.print = NA, nfeatures.print = NA, verbose=FALSE)
  }
  cat("\n")
  
  #Find pairwise anchors and keep distance information
  if (is.null(reference)) {
    ref.anchors <- FindIntegrationAnchors.wdist(object.list, dims = dims, k.anchor = k.anchor, anchor.features=anchor.features,
                                                normalization.method = normalization.method,
                                                scale=scale, reduction=reduction, assay=assay, k.score=k.score, verbose=verbose)
  } else {
    ref.anchors <- FindIntegrationAnchors.wdist(object.list, reference=reference, dims = dims, k.anchor = k.anchor, anchor.features=anchor.features,
                                                normalization.method = normalization.method,
                                                scale=scale, reduction=reduction, assay=assay, k.score=k.score, verbose=verbose)
  }
  
  #average reciprocal distances
  mat <- ref.anchors@anchors[,c("dist1.2","dist2.1")]
  ref.anchors@anchors['dist.mean'] <- apply(mat, 1, mean)
  
  if (!is.null(cell.labels)) {
     ref.anchors <- inconsistent_anchors(ref.anchors, cell.labels,
                                         accept_rate_ss=accept_rate_ss, quantile_ss=quantile_ss)
  }
  
  return(ref.anchors)
}

#' FilterAnchors.STACAS
#'
#' Re-weighs and filters integration anchors based on pairwise rPCA distance and consistency in cell annotation.
#' These are calculated using \code{FindAnchors.STACAS}.
#'
#' @param ref.anchors A set of anchors calculated using \code{FindAnchors.STACAS}, containing the pairwise distances between anchors.
#' @param alpha Weight on rPCA distance for rescoring (between 0 and 1).
#' @param semi_supervised Use consistency between cell type labels to remove anchors.
#' @param dist.pct Center of logistic function, based on quantile value of rPCA distance distribution
#' @param dist.scale.factor Scale factor for logistic function (multiplied by SD of rPCA distance distribution)

#' @return A new anchor object with reweighted anchors scores, and optionally filtered by consistency of cell type annotation
#' @export
#' 

FilterAnchors.STACAS <- function(ref.anchors,
                                 alpha = 0.5,
                                 semi_supervised = FALSE,
                                 dist.pct = 0.5,
                                 dist.scale.factor = 2)
{
  df <- ref.anchors@anchors
  knn_score <- df$score
  epsilon <- 10^(-10)
  
  dist.mean.center = quantile(df$dist.mean,dist.pct)
  dist.mean.scale = sd(df$dist.mean,na.rm = T)/dist.scale.factor

  squash <- logistic(x = df$dist.mean, invert = T, 
                     center = dist.mean.center, 
                     scale = dist.mean.scale)
  df$score <-  alpha*squash + (1-alpha)*knn_score
  
  df$score[df$score < epsilon ] <- epsilon 
 
  if (semi_supervised) {
    if ("flag" %in% colnames(df)) {
       df <- df[df$flag==TRUE,]
    } else {
       warning("Cannot find 'flag' column in anchor object. Did you run FindAnchors.STACAS with cell.labels?
               Skipping semi-supervised alignment.")
    }
  }
  ref.anchors@anchors <- df
  return(ref.anchors)
}


#' Integration tree generation 
#'
#' Build an integration tree by clustering samples in a hierarchical manner. Cumulative scoring among anchor pairs will be used as pairwise similarity criteria of samples.
#' 
#' @param anchorset scored anchors  obtained from \code{scoring_anchors} STACAS function
#' @param hclust.method custering method to be used (complete, average, single, ward, etc) 
#' @param usecol colum name to be used to compute sample similarity. Default "score"
#' @param dist.hard.thr distance threshold used to mimic original STACAS behaviour
#' @param method agregation method to be used among anchors for sample similarity computation. Default: cum.sum
#' @param plot logical indicating if dendrogram must be ploted

#' @return A n integration tree to be passed to the integration function.
#' @export

SampleTree.STACAS <- function (
    anchorset,
    obj.names = NULL,
    hclust.method = c("average","single","complete"),
    usecol = "score",
    method = c("weight.sum","counts"),
    plot = TRUE
) {
  
  hclust.method <- hclust.method[1]
  method <- method[1]
  
  object.list <- slot(object = anchorset, name = "object.list")
  reference.objects <- slot(object = anchorset, name = "reference.objects")
  anchors <- slot(object = anchorset, name = "anchors")
  
  objects.ncell <- sapply(X = object.list, FUN = ncol)
  offsets <- slot(object = anchorset, name = "offsets")
  
  usecol = match.arg(arg = usecol ,choices = c("dist.mean","score"))
  method = match.arg(arg = method ,choices = c("weight.sum","counts"))
  
  ## Compute similarity matrix between datasets
  similarity.matrix <- weighted.Anchors.STACAS(
    anchor.df = anchors,
    offsets = offsets,
    obj.lengths = objects.ncell,
    usecol = usecol,
    method = method
  )
  
  similarity.matrix <- similarity.matrix[reference.objects, reference.objects]
  similarity.matrix <- similarity.matrix+10^-6 #avoid 1/0
  
  if (!is.null(obj.names)) { 
    rownames(similarity.matrix) <- obj.names
    colnames(similarity.matrix) <- obj.names
  }
  #minor change: use 1- simil instead of 1/simil
  #distance.matrix <- as.dist(m = 1 / similarity.matrix)
  distance.matrix <- as.dist(m = 1 - similarity.matrix)
  
  if(plot){
    plot(hclust(d = distance.matrix,method = hclust.method))
  }
  sample.tree <- hclust(d = distance.matrix,method = hclust.method)$merge
  sample.tree <- AdjustSampleTree.Seurat(x = sample.tree, reference.objects = reference.objects)
  
  #Sum of anchors between sets
  nanch <- list()
  names(x = object.list) <- as.character(-(1:length(x = object.list)))
  for (i in 1:length(object.list)) {
    nanch[[as.character(-i)]] <- strength_function(subset(anchors,dataset1==i),method = method)
  }
  
  for (r in 1:nrow(sample.tree)) {
    pair <- sample.tree[r, ]
    
    length1 <- nanch[[as.character(pair[1])]]
    length2 <- nanch[[as.character(pair[2])]]

    if (length2 > length1) {
      pair <- rev(pair)
      sample.tree[r, ] <- pair
    }
    
    nanch[[as.character(r)]] <- length1 + length2  #cumulative (weighted) # of anchors
  }
  return(sample.tree)
}

#' PlotAnchors.STACAS
#'
#' Plot distribution of rPCA distances between pairs of datasets
#'
#' @param ref.anchors A set of anchors calculated using \code{FindAnchors.STACAS}, containing the pairwise distances between anchors.
#' @param obj.names Vector of object names, one for each dataset in ref.anchors
#' @param dist.pct Quantile of rPCA distance distribution
#' @return A plot of the distribution of rPCA distances
#' @export
#' 

PlotAnchors.STACAS <- function(
    ref.anchors = NULL,
    obj.names = NULL,
    dist.pct = 0.5
) {
  anchortab <- ref.anchors@anchors
  
  levs <- levels(as.factor(anchortab$dataset1))
  if(is.null(obj.names)) {
    obj.names <- levs
  }
  
  if(length(obj.names) != length(levs)) {
    stop("If you provide dataset names, they must be as many as the levels in the anchor set")
  }
  
  if(dist.pct<0 | dist.pct>1) {
    stop("Variable dist.pct must be a real number between 0 and 1")
  }
  
  dist.thr = quantile(anchortab$dist.mean, dist.pct)
  
  ###Make distribution plots
  anchortab.toprint <- anchortab[]
  
  a.tmp <- anchortab.toprint
  for (r in 1:dim(anchortab.toprint)[1]) {
    anchortab.toprint[r,"dataset1"] <- obj.names[a.tmp[r,"dataset1"]]
    anchortab.toprint[r,"dataset2"] <- obj.names[a.tmp[r,"dataset2"]]
  }
  rm(a.tmp)
  
  my.colors=rainbow(length(levs), s=0.7)
  names(my.colors) <- obj.names
  pll=list()
  for (ds1 in 1:length(levs)) {
    data = subset(anchortab.toprint, subset=dataset1==obj.names[ds1])
    
    pll[[ds1]] <- ggplot(data=data, aes(x=dist.mean, y=dataset2, fill=dataset2)) +
      geom_density_ridges(alpha=0.4, scale=0.9) +
      scale_fill_manual(values=my.colors) +
      theme_ridges() +
      theme(legend.position = "none") +
      geom_vline(xintercept = dist.thr, linetype="dashed", size=0.75) +
      ggtitle(sprintf("%s - thr=%.3f", obj.names[ds1], dist.thr))
    
  }
  return(pll)
}


#Wrapper to be updated

Run.STACAS <- function (
    object.list = NULL,
    assay = NULL,
    reference = NULL,
    anchor.features = 500,
    dims = 1:10,
    normalization.method = c("LogNormalize", "SCT"),
    dist.pct = 0.8,
    k.anchor = 5,
    k.score = 30,
    plot.file = "anchor.dist.mean.pairwise.png",
    verbose = TRUE
) {
  normalization.method <- match.arg(arg = normalization.method)
  
  ref.anchors <- FindAnchors.STACAS(object.list, dims=dims, k.anchor=k.anchor, k.score=k.score, assay=assay,
                                    normalization.method = normalization.method,
                                    reference=reference, anchor.features=anchor.features, verbose=verbose)
  
  if (is.null(names(object.list))) {
    names <- 1:length(object.list)
  } else {
    names <- names(object.list)
  }
  
  plots <- PlotAnchors.STACAS(ref.anchors, obj.names=names,  dist.pct = dist.pct)
  
  g.cols <- 3
  g.rows <- as.integer((length(plots)+2)/g.cols)
  g <- do.call("arrangeGrob", c(plots, ncol=g.cols, nrow=g.rows))
  ggsave(plot.file, plot=g, width = 6*g.cols, height = 6*g.rows)
  
  plot.new()
  grid.draw(g)
  
  ref.anchors.filtered <- FilterAnchors.STACAS(ref.anchors,  dist.pct = dist.pct)
  
  return(ref.anchors.filtered)
}