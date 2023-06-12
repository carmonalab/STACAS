#' Find integration anchors using STACAS
#'
#' This function computes anchors between datasets for single-cell data integration. It is based on the Seurat function
#' \code{FindIntegrationAnchors}, but is optimized for integration of heterogenous data sets containing only 
#' partially overlapping cells subsets. It also computes a measure of distance between candidate anchors (rPCA), 
#' which is combined with the Seurat's anchor weight by the factor \code{alpha}. Prior knowledge about
#' cell types can optionally be provided to guide anchor finding.
#' Give this information in the \code{cell.labels} metadata column. This annotation level, which can be incomplete
#' (set to NA for cells of unknown type), is used to penalize anchor pairs with inconsistent annotation.
#' The set of anchors returned by this function can then be passed to \code{IntegrateData.STACAS}
#' for dataset integration.
#'
#' @param object.list A list of Seurat objects. Anchors will be determined between pairs of objects, 
#' and can subsequently be used for Seurat dataset integration.
#' @param assay A vector containing the assay to use for each Seurat object in object.list.
#' If not specified, uses the default assay.
#' @param reference A vector specifying the object/s to be used as a reference
#' during integration. If NULL (default), all pairwise anchors are found (no
#' reference/s). If not NULL, the corresponding objects in \code{object.list}
#' will be used as references. When using a set of specified references, anchors
#' are first found between each query and each reference. The references are
#' then integrated through pairwise integration. Each query is then mapped to
#' the integrated reference.
#' @param max.seed.objects Number of objects to use as seeds to build the integration tree.
#' Automatically chooses the largest max.seed.objects datasets;
#' the remaining datasets will be added sequentially to the reference.
#' @param anchor.features Can be either: \itemize{
#'   \item{A numeric value. This will call \code{FindVariableFeatures.STACAS} to identify \code{anchor.features}
#'       that are consistently variable across datasets}
#'   \item{A pre-calculated vector of integration features to be used for anchor search.}}
#' @param genesBlockList  If \code{anchor.features} is numeric, \code{genesBlockList}
#'     optionally takes a (list of) vectors of gene names. These genes will be
#'     removed from the integration features. If set to "default",
#'     STACAS uses its internal list \code{data("genes.blocklist")}.
#'     This is useful to mitigate effect of genes associated with technical
#'     artifacts or batch effects (e.g. mitochondrial, heat-shock response).
#' @param dims The number of dimensions used for PCA reduction
#' @param k.anchor The number of neighbors to use for identifying anchors
#' @param k.score The number of neighbors to use for scoring anchors
#' @param alpha Weight on rPCA distance for rescoring (between 0 and 1).
#' @param anchor.coverage Center of logistic function, based on quantile value of rPCA distance distribution
#' @param correction.scale Scale factor for logistic function (multiplied by SD of rPCA distance distribution)
#' @param cell.labels A metadata column name, storing cell type annotations. These will be taken into account
#' for semi-supervised alignment (optional). Note that not all cells need to be annotated - please set
#' unannotated cells as NA or 'unknown' for this column. Cells with NA or 'unknown' cell labels will not be
#' penalized in semi-supervised alignment.
#' @param label.confidence How much you trust the provided cell labels (from 0 to 1).
#' @param seed Random seed for probabilistic anchor acceptance
#' @param verbose Print all output
#' 
#' @return Returns an AnchorSet object, which can be passed to \code{IntegrateData.STACAS}
#' @import Seurat
#' @importFrom pbapply pblapply
#' @export

FindAnchors.STACAS <- function (
  object.list = NULL,
  assay = NULL,
  reference = NULL,
  max.seed.objects = 10,
  anchor.features = 1000,
  genesBlockList = "default",
  dims = 30,
  k.anchor = 5,
  k.score = 30,
  alpha=0.8,
  anchor.coverage = 0.5,
  correction.scale = 2,  
  cell.labels = NULL,
  label.confidence = 1,
  seed = 123,
  verbose = TRUE
) {
  
  if (label.confidence<0 | label.confidence>1) {
    stop("label.confidence must be a number between 0 and 1")
  }
  if (anchor.coverage<0 | anchor.coverage>1) {
    stop("anchor.coverage must be a number between 0 and 1")
  }
  if (alpha<0 | alpha>1) {
    stop("alpha must be a number between 0 and 1")
  }
  
  if (length(dims)==1 & is.numeric(dims)) {
    dims.vec <- 1:dims
  } else if (is.numeric(dims)) {
    dims.vec <- dims
    dims <- max(dims)
  } else {
    stop("unsupported type for 'dims' parameter")
  }
  nobj <- length(object.list)
  
  #default assay, or user-defined assay
  if (!is.null(assay)) {
    if (length(assay) != nobj) {
      stop("If specifying the assay, please specify one assay per object in the object.list")
    }
    object.list <- sapply(
      X = 1:length(object.list),
      FUN = function(x) {
        DefaultAssay(object.list[[x]]) <- assay[x]
        return(object.list[[x]])
      }
    )
  } else {
    assay <- sapply(object.list, DefaultAssay)
  }

  #Calculate anchor genes
  if (is.numeric(anchor.features)) {
    n.this <- anchor.features
    if (verbose) {
      message("Computing ", anchor.features, " integration features")
    }
    check.genes(object.list)
    
    object.list <- lapply(object.list, function(x) {
      FindVariableFeatures.STACAS(x, nfeat = n.this, genesBlockList=genesBlockList)
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
  for (i in 1:nobj) {
    object.list[[i]] <- ScaleData(object.list[[i]], assay=assay[i], model.use="linear",
                                  do.center=FALSE, do.scale=FALSE,
                                  features = anchor.features, verbose=FALSE)
    if (verbose) {
      cat(paste0(" ",i,"/",nobj))
    }
    object.list[[i]] <- RunPCA(object.list[[i]], features = anchor.features, npcs=dims,
                               ndims.print = NA, nfeatures.print = NA, verbose=FALSE)
  }
  if (verbose) {
    cat("\nFinding integration anchors...\n")
  }
  
  #With a large list of input objects, select the top max.seed.objects to build a reference
  if (is.null(reference) & nobj > max.seed.objects) {
    sizes <- unlist(lapply(object.list, ncol))
    names(sizes) <- seq_along(sizes)
    sizes <- sort(sizes, decreasing = T)[1:max.seed.objects]
    reference <- sort(as.numeric(names(sizes)))
  }
  
  #Find pairwise anchors and keep distance information
  ref.anchors <- FindIntegrationAnchors.wdist(object.list, dims = dims.vec, k.anchor = k.anchor,
                                              anchor.features=anchor.features,
                                              reference = reference,
                                              assay=assay, k.score=k.score, verbose=verbose)
  
  #store Seurat knn consistency score
  ref.anchors@anchors['knn.score'] <- ref.anchors@anchors['score']
  #average reciprocal distances
  mat <- ref.anchors@anchors[,c("dist1.2","dist2.1")]
  ref.anchors@anchors['dist.mean'] <- apply(mat, 1, mean)
  
  ref.anchors <- reweight_anchors(ref.anchors, alpha=alpha,
                                  dist.pct=anchor.coverage,
                                  dist.scale.factor=correction.scale)
  if (!is.null(cell.labels)) {
     ref.anchors <- inconsistent_anchors(ref.anchors, cell.labels, seed=seed,
                                         label.confidence=label.confidence,
                                         quantile_ss=0)
  }
  
  return(ref.anchors)
}

#' Integration tree generation 
#'
#' Build an integration tree by clustering samples in a hierarchical manner. Cumulative scoring among anchor pairs will be used as pairwise similarity criteria of samples.
#' 
#' @param anchorset Scored anchorsobtained from \code{FindAnchors.STACAS} and \code{FilterAnchors.STACAS} function
#' @param hclust.method Clustering method to be used (single, complete, average, ward) 
#' @param usecol Column name to be used to compute sample similarity. Default "score"
#' @param obj.names Option vector of names for objects in anchorset
#' @param method Aggregation method to be used among anchors for sample similarity computation. Default: weight.sum
#' @param semisupervised Whether to use cell type label information (if available)
#' @param plot Logical indicating if dendrogram must be plotted
#' @return An integration tree to be passed to the integration function.
#' @import Seurat
#' @importFrom stats hclust
#' @export

SampleTree.STACAS <- function (
    anchorset,
    obj.names = NULL,
    hclust.method = c("single","complete","ward.D2","average"),
    usecol = c("score","dist.mean"),
    method = c("weight.sum","counts"),
    semisupervised = TRUE,
    plot = TRUE
) {
  
  hclust.method <- hclust.method[1]
  method <- method[1]
  usecol <- usecol[1]
  
  object.list <- slot(object = anchorset, name = "object.list")
  reference.objects <- slot(object = anchorset, name = "reference.objects")
  anchors <- slot(object = anchorset, name = "anchors")
  objects.ncell <- sapply(X = object.list, FUN = ncol)
  offsets <- slot(object = anchorset, name = "offsets")
  
  #Single reference
  if (length(reference.objects) == 1) {
    return(reference.objects)
  }
  
  if (semisupervised & "Retain_ss" %in% colnames(anchors)) {
    anchors <- anchors[anchors$Retain_ss==TRUE,]
  }

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
  
  names <- obj.names[reference.objects]
  
  if (!is.null(names)) { 
    rownames(similarity.matrix) <- names
    colnames(similarity.matrix) <- names
  }
  distance.matrix <- as.dist(m = 1 - similarity.matrix)
  
  if(plot & length(reference.objects)>2) {
    plot(hclust(d = distance.matrix,method = hclust.method))
  }
  sample.tree <- hclust(d = distance.matrix,method = hclust.method)$merge
  sample.tree <- AdjustSampleTree.Seurat(x = sample.tree, reference.objects = reference.objects)
  
  #Sum of anchors between sets
  nanch <- list()
  names(object.list) <- as.character(-(1:length(object.list)))
  for (i in reference.objects) {
    nanch[[as.character(-i)]] <- strength_function(subset(anchors,dataset1==i),
                                                   method = method,
                                                   usecol=usecol)
  }
  
  #Which is the most connected dataset?
  base <- which.max(nanch)
  if (!is.null(names)) { 
    base <- names[base]
  }    
  message(sprintf("Building integration tree with base dataset: %s", base))
  
  for (r in 1:nrow(sample.tree)) {
    pair <- sample.tree[r, ]
    
    w1 <- nanch[[as.character(pair[1])]]
    w2 <- nanch[[as.character(pair[2])]]

    if (w2 > w1) {
      pair <- rev(pair)
      sample.tree[r, ] <- pair
    }
    
 #   nanch[[as.character(r)]] <- w1 + w2  #cumulative (weighted) # of anchors
    nanch[[as.character(r)]] <- max(w1,w2)  #keep weight of more connected of the pair 
  }
  return(sample.tree)
}

#' FindVariableFeatures.STACAS
#'
#' Select highly variable genes (HVG) from an expression matrix. Genes from a blocklist
#' (e.g. cell cycling genes, mitochondrial genes) can be excluded from the list of
#' variable genes, as well as genes with very low or very high average expression
#'
#' @param obj A Seurat object containing an expression matrix
#' @param nfeat Number of top HVG to be returned
#' @param genesBlocklist Optionally takes a list of vectors of gene names. These genes will be removed from initial HVG set. If set to "default",
#'     STACAS uses its internal list \code{data("genes.blocklist")}.
#'     This is useful to mitigate effect of genes associated with technical artifacts or batch effects
#'     (e.g. mitochondrial, heat-shock response). 
#'     If set to `NULL` no genes will be excluded
#' @param min.exp Minimum average normalized expression for HVG. If lower, the gene will be excluded
#' @param max.exp Maximum average normalized expression for HVG. If higher, the gene will be excluded
#' @return Returns a list of highly variable genes
#' @import Seurat
#' @export FindVariableFeatures.STACAS
#' 
FindVariableFeatures.STACAS <- function(
    obj,
    nfeat=1500,
    genesBlockList="default",
    min.exp=0.01,
    max.exp=3)
{
  
  assay <- DefaultAssay(obj)
  #Calculate a fixed number of HVG, then filtered to nfeat at the end
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = 10000, verbose=F)
  
  varfeat <- VariableFeatures(obj)
  
  if (is.list(genesBlockList)) {
    genes.block <- unlist(genesBlockList) #user-provided list
  } else if (is.vector(genesBlockList)) {
      if (genesBlockList[1] == "default") {
        genes.block <- unlist(get.blocklist(obj))  #default list
      } else {
        genes.block <- genesBlockList #user-provided vector
      }
  } else {
    genes.block <- NULL # No excluded genes
  }
  
  varfeat <- setdiff(varfeat, genes.block)
  
  #Also remove genes that are very poorly or always expressed (=not really variable genes)
  means <- apply(obj@assays[[assay]]@data[varfeat,], 1, mean)
  removeGenes2 <- names(means[means<min.exp | means>max.exp])
  
  varfeat <- setdiff(varfeat, removeGenes2)
  n <- min(length(varfeat), nfeat)
  
  VariableFeatures(obj) <- varfeat[1:n]
  
  return(obj)
}  

#' PlotAnchors.STACAS
#'
#' Plot distribution of rPCA distances between pairs of datasets
#'
#' @param ref.anchors A set of anchors calculated using \code{FindAnchors.STACAS}, containing the pairwise distances between anchors.
#' @param obj.names Vector of object names, one for each dataset in ref.anchors
#' @param anchor.coverage Quantile of rPCA distance distribution
#' @return A plot of the distribution of rPCA distances
#' @import ggridges
#' @import ggplot2
#' @importFrom grDevices rainbow
#' @export
#' 

PlotAnchors.STACAS <- function(
    ref.anchors = NULL,
    obj.names = NULL,
    anchor.coverage = 0.5
) {
  anchortab <- ref.anchors@anchors
  
  levs <- levels(as.factor(anchortab$dataset1))
  if(is.null(obj.names)) {
    obj.names <- levs
  }
  
  if(length(obj.names) != length(levs)) {
    stop("If you provide dataset names, they must be as many as the levels in the anchor set")
  }
  
  if(anchor.coverage<0 | anchor.coverage>1) {
    stop("Variable anchor.coverage must be a real number between 0 and 1")
  }
  
  dist.thr = quantile(anchortab$dist.mean, anchor.coverage)
  
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

#' IntegrateData.STACAS
#'
#' Integrate a list of datasets using STACAS anchors. Based on the \code{IntegrateData} function from Seurat.
#' This function requires that you have calculated a set of integration anchors using \code{FindAnchors.STACAS}.
#' To perform semi-supervised integration, run \code{FindAnchors.STACAS} with cell type annotations labels.
#' Integration anchors with inconsistent cell type will be excluded from integration, providing an
#' integrated space that is partially guided by prior information.
#'
#' @param anchorset A set of anchors calculated using \code{FindAnchors.STACAS}
#' @param new.assay.name Assay to store the integrated data
#' @param features.to.integrate Which genes to include in the corrected integrated space (def. variable genes)
#' @param dims Number of dimensions for local anchor weighting
#' @param k.weight Number of neighbors for local anchor weighting. Set \code{k.weight="max"} to disable local weighting
#' @param sample.tree Specify the order of integration. See \code{SampleTree.STACAS} to calculate an integration tree.
#' @param hclust.method Clustering method for integration tree (single, complete, average, ward) 
#' @param semisupervised Whether to use cell type label information (if available)
#' @param verbose Print progress bar and output
#' @return Returns a \code{Seurat} object with a new integrated Assay, with batch-corrected expression values
#' @import Seurat
#' @importFrom pbapply pblapply
#' @export IntegrateData.STACAS

IntegrateData.STACAS <- function(
    anchorset,
    new.assay.name = "integrated",
    features.to.integrate = NULL,
    dims = 30,
    k.weight = 100,
    sample.tree = NULL,
    hclust.method = c("single","complete","ward.D2","average"),
    semisupervised = TRUE,
    verbose = TRUE
) {
  
  if (length(dims)==1 & is.numeric(dims)) {
    dims.vec <- 1:dims
  } else if (is.numeric(dims)) {
    dims.vec <- dims
    dims <- max(dims)
  } else {
    stop("unsupported type for 'dims' parameter")
  }

  # default integration tree
  if (is.null(sample.tree)) {
    sample.tree <- SampleTree.STACAS(
      anchorset = anchorset,
      semisupervised = semisupervised,
      hclust.method = hclust.method[1],
      plot = FALSE
    )
  }
  
  if (semisupervised & "Retain_ss" %in% colnames(anchorset@anchors)) {
    anchorset@anchors <- anchorset@anchors[anchorset@anchors$Retain_ss==TRUE,]
  }
  
  reference.datasets <- slot(object = anchorset, name = 'reference.objects')
  nobj <- length(anchorset@object.list)
  anchors <- slot(object = anchorset, name = 'anchors')
  features <- slot(object = anchorset, name = "anchor.features")
  
  unintegrated <- suppressWarnings(expr = merge(
    x = anchorset@object.list[[1]],
    y = anchorset@object.list[2:nobj]
  ))

  if (!is.null(x = features.to.integrate)) {
    features.to.integrate <- intersect(
      x = features.to.integrate,
      y = Reduce(
        f = intersect,
        x = lapply(
          X = anchorset@object.list,
          FUN = rownames
        )
      )
    )
  }
  
  # perform pairwise integration
  reference.integrated <- PairwiseIntegrateReference.STACAS(
    anchorset = anchorset,
    new.assay.name = new.assay.name,
    features = features,
    features.to.integrate = features.to.integrate,
    dims = dims.vec,
    k.weight = k.weight,
    sample.tree = sample.tree,
    preserve.order = TRUE,
    verbose = verbose
  )
  
  if (length(reference.datasets) == nobj) {
    DefaultAssay(reference.integrated) <- new.assay.name
    VariableFeatures(reference.integrated) <- features
    
    reference.integrated[["FindIntegrationAnchors"]] <- slot(anchorset, name = "command")
    reference.integrated <- suppressWarnings(LogSeuratCommand(reference.integrated))
    
    return(reference.integrated)
  } else {
    active.assay <- DefaultAssay(anchorset@object.list[reference.datasets][[1]])
    reference.integrated[[active.assay]] <- NULL
    reference.integrated[[active.assay]] <- CreateAssayObject(
      data = GetAssayData(
        object = reference.integrated[[new.assay.name]],
        slot = 'data'
      ),
      check.matrix = FALSE
    )
    DefaultAssay(reference.integrated) <- active.assay
    reference.integrated[[new.assay.name]] <- NULL
    VariableFeatures(reference.integrated) <- features
    
    anchorset@object.list <- lapply(anchorset@object.list, function(x) {
      if (!active.assay %in% Assays(x)) {
        suppressWarnings(x[[active.assay]] <- x[[DefaultAssay(x)]])
      }
      x
    })
    
    integrated.data <- MapQueryData.local(
      anchorset = anchorset,
      reference = reference.integrated,
      new.assay.name = new.assay.name,
      features = features,
      features.to.integrate = features.to.integrate,
      dims = dims.vec,
      k.weight = k.weight,
      verbose = verbose
    )
    rm(reference.integrated)
    
    # Construct final assay object
    unintegrated[[new.assay.name]] <- CreateAssayObject(
      data = integrated.data,
      check.matrix = FALSE
    )
    rm(integrated.data)
    
    unintegrated <- SetIntegrationData(
      object = unintegrated,
      integration.name = "Integration",
      slot = "anchors",
      new.data = anchors
    )
    unintegrated <- SetIntegrationData(
      object = unintegrated,
      integration.name = "Integration",
      slot = "sample.tree",
      new.data = sample.tree
    )
    DefaultAssay(unintegrated) <- new.assay.name
    VariableFeatures(unintegrated) <- features
    unintegrated[["FindIntegrationAnchors"]] <- slot(anchorset, name = "command")
    unintegrated <- suppressWarnings(LogSeuratCommand(unintegrated))
    return(unintegrated)
  }
}


#' Run the STACAS integration pipeline
#'
#' This function is a wrapper for running the several steps required to integrate single-cell
#' datasets using STACAS: 1) Finding integration anchors; 2) Calculating the sample tree for
#' the order of dataset integration; 3) Dataset batch effect correction and integration
#'
#' @param object.list A list of Seurat objects. Anchors will be determined between pairs of objects, 
#' and can subsequently be used for Seurat dataset integration.
#' @param assay A vector containing the assay to use for each Seurat object in object.list.
#' If not specified, uses the default assay.
#' @param new.assay.name Assay to store the integrated data
#' @param reference A vector specifying the object/s to be used as a reference
#' during integration. If NULL (default), all pairwise anchors are found (no
#' reference/s). If not NULL, the corresponding objects in \code{object.list}
#' will be used as references. When using a set of specified references, anchors
#' are first found between each query and each reference. The references are
#' then integrated through pairwise integration. Each query is then mapped to
#' the integrated reference.
#' @param max.seed.objects Number of objects to use as seeds to
#' build the integration tree. Automatically chooses the largest max.seed.objects datasets;
#' the remaining datasets will be added sequentially to the reference.
#' @param anchor.features Can be either: \itemize{
#'   \item{A numeric value. This will call \code{Seurat::SelectIntegrationFeatures} to identify \code{anchor.features}
#'       genes for anchor finding.}
#'   \item{A pre-calculated vector of integration features to be used for anchor search.}}
#' @param genesBlockList  If \code{anchor.features} is numeric, \code{genesBlockList} optionally takes a list of vectors of
#'     gene names. These genes will be removed from the integration features. If set to "default",
#'     STACAS uses its internal list \code{data("genes.blocklist")}.
#'     This is useful to mitigate effect of genes associated with technical artifacts or batch effects
#'     (e.g. mitochondrial, heat-shock response). 
#' @param dims The number of dimensions used for PCA reduction
#' @param k.anchor The number of neighbors to use for identifying anchors
#' @param k.score The number of neighbors to use for scoring anchors
#' @param k.weight Number of neighbors for local anchor weighting. Set \code{k.weight="max"} to disable local weighting
#' @param alpha Weight on rPCA distance for rescoring (between 0 and 1).
#' @param anchor.coverage Center of logistic function, based on quantile value of rPCA distance distribution
#' @param correction.scale Scale factor for logistic function (multiplied by SD of rPCA distance distribution)
#' @param cell.labels A metadata column name, storing cell type annotations. These will be taken into account
#' for semi-supervised alignment (optional). Cells annotated as NA or NULL will not be penalized in semi-supervised
#' alignment
#' @param label.confidence How much you trust the provided cell labels (from 0 to 1).
#' @param hclust.method Clustering method for integration tree (single, complete, average, ward) 
#' @param seed Random seed for probabilistic anchor acceptance
#' @param verbose Print all output
#' 
#' @return Returns a \code{Seurat} object with a new integrated Assay. Also, centered, scaled variable features data are returned in the scale.data slot, and the pca of these batch-corrected scale data in the pca `reduction` slot 
#' @import Seurat
#' @export

Run.STACAS <- function (
    object.list = NULL,
    assay = NULL,
    new.assay.name = "integrated",
    reference = NULL,
    max.seed.objects = 10,
    anchor.features = 1000,
    genesBlockList = "default",
    dims = 30,
    k.anchor = 5,
    k.score = 30,
    k.weight = 100,
    alpha=0.8,
    anchor.coverage = 0.5,
    correction.scale = 2,  
    cell.labels = NULL,
    label.confidence = 1,
    hclust.method = c("single","complete","ward.D2","average"),
    seed = 123,
    verbose = FALSE
) {
  
  if (length(dims)==1 & is.numeric(dims)) {
    dims.vec <- 1:dims
  } else if (is.numeric(dims)) {
    dims.vec <- dims
    dims <- max(dims)
  } else {
    stop("unsupported type for 'dims' parameter")
  }
  
  # 1. Find anchors
  stacas_anchors <- FindAnchors.STACAS(object.list, assay=assay,
                                    anchor.features=anchor.features,
                                    reference=reference,
                                    max.seed.objects = max.seed.objects,
                                    genesBlockList=genesBlockList, dims=dims,
                                    k.anchor=k.anchor, k.score=k.score,
                                    alpha=alpha, anchor.coverage=anchor.coverage,
                                    correction.scale=correction.scale,
                                    cell.labels=cell.labels, seed=seed,
                                    label.confidence=label.confidence, verbose=verbose)
  
  rm(object.list)
  
  if (is.null(cell.labels)) {
    semisupervised <- FALSE
  } else {
    semisupervised <- TRUE
  }
  
  # 2. Integration tree
  tree <- SampleTree.STACAS(
    anchorset = stacas_anchors,
    semisupervised = semisupervised,
    hclust.method = hclust.method[1],
    plot = FALSE
  )
  
  # 3. Integrate datasets
  integrated <- IntegrateData.STACAS(stacas_anchors, dims=dims, sample.tree=tree,
                                     new.assay.name = new.assay.name,
                                     k.weight = k.weight, semisupervised = semisupervised,
                                     features.to.integrate=stacas_anchors@anchor.features)
  
  rm(stacas_anchors)
  
  # 4. Calculate batch-corrected PCA space
  integrated <- ScaleData(integrated)
  integrated <- RunPCA(integrated, npcs=dims)
  
  return(integrated)
}

#' Standardize gene symbols
#'
#' Converts gene names of a Seurat single-cell object to a dictionary of
#' standard symbols. This function is useful prior to integration of datasets
#' from different studies, where gene names may be inconsistent.
#' 
#' @param obj A Seurat object
#' @param EnsemblGeneTable A data frame of gene name mappings. This should have
#'     the format of \href{https://www.ensembl.org/info/data/biomart/index.html}{Ensembl BioMart tables}
#'     with fields "Gene.name", "Gene.Synonym" and "Gene stable ID". See also
#'     the default conversion table in STACAS with \code{data(EnsemblGeneTable)}
#' @param EnsemblGeneFile If \code{EnsemblGeneTable==NULL}, read a gene mapping
#'     table from this file
#'     
#' @return Returns a Seurat object with standard gene names. Genes not found in
#'     the standard list are removed. Synonyms are accepted when
#'     the conversion is not ambiguous.
#' @import Seurat
#' @import R.utils
#' @importFrom data.table fread
#' @export

StandardizeGeneSymbols = function(obj, EnsemblGeneTable=NULL, EnsemblGeneFile=NULL){
  
  assay <- DefaultAssay(obj)
  #If file is given
  if (is.null(EnsemblGeneTable)) {
    if (is.null(EnsemblGeneFile)) {
      stop("Please provide EnsemblID table or file")
    }
    EnsemblGeneTable <- fread(EnsemblGeneFile)
  } 
  
  #Translate Ensembl IDs if necessary
  ens.format <- FALSE
  genes.in <- rownames(obj)
  ngenes <- length(genes.in)
  
  ens.count <- length(intersect(genes.in, EnsemblGeneTable[["Gene stable ID"]]))
  gname.count <- length(intersect(genes.in, EnsemblGeneTable[["Gene name"]]))
  
  max <- max(ens.count, gname.count)
  if (max < length(genes.in)/2) {
    warning("Over 50% of genes in input object not found in reference gene table")
  }
  if (ens.count > gname.count) {
    ens.format <- TRUE
  }
  
  if (ens.count > gname.count) {  #Input object has Ensembl IDs
    genes.tr <- EnsemblGeneTable[["Gene name"]][match(genes.in, EnsemblGeneTable[["Gene stable ID"]])]
    names(genes.tr) <- genes.in
    
    genes.tr <- genes.tr[!is.na(genes.tr) & genes.tr != ""]
  } else {
    genes.tr <- genes.in
    names(genes.tr) <- genes.in
  }
  
  ###### 1. First match dictionary 
  geneRef_dict <- EnsemblGeneTable[["Gene name"]]
  names(geneRef_dict) <- EnsemblGeneTable[["Gene Synonym"]]
  geneRef_dict <- geneRef_dict[!is.null(names(geneRef_dict))]
  
  message(paste("Number of genes in input object:", ngenes))
  genesAllowList1 <- genes.tr[!is.na(genes.tr) & genes.tr != "" &
                                genes.tr %in% EnsemblGeneTable[["Gene name"]]] #keep genes with standard Gene.name
  l <- length(genesAllowList1)
  
  message(sprintf("Number of genes with standard symbols: %i (%.2f%%)", l, l/ngenes*100))
  
  if (l < ngenes & !ens.format){
    message(paste("Examples of non-standard Gene.names:"))
    message(paste(head(genes.tr[ !genes.tr %in% EnsemblGeneTable[["Gene name"]] ])))
  }
  
  ###### 2. Search among synonyms
  genesAllowList2 <- genes.tr[!genes.tr %in% EnsemblGeneTable[["Gene name"]] & 
                                genes.tr %in% EnsemblGeneTable[["Gene Synonym"]]] # keep genes with accepted Gene.name synonym
  genesAllowList2.gn <- geneRef_dict[genesAllowList2] # translate Gene.Synonym to standard Gene.name
  
  message(paste("Additional number of genes with accepted Gene name synonym: ",length(genesAllowList2.gn)))
  
  #Names of genesAllowList contain IDs in matrix - elements contain the new names
  genesAllowList <- c(genesAllowList1,genesAllowList2.gn)
  
  ###### 3. Check for duplicates
  is.dup <- duplicated(genesAllowList)
  genesAllowList <- genesAllowList[!is.dup]
  message(sprintf("Number of duplicated Gene.name: %i (%.2f%%)", sum(is.dup), sum(is.dup)/ngenes*100))
  
  l <- length(genesAllowList)
  message(sprintf("Final number of genes: %i (%.2f%%)", l, l/ngenes*100))
  
  ###### 4. Subset matrix for allowed genes, and translate names
  rows.select <- rownames(obj@assays[[assay]]@counts)[rownames(obj@assays[[assay]]@counts) %in% names(genesAllowList)]
  obj <- obj[rows.select, ]
  rownames(obj@assays[[assay]]@data) <- unname(genesAllowList[rows.select])
  rownames(obj@assays[[assay]]@counts) <- unname(genesAllowList[rows.select])
  return(obj)
}

