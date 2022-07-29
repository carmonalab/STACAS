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
#'   \item{A numeric value. This will call \code{FindVariableFeatures.STACAS} to identify \code{anchor.features}
#'       that are consistently variable across datasets}
#'   \item{A pre-calculated vector of integration features to be used for anchor search.}}
#' @param genesBlockList  If \code{anchor.features} is numeric, \code{genesBlockList} optionally takes a list of vectors of
#'     gene names. These genes will be removed from the integration features. If set to "default",
#'     STACAS uses its internal list \code{data("genes.blocklist")}.
#'     This is useful to mitigate effect of genes associated with technical artifacts or batch effects
#'     (e.g. mitochondrial, heat-shock response). 
#' @param dims The number of dimensions used for PCA reduction
#' @param normalization.method Which normalization method was used to prepare the data - either LogNormalize (default) or SCT
#' @param k.anchor The number of neighbors to use for identifying anchors
#' @param k.score The number of neighbors to use for scoring anchors
#' @param alpha Weight on rPCA distance for rescoring (between 0 and 1).
#' @param anchor.coverage Center of logistic function, based on quantile value of rPCA distance distribution
#' @param correction.scale Scale factor for logistic function (multiplied by SD of rPCA distance distribution)
#' @param cell.labels A metadata column name, storing cell type annotations. These will be taken into account
#' for semi-supervised alignment (optional). Cells annotated as NA or NULL will not be penalized in semi-supervised
#' alignment
#' @param label.confidence How much you trust the provided cell labels (from 0 to 1).
#' @param future.maxSize For multi-core functionality, maximum allowed total size (in Gb) of global variables.
#'      To be incremented if required by \code{future.apply}
#' @param seed Random seed for probabilistic anchor acceptance
#' @param verbose Print all output
#' 
#' @return Returns an AnchorSet object, which can be passed to \code{IntegrateData.STACAS}
#' @import Seurat
#' @importFrom future nbrOfWorkers
#' @importFrom future.apply future_lapply
#' @importFrom pbapply pblapply
#' @export

FindAnchors.STACAS <- function (
  object.list = NULL,
  assay = NULL,
  reference = NULL,
  anchor.features = 1000,
  genesBlockList = "default",
  dims = 1:30,
  normalization.method = c("LogNormalize", "SCT"),
  k.anchor = 5,
  k.score = 30,
  alpha=0.8,
  anchor.coverage = 0.5,
  correction.scale = 2,  
  cell.labels = NULL,
  label.confidence = 1,
  future.maxSize = 16,
  seed = 123,
  verbose = FALSE
) {
  
  normalization.method <- match.arg(arg = normalization.method)
  
  if (label.confidence<0 | label.confidence>1) {
    stop("label.confidence must be a number between 0 and 1")
  }
  if (anchor.coverage<0 | anchor.coverage>1) {
    stop("anchor.coverage must be a number between 0 and 1")
  }
  if (alpha<0 | alpha>1) {
    stop("alpha must be a number between 0 and 1")
  }
  
  options(future.globals.maxSize= future.maxSize*1000*1024^2)
  
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

  genes.conserved <- check.genes(object.list)
  
  #Calculate anchor genes
  if (is.numeric(anchor.features)) {
    #Genes to exclude from variable features
    if (is.null(genesBlockList)) { 
      genes.block <- NULL #no excluded genes
    } else if (class(genesBlockList) == "list") {
      genes.block <- genesBlockList #user-provided list
    } else {
      genes.block <- get.blocklist(object.list[[1]])  #default list
    }

    n.this <- anchor.features
    if (verbose) {
      message("Computing ", anchor.features, " integration features")
    }
    object.list <- lapply(object.list, function(x) {
      FindVariableFeatures.STACAS(x, nfeat = n.this, genesBlockList=genes.block)
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
                                  features = anchor.features, verbose=FALSE)
    cat(paste0(" ",i,"/",length(object.list)))
    object.list[[i]] <- RunPCA(object.list[[i]], features = anchor.features, ndims.print = NA, nfeatures.print = NA, verbose=FALSE)
  }
  cat("\nFinding integration anchors...\n")
  
  #Find pairwise anchors and keep distance information
  if (is.null(reference)) {
    ref.anchors <- FindIntegrationAnchors.wdist(object.list, dims = dims, k.anchor = k.anchor, anchor.features=anchor.features,
                                                normalization.method = normalization.method,
                                                assay=assay, k.score=k.score, verbose=verbose)
  } else {
    ref.anchors <- FindIntegrationAnchors.wdist(object.list, reference=reference, dims = dims, k.anchor = k.anchor, anchor.features=anchor.features,
                                                normalization.method = normalization.method,
                                                assay=assay, k.score=k.score, verbose=verbose)
  }
  
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
#' @param hclust.method Clustering method to be used (complete, average, single, ward) 
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
    hclust.method = c("ward.D2","average","single","complete"),
    usecol = c("score","dist.mean"),
    method = c("weight.sum","counts"),
    semisupervised = TRUE,
    plot = TRUE
) {
  
  hclust.method <- hclust.method[1]
  method <- method[1]
  usecol = usecol[1]
  
  object.list <- slot(object = anchorset, name = "object.list")
  reference.objects <- slot(object = anchorset, name = "reference.objects")
  anchors <- slot(object = anchorset, name = "anchors")
  objects.ncell <- sapply(X = object.list, FUN = ncol)
  offsets <- slot(object = anchorset, name = "offsets")
  
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
  
  if (!is.null(obj.names)) { 
    rownames(similarity.matrix) <- obj.names
    colnames(similarity.matrix) <- obj.names
  }
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
    nanch[[as.character(-i)]] <- strength_function(subset(anchors,dataset1==i),
                                                   method = method,
                                                   usecol=usecol)
  }
  
  #Which is the most connected dataset?
  base <- which.max(nanch)
  if (!is.null(obj.names)) { 
    base <- obj.names[base]
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

#' Find variable features for STACAS
#'
#' Select highly variable genes (HVG) from an expression matrix. Genes from a blocklist
#' (e.g. cell cycling genes, mitochondrial genes) can be excluded from the list of
#' variable genes, as well as genes with very low or very high average expression.
#'
#' @param obj A Seurat object containing an expression matrix
#' @param nfeat Number of top HVG to be returned
#' @param genesBlocklist Optionally takes a list of vectors of gene names. These genes will be removed from initial HVG set. If set to "default",
#'     STACAS uses its internal list \code{data("genes.blocklist")}.
#'     This is useful to mitigate effect of genes associated with technical artifacts or batch effects
#'     (e.g. mitochondrial, heat-shock response). 
#'     If set to `NULL` no genes will be excluded
#' @param min.exp Minimum average normalized expression variable for HVG. If lower, the gene will be excluded
#' @param max.exp Maximum average normalized expression variable for HVG. If higher, the gene will be excluded
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
  
  #Calculate a fixed number of HVG, then filtered to nfeat at the end
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = 10000, verbose=F)
  
  varfeat <- obj@assays$RNA@var.features
  
  if (is.null(genesBlockList)) { 
    genes.block <- NULL #no excluded genes
  } else if (class(genesBlockList) == "list") {
    genes.block <- genesBlockList #user-provided list
  } else {
    genes.block <- get.blocklist(object.list[[1]])  #default list
  }
  
  varfeat <- setdiff(varfeat, unlist(genes.block))
  
  #Also remove genes that are very poorly or always expressed (=not really variable genes)
  means <- apply(obj@assays$RNA@data[varfeat,], 1, mean)
  removeGenes2 <- names(means[means<min.exp | means>max.exp])
  
  varfeat <- setdiff(varfeat, removeGenes2)
  n <- min(length(varfeat), nfeat)
  
  obj@assays$RNA@var.features <- varfeat[1:n]
  
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
#'
#' @param anchorset A set of anchors calculated using \code{FindAnchors.STACAS} and filtered using \code{FilterAnchors.STACAS}
#' @param new.assay.name Assay to store the integrated data
#' @param normalization.method How was the data normalized?
#' @param features.to.integrate Which genes to include in the corrected integrated space (def. variable genes)
#' @param dims Number of dimensions for local anchor weighting
#' @param k.weight Number of neighbors for local anchor weighting. Set \code{k.weight="max"} to disable local weighting
#' @param sample.tree Specify the order of integration. See \code{SampleTree.STACAS} to calculate an integration tree.
#' @param semisupervised Whether to use cell type label information (if available)
#' @param verbose Print progress bar and output
#' @return Returns a \code{Seurat} object with a new integrated Assay. If normalization.method = "LogNormalize", the integrated data is returned to the data slot and can be treated as log-normalized, corrected data. If normalization.method = "SCT", the integrated data is returned to the scale.data slot and can be treated as centered, corrected Pearson residuals.
#' @import Seurat
#' @export IntegrateData.STACAS

IntegrateData.STACAS <- function(
    anchorset,
    new.assay.name = "integrated",
    normalization.method = c("LogNormalize", "SCT"),
    features.to.integrate = NULL,
    dims = 1:30,
    k.weight = 100,
    sample.tree = NULL,
    semisupervised = TRUE,
    verbose = TRUE
) {

  normalization.method <- match.arg(arg = normalization.method)
  reference.datasets <- slot(object = anchorset, name = 'reference.objects')
  
  if (semisupervised & "Retain_ss" %in% colnames(anchorset@anchors)) {
    anchorset@anchors <- anchorset@anchors[anchorset@anchors$Retain_ss==TRUE,]
  }
  
  object.list <- slot(object = anchorset, name = 'object.list')
  anchors <- slot(object = anchorset, name = 'anchors')
  ref <- object.list[reference.datasets]
  features <- slot(object = anchorset, name = "anchor.features")
  unintegrated <- suppressWarnings(expr = merge(
    x = object.list[[1]],
    y = object.list[2:length(x = object.list)]
  ))
  if (!is.null(x = features.to.integrate)) {
    features.to.integrate <- intersect(
      x = features.to.integrate,
      y = Reduce(
        f = intersect,
        x = lapply(
          X = object.list,
          FUN = rownames
        )
      )
    )
  }
  if (normalization.method == "SCT") {
    model.list <- list()
    for (i in 1:length(x = object.list)) {
      assay <- DefaultAssay(object = object.list[[i]])
      if (length(x = setdiff(x = features.to.integrate, y = features)) != 0) {
        object.list[[i]] <- GetResidual(
          object = object.list[[i]],
          features = setdiff(x = features.to.integrate, y = features),
          verbose = verbose
        )
      }
      model.list[[i]] <- slot(object = object.list[[i]][[assay]], name = "SCTModel.list")
      object.list[[i]][[assay]] <- suppressWarnings(expr = CreateSCTAssayObject(
        data = GetAssayData(
          object = object.list[[i]],
          assay = assay,
          slot = "scale.data")
      )
      )
    }
    model.list <- unlist(x = model.list)
    slot(object = anchorset, name = "object.list") <- object.list
  }
  # perform pairwise integration of reference objects
  reference.integrated <- PairwiseIntegrateReference.STACAS(
    anchorset = anchorset,
    new.assay.name = new.assay.name,
    normalization.method = normalization.method,
    features = features,
    features.to.integrate = features.to.integrate,
    dims = dims,
    k.weight = k.weight,
    sample.tree = sample.tree,
    preserve.order = TRUE,
    verbose = verbose
  )
  
  # set SCT model
  if (normalization.method == "SCT") {
    if (is.null(x = Tool(object = reference.integrated, slot = "Integration"))) {
      reference.sample <- slot(object = anchorset, name = "reference.objects")
    } else {
      reference.sample <- SampleIntegrationOrder(
        tree = slot(
          object = reference.integrated,
          name = "tools"
        )$Integration@sample.tree
      )[1]
    }
    reference.cells <- Cells(x = object.list[[reference.sample]])
    reference.model <- NULL
    if (length(x = model.list) > 0) {
      reference.model <- sapply(X = model.list, FUN = function(model) {
        reference.check <- FALSE
        model.cells <- Cells(x = model)
        if (length(x = model.cells) > 0 &
            length(x = setdiff(x = model.cells, y = reference.cells)) == 0) {
          reference.check <- TRUE
        }
        return(reference.check)
      }
      )
      reference.model <- model.list[[which(reference.model)]]
    }
  }
  
  if (length(x = reference.datasets) == length(x = object.list)) {
    if (normalization.method == "SCT") {
      reference.integrated[[new.assay.name]] <- CreateSCTAssayObject(
        data = GetAssayData(object = reference.integrated, assay = new.assay.name, slot = "data"),
        scale.data = ScaleData(
          object = GetAssayData(object = reference.integrated, assay = new.assay.name, slot = "scale.data"),
          do.scale = FALSE,
          do.center = TRUE,
          verbose = FALSE),
        SCTModel.list = reference.model
      )
      levels(x =  reference.integrated[[new.assay.name]]) <- "refmodel"
      reference.integrated[[assay]] <- unintegrated[[assay]]
    }
    DefaultAssay(object = reference.integrated) <- new.assay.name
    VariableFeatures(object = reference.integrated) <- features
    reference.integrated[["FindIntegrationAnchors"]] <- slot(object = anchorset, name = "command")
    reference.integrated <- suppressWarnings(LogSeuratCommand(object = reference.integrated))
    return(reference.integrated)
  } else {
    active.assay <- DefaultAssay(object = ref[[1]])
    reference.integrated[[active.assay]] <- NULL
    # TODO: restore once check.matrix is in SeuratObject
    # reference.integrated[[active.assay]] <- CreateAssayObject(
    #   data = GetAssayData(
    #     object = reference.integrated[[new.assay.name]],
    #     slot = 'data'
    #   ),
    #   check.matrix = FALSE
    # )
    reference.integrated[[active.assay]] <- CreateAssayObject(
      data = GetAssayData(
        object = reference.integrated[[new.assay.name]],
        slot = 'data'
      )
    )
    DefaultAssay(object = reference.integrated) <- active.assay
    reference.integrated[[new.assay.name]] <- NULL
    VariableFeatures(object = reference.integrated) <- features
    # Extract the query objects (if any) and map to reference
    integrated.data <- MapQueryData(
      anchorset = anchorset,
      reference = reference.integrated,
      new.assay.name = new.assay.name,
      normalization.method = normalization.method,
      features = features,
      features.to.integrate = features.to.integrate,
      dims = dims,
      k.weight = k.weight,
      weight.reduction = weight.reduction,
      sd.weight = sd.weight,
      preserve.order = TRUE,
      eps = eps,
      verbose = verbose
    )
    
    # Construct final assay object
    # TODO: restore once check.matrix is in SeuratObject
    # integrated.assay <- CreateAssayObject(
    #   data = integrated.data,
    #   check.matrix = FALSE
    # )
    integrated.assay <- CreateAssayObject(
      data = integrated.data
    )
    if (normalization.method == "SCT") {
      integrated.assay <- CreateSCTAssayObject(
        data =  integrated.data,
        scale.data = ScaleData(
          object = integrated.data,
          do.scale = FALSE,
          do.center = TRUE,
          verbose = FALSE),
        SCTModel.list = reference.model
      )
      levels(x = integrated.assay) <- "refmodel"
    }
    unintegrated[[new.assay.name]] <- integrated.assay
    unintegrated <- SetIntegrationData(
      object = unintegrated,
      integration.name = "Integration",
      slot = "anchors",
      new.data = anchors
    )
    if (!is.null(x = Tool(object = reference.integrated, slot = "Integration"))) {
      sample.tree <- GetIntegrationData(
        object = reference.integrated,
        integration.name = "Integration",
        slot = "sample.tree"
      )
    }
    unintegrated <- SetIntegrationData(
      object = unintegrated,
      integration.name = "Integration",
      slot = "sample.tree",
      new.data = sample.tree
    )
    DefaultAssay(object = unintegrated) <- new.assay.name
    VariableFeatures(object = unintegrated) <- features
    unintegrated[["FindIntegrationAnchors"]] <- slot(object = anchorset, name = "command")
    unintegrated <- suppressWarnings(LogSeuratCommand(object = unintegrated))
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
#' @param reference  A vector specifying the (indices of the) objects to be used as a reference during integration.
#' If NULL (default), all pairwise anchors are found.
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
#' @param normalization.method Which normalization method was used to prepare the data - either LogNormalize (default) or SCT
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
#' @param seed Random seed for probabilistic anchor acceptance
#' @param verbose Print all output
#' 
#' @return Returns a \code{Seurat} object with a new integrated Assay. Also, centered, scaled variable features data are returned in the scale.data slot, and the pca of these batch-corrected scale data in the pca `reduction` slot 
#' @import Seurat
#' @export

Run.STACAS <- function (
    object.list = NULL,
    assay = NULL,
    reference = NULL,
    anchor.features = 1000,
    genesBlockList = "default",
    dims = 1:30,
    normalization.method = c("LogNormalize", "SCT"),
    k.anchor = 5,
    k.score = 30,
    k.weight = 100,
    alpha=0.8,
    anchor.coverage = 0.5,
    correction.scale = 2,  
    cell.labels = NULL,
    label.confidence = 1,
    seed = 123,
    verbose = FALSE
) {
  
  # 1. Find anchors
  stacas_anchors <- FindAnchors.STACAS(object.list, assay=assay, reference=reference,
                                    anchor.features=anchor.features,
                                    genesBlockList=genesBlockList, dims=dims,
                                    normalization.method=normalization.method,
                                    k.anchor=k.anchor, k.score=k.score,
                                    alpha=alpha, anchor.coverage=anchor.coverage,
                                    correction.scale=correction.scale,
                                    cell.labels=cell.labels, seed=seed,
                                    label.confidence=label.confidence, verbose=verbose)
  
  if (is.null(cell.labels)) {
    semisupervised <- FALSE
  } else {
    semisupervised <- TRUE
  }
  
  # 2. Integration tree
  tree <- SampleTree.STACAS(
    anchorset = stacas_anchors,
    semisupervised = semisupervised,
    plot = FALSE
  )
  
  # 3. Integrate datasets
  integrated <- IntegrateData.STACAS(stacas_anchors, dims=dims, sample.tree=tree,
                                     k.weight = k.weight, semisupervised = semisupervised,
                                     normalization.method=normalization.method,
                                     features.to.integrate=stacas_anchors@anchor.features)
  
  # 4. Calculate batch-corrected PCA space
  normalization.method <- match.arg(arg = normalization.method)
  if(normalization.method == "LogNormalize"){
    integrated <- ScaleData(integrated)
  }
  integrated <- RunPCA(integrated, npcs=max(dims))
  
  return(integrated)
}


