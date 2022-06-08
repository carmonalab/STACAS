# Logistic transformation
logistic <- function(x,scale = NULL, center = NULL, invert = F){
  a <- ifelse(invert,1, -1)
  if(is.null(scale))  scale = sd(x)
  if(is.null(center)) center = mean(x)
  sigm <- 1/(1+ exp(a*(x-center)/scale))
  return(sigm)    
}

# Mark anchor as inconsistent based on cell-type labeling 
# Anchors will be consider as inconsistent if they are linking cells classified as different cell populations. 
inconsistent_anchors <- function(anchors,
                                 cell.labels=NULL,
                                 accept_rate_ss=0.5,
                                 quantile_ss=0.5){
  
  if (! cell.labels %in% colnames(anchors@object.list[[1]]@meta.data)) {
    stop(sprintf("Please specify a valid metadata column with label annotations (cell.labels). %s not found", cell.labels))
  }
  
  labels <- data.frame()
  for(i in anchors@reference.objects){
    x <-anchors@object.list[[i]]
    labels <- rbind(labels, data.frame(dataset = i, ncell = 1:ncol(x), label = x@meta.data[,cell.labels]))
  }
  
  labels[,"dataset_cell"] <- paste0(labels[,"dataset"], "_", labels[,"ncell"])
  
  nas <- is.na(labels[,"label"]) | (labels[, "label"]== "NaN") | is.null(labels[, "label"])
  labels[nas, "label"] <- "unknown"
  
  df <- anchors@anchors
  df[,"dataset1_cell1"] <- paste0(df[,"dataset1"], "_", df[,"cell1"])
  df[,"dataset2_cell2"] <- paste0(df[,"dataset2"], "_", df[,"cell2"])
  
  order1 <- match(df$dataset1_cell1, labels$dataset_cell)
  df[,"Label_1"] <- labels[order1, "label"]
  
  order2 <- match(df$dataset2_cell2, labels$dataset_cell)
  df[,"Label_2"] <- labels[order2, "label"]
  
  #Flag consistency of labels
  match_ok <- df$Label_1 == df$Label_2 | df$Label_1 == "unknown" | df$Label_2 == "unknown"
  
  df$flag <- match_ok
  
  #Apply Boltzmann-based rejection based on deterministic label agreement
  if (quantile_ss < 1) {
    df <- boltzmann_based_rejection(anchors=df, accept_rate = accept_rate_ss, q = quantile_ss)
  }
  anchors@anchors <- df
  return(anchors)
}

# Boltzmann-based rejection of inconsistent anchors
# q: quantile score of non-inconsistent anchors
# accept_rate: probability of accepting anchors above the quantile score

boltzmann_based_rejection <- function(anchors, accept_rate = 0.5, q = 0.5){
  
  e0 <- quantile(anchors$score, q)                  #scale is given by 'q' quantile score
  DE <- anchors[anchors$flag==FALSE, "score"] - e0   #DE is given by the actual score of inconsistent anchors.   
  kT <- quantile(DE[DE>0],1-accept_rate)
  rejection_prob <- vapply(exp(-DE/kT), function(x){min(1,x)}, numeric(1))
  
  accept <- rejection_prob**(accept_rate) < runif(length(DE), 0, 1)
  
  #Add this column
  flag_bz <- anchors$flag
  flag_bz[anchors$flag==FALSE] <- accept
  anchors['flag'] <- flag_bz
  
  return(anchors)
}

# To be used in the distance matrix among datasets (when building SampleTree integration)
strength_function <- function(anchor.df,
                              method = "weight.sum", 
                              usecol = "dist.mean"
                              ){
  method = match.arg(arg = method ,choices = c("counts","mean.score","weight.sum"))
  usecol = match.arg(arg = usecol ,choices = c("dist.mean","score"))
  
  if(method== "counts") {
    strength <- nrow(anchor.df)
  }
  
  if(method== "mean.score") {
    if(usecol == "dist.mean")    strength <- 1-anchor.df[,usecol]%>%mean()
    if(usecol == "score")    strength <- anchor.df[,usecol]%>%mean()
  }

  if(method== "weight.sum") {
    if(usecol == "dist.mean")  strength <- (1-anchor.df[,usecol])%>%sum()
    if(usecol == "score")  strength <- (anchor.df[,usecol])%>%sum()
  }
  return(strength)    
}

#aux function to determine anchor threshold from pairwise anchor score percentile
threshold_from_percentile <- function(anchors, dist.pct=0.8) {
  i=1
  min.pcntile <- vector()
  for (l1 in lev1) {
    j=1
    pcntile <- vector()
    for (l2 in lev2) {
      v <- subset(ref.anchors@anchors, subset=(dataset1==l1 & dataset2==l2))$dist.mean
      if(length(v)>0) {
        pcntile[j] <- sort(v)[dist.pct*length(v)]
        j=j+1
      }
    }
    min.pcntile[i] <- min(pcntile)
    i=i+1
  }
  dist.thr <- min(min.pcntile)
  return(dist.thr)
}

# Compute similarity matrix betwen datasets based on different criteria
weighted.Anchors.STACAS <- function(
  anchor.df,
  offsets,
  obj.lengths,
  usecol = "score",
  method = "weight.sum"
) {
  usecol = match.arg(arg = usecol ,choices = c("dist.mean","score"))
  method = match.arg(arg = method ,choices = c("weight.sum","counts"))
  
  similarity.matrix <- matrix(data = 0, ncol = length(x = offsets), nrow = length(x = offsets))
  similarity.matrix[upper.tri(x = similarity.matrix, diag = TRUE)] <- NA
  total.cells <- sum(obj.lengths)
  offsets <- c(offsets, total.cells)
  for (i in 1:nrow(x = similarity.matrix)){
    for (j in 1:ncol(x = similarity.matrix)){
      if (!is.na(x = similarity.matrix[i, j])){
        relevant.rows <- anchor.df[(anchor.df$dataset1 %in% c(i, j)) & (anchor.df$dataset2 %in% c(i, j)), ]
        score <- strength_function(relevant.rows, method = method, usecol = usecol)/2
        ncell <- min(obj.lengths[[i]], obj.lengths[[j]])
        similarity.matrix[i, j] <- score / ncell
      }
    }
  }
  return(similarity.matrix)
}

select.variable.genes = function(obj, nfeat=1500, blacklist=NULL, min.exp=0.01, max.exp=3){
  
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = nfeat)
  
  varfeat <- obj@assays$RNA@var.features
  
  if (!is.null(blacklist)) {
    removeGenes1 <- varfeat[varfeat %in% unlist(blacklist)]
    varfeat <- setdiff(varfeat, removeGenes1)
  }
  #Also remove genes that are very poorly or always expressed (=not really variable genes)
  means <- apply(obj@assays$RNA@data[varfeat,], 1, mean)
  removeGenes2 <- names(means[means<min.exp | means>max.exp])
 
  varfeat <- setdiff(varfeat, removeGenes2 )
  obj@assays$RNA@var.features <- varfeat
  
  return(obj)
}  

get.blocklist = function(obj) {
  data("genes.blocklist")
  unlist(genes.blocklist$Mm)
  
  mm <- intersect(unlist(genes.blocklist$Mm), rownames(obj))
  hs <- intersect(unlist(genes.blocklist$Hs), rownames(obj))
  
  if (length(mm) > length(hs)) {
    return(genes.blocklist$Mm)
  } else {
    return(genes.blocklist$Hs)
  }
}

# Count anchors between data sets
# IMPORTED from Seurat 3.2.1
CountAnchors.Seurat <- function(
    anchor.df,
    offsets,
    obj.lengths
) {
  similarity.matrix <- matrix(data = 0, ncol = length(x = offsets), nrow = length(x = offsets))
  similarity.matrix[upper.tri(x = similarity.matrix, diag = TRUE)] <- NA
  total.cells <- sum(obj.lengths)
  offsets <- c(offsets, total.cells)
  for (i in 1:nrow(x = similarity.matrix)){
    for (j in 1:ncol(x = similarity.matrix)){
      if (!is.na(x = similarity.matrix[i, j])){
        relevant.rows <- anchor.df[(anchor.df$dataset1 %in% c(i, j)) & (anchor.df$dataset2 %in% c(i, j)), ]
        score <- nrow(x = relevant.rows)
        ncell <- min(obj.lengths[[i]], obj.lengths[[j]])
        similarity.matrix[i, j] <- score / ncell
      }
    }
  }
  return(similarity.matrix)
}

# Adjust sample tree to only include given reference objects
# IMPORTED from Seurat 3.2.1
AdjustSampleTree.Seurat <- function(x, reference.objects) {
  for (i in 1:nrow(x = x)) {
    obj.id <- -(x[i, ])
    if (obj.id[[1]] > 0) {
      x[i, 1] <- -(reference.objects[[obj.id[[1]]]])
    }
    if (obj.id[[2]] > 0) {
      x[i, 2] <- -(reference.objects[[obj.id[[2]]]])
    }
  }
  return(x)
}

###Modified function to return anchor distances together with anchor pairs - distances can be used to filter anchors at a later stage
FindIntegrationAnchors.wdist <- function(
    object.list = NULL,
    assay = NULL,
    reference = NULL,
    anchor.features = 2000,
    scale = TRUE,
    normalization.method = c("LogNormalize", "SCT"),
    sct.clip.range = NULL,
    reduction = c("cca", "rpca"),
    l2.norm = TRUE,
    dims = 1:10,
    k.anchor = 5,
    k.filter = 200,
    k.score = 30,
    max.features = 200,
    nn.method = "rann",
    eps = 0,
    verbose = TRUE
) {
  normalization.method <- match.arg(arg = normalization.method)
  reduction <- match.arg(arg = reduction)
  if (reduction == "rpca") {
    reduction <- "pca"
  }
  my.lapply <- ifelse(
    test = verbose && future::nbrOfWorkers() == 1,
    yes = pbapply::pblapply,
    no = future.apply::future_lapply
  )
  object.ncells <- sapply(X = object.list, FUN = function(x) dim(x = x)[2])
  if (any(object.ncells <= max(dims))) {
    bad.obs <- which(x = object.ncells <= max(dims))
    stop("Max dimension too large: objects ", paste(bad.obs, collapse = ", "),
         " contain fewer than ", max(dims), " cells. \n Please specify a",
         " maximum dimensions that is less than the number of cells in any ",
         "object (", min(object.ncells), ").")
  }
  if (!is.null(x = assay)) {
    if (length(x = assay) != length(x = object.list)) {
      stop("If specifying the assay, please specify one assay per object in the object.list")
    }
    object.list <- sapply(
      X = 1:length(x = object.list),
      FUN = function(x) {
        DefaultAssay(object = object.list[[x]]) <- assay[x]
        return(object.list[[x]])
      }
    )
  } else {
    assay <- sapply(X = object.list, FUN = DefaultAssay)
  }
  object.list <- Seurat:::CheckDuplicateCellNames(object.list = object.list)
  
  slot <- "data"
  if (normalization.method == "SCT") {
    slot <- "scale.data"
    scale <- FALSE
    if (is.numeric(x = anchor.features)) {
      stop("Please specify the anchor.features to be used. The expected ",
           "workflow for integratinge assays produced by SCTransform is ",
           "SelectIntegrationFeatures -> PrepSCTIntegration -> ",
           "FindIntegrationAnchors.")
    }
    sct.check <- sapply(
      X = 1:length(x = object.list),
      FUN = function(x) {
        sct.cmd <- grep(
          pattern = 'PrepSCTIntegration',
          x = Command(object = object.list[[x]]),
          value = TRUE
        )
        # check assay has gone through PrepSCTIntegration
        if (!any(grepl(pattern = "PrepSCTIntegration", x = Command(object = object.list[[x]]))) ||
            Command(object = object.list[[x]], command = sct.cmd, value = "assay") != assay[x]) {
          stop("Object ", x, " assay - ", assay[x], " has not been processed ",
               "by PrepSCTIntegration. Please run PrepSCTIntegration prior to ",
               "FindIntegrationAnchors if using assays generated by SCTransform.", call. = FALSE)
        }
        # check that the correct features are being used
        if (all(Command(object = object.list[[x]], command = sct.cmd, value = "anchor.features") != anchor.features)) {
          stop("Object ", x, " assay - ", assay[x], " was processed using a ",
               "different feature set than in PrepSCTIntegration. Please rerun ",
               "PrepSCTIntegration with the same anchor.features for all objects in ",
               "the object.list.", call. = FALSE)
        }
      }
    )
  }
  if (is.numeric(x = anchor.features) && normalization.method != "SCT") {
    if (verbose) {
      message("Computing ", anchor.features, " integration features")
    }
    
    #NOTE: improve on the following function to exclude blacklisted genes (mito, ribo, cycling), as well as
    #genes that are expressed in very few or all cells
    anchor.features <- SelectIntegrationFeatures(
      object.list = object.list,
      nfeatures = anchor.features,
      assay = assay
    )
  }
  if (scale) {
    if (verbose) {
      message("Scaling features for provided objects")
    }
    object.list <- my.lapply(
      X = object.list,
      FUN = function(object) {
        ScaleData(object = object, features = anchor.features, verbose = FALSE)
      }
    )
  }
  nn.reduction <- reduction
  # if using pca, only need to compute the internal neighborhood structure once
  # for each dataset
  internal.neighbors <- list()
  if (nn.reduction == "pca") {
    k.filter <- NA
    if (verbose) {
      message("Computing within dataset neighborhoods")
    }
    k.neighbor <- max(k.anchor, k.score)
    internal.neighbors <- my.lapply(
      X = 1:length(x = object.list),
      FUN = function(x) {
        Seurat:::NNHelper(
          data = Embeddings(object = object.list[[x]][[nn.reduction]])[, dims],
          k = k.neighbor + 1,
          method = nn.method,
          eps = eps
        )
      }
    )
  }
  
  # determine pairwise combinations
  combinations <- expand.grid(1:length(x = object.list), 1:length(x = object.list))
  combinations <- combinations[combinations$Var1 < combinations$Var2, , drop = FALSE]
  # determine the proper offsets for indexing anchors
  objects.ncell <- sapply(X = object.list, FUN = ncol)
  offsets <- as.vector(x = cumsum(x = c(0, objects.ncell)))[1:length(x = object.list)]
  if (is.null(x = reference)) {
    # case for all pairwise, leave the combinations matrix the same
    if (verbose) {
      message("Finding all pairwise anchors")
    }
  } else {
    reference <- unique(x = sort(x = reference))
    if (max(reference) > length(x = object.list)) {
      stop('Error: requested reference object ', max(reference), " but only ",
           length(x = object.list), " objects provided")
    }
    # modify the combinations matrix to retain only R-R and R-Q comparisons
    if (verbose) {
      message("Finding anchors between all query and reference datasets")
      ok.rows <- (combinations$Var1 %in% reference) | (combinations$Var2 %in% reference)
      combinations <- combinations[ok.rows, ]
    }
  }
  # determine all anchors
  plot.list=list()
  all.anchors <- my.lapply(
    X = 1:nrow(x = combinations),
    FUN = function(row) {
      i <- combinations[row, 1]
      j <- combinations[row, 2]
      object.1 <- DietSeurat(
        object = object.list[[i]],
        assays = assay[i],
        features = anchor.features,
        counts = FALSE,
        scale.data = TRUE,
        dimreducs = reduction
      )
      object.2 <- DietSeurat(
        object = object.list[[j]],
        assays = assay[j],
        features = anchor.features,
        counts = FALSE,
        scale.data = TRUE,
        dimreducs = reduction
      )
      # suppress key duplication warning
      suppressWarnings(object.1[["ToIntegrate"]] <- object.1[[assay[i]]])
      DefaultAssay(object = object.1) <- "ToIntegrate"
      if (reduction %in% Reductions(object = object.1)) {
        slot(object = object.1[[reduction]], name = "assay.used") <- "ToIntegrate"
      }
      object.1 <- DietSeurat(object = object.1, assays = "ToIntegrate", scale.data = TRUE, dimreducs = reduction)
      suppressWarnings(object.2[["ToIntegrate"]] <- object.2[[assay[j]]])
      DefaultAssay(object = object.2) <- "ToIntegrate"
      if (reduction %in% Reductions(object = object.2)) {
        slot(object = object.2[[reduction]], name = "assay.used") <- "ToIntegrate"
      }
      object.2 <- DietSeurat(object = object.2, assays = "ToIntegrate", scale.data = TRUE, dimreducs = reduction)
      object.pair <- switch(
        EXPR = reduction,
        'cca' = {
          object.pair <- RunCCA(
            object1 = object.1,
            object2 = object.2,
            assay1 = "ToIntegrate",
            assay2 = "ToIntegrate",
            features = anchor.features,
            num.cc = max(dims),
            renormalize = FALSE,
            rescale = FALSE,
            verbose = verbose
          )
          if (l2.norm){
            object.pair <- L2Dim(object = object.pair, reduction = reduction)
            reduction <- paste0(reduction, ".l2")
            nn.reduction <- reduction
          }
          reduction.2 <- character()
          object.pair
        },
        'pca' = {
          common.features <- intersect(
            x = rownames(x = Loadings(object = object.1[["pca"]])),
            y = rownames(x = Loadings(object = object.2[["pca"]]))
          )
          object.pair <- merge(x = object.1, y = object.2, merge.data = TRUE)
          
          projected.embeddings.1<- t(x = GetAssayData(object = object.1, slot = "scale.data")[common.features, ]) %*%
            Loadings(object = object.2[["pca"]])[common.features, ]
          
          object.pair[['projectedpca.1']] <- CreateDimReducObject(
            embeddings = rbind(projected.embeddings.1, Embeddings(object = object.2[["pca"]])),
            assay = DefaultAssay(object = object.1),
            key = "projectedpca1_"
          )
          projected.embeddings.2 <- t(x = GetAssayData(object = object.2, slot = "scale.data")[common.features, ]) %*%
            Loadings(object = object.1[["pca"]])[common.features, ]
          object.pair[['projectedpca.2']] <- CreateDimReducObject(
            embeddings = rbind(projected.embeddings.2, Embeddings(object = object.1[["pca"]])),
            assay = DefaultAssay(object = object.2),
            key = "projectedpca2_"
          )
          object.pair[["pca"]] <- CreateDimReducObject(
            embeddings = rbind(
              Embeddings(object = object.1[["pca"]]),
              Embeddings(object = object.2[["pca"]])),
            assay = DefaultAssay(object = object.1),
            key = "pca_"
          )
          
          reduction <- "projectedpca.1"
          reduction.2 <- "projectedpca.2"
          if (l2.norm){
            slot(object = object.pair[["projectedpca.1"]], name = "cell.embeddings") <- Sweep(
              x = Embeddings(object = object.pair[["projectedpca.1"]]),
              MARGIN = 2,
              STATS = apply(X = Embeddings(object = object.pair[["projectedpca.1"]]), MARGIN = 2, FUN = sd),
              FUN = "/"
            )
            slot(object = object.pair[["projectedpca.2"]], name = "cell.embeddings") <- Sweep(
              x = Embeddings(object = object.pair[["projectedpca.2"]]),
              MARGIN = 2,
              STATS = apply(X = Embeddings(object = object.pair[["projectedpca.2"]]), MARGIN = 2, FUN = sd),
              FUN = "/"
            )
            object.pair <- L2Dim(object = object.pair, reduction = "projectedpca.1")
            object.pair <- L2Dim(object = object.pair, reduction = "projectedpca.2")
            reduction <- paste0(reduction, ".l2")
            reduction.2 <- paste0(reduction.2, ".l2")
          }
          object.pair
        },
        stop("Invalid reduction parameter. Please choose either cca or rpca")
      )
      
      internal.neighbors <- internal.neighbors[c(i, j)]
      
      anchors <- FindAnchors.wdist(
        object.pair = object.pair,
        assay = c("ToIntegrate", "ToIntegrate"),
        slot = slot,
        cells1 = colnames(x = object.1),
        cells2 = colnames(x = object.2),
        internal.neighbors = internal.neighbors,
        reduction = reduction,
        reduction.2 = reduction.2,
        nn.reduction = nn.reduction,
        dims = dims,
        k.anchor = k.anchor,
        k.filter = k.filter,
        k.score = k.score,
        max.features = max.features,
        nn.method = nn.method,
        eps = eps,
        verbose = verbose
      )
      anchors[, 1] <- anchors[, 1] + offsets[i]
      anchors[, 2] <- anchors[, 2] + offsets[j]
      
      return(anchors)
    }
  )
  
  all.anchors <- do.call(what = 'rbind', args = all.anchors)
  
  all.anchors <- rbind(all.anchors, all.anchors[, c(2, 1, 3, 5, 4)])  ##keep distance information
  
  all.anchors <- AddDatasetID.2(anchor.df = all.anchors, offsets = offsets, obj.lengths = objects.ncell)  ##add dataset IDs
  command <- LogSeuratCommand(object = object.list[[1]], return.command = TRUE)
  anchor.set <- new(Class = "IntegrationAnchorSet",
                    object.list = object.list,
                    reference.objects = reference %||% seq_along(object.list),
                    anchors = all.anchors,
                    offsets = offsets,
                    anchor.features = anchor.features,
                    command = command
  )
  return(anchor.set)
}

#First (otherwise second) function
`%||%` <- function(lhs, rhs) {
  if (!is.null(x = lhs)) {
    return(lhs)
  } else {
    return(rhs)
  }
}

##Modified Sweep function (from Seurat 3.1)
Sweep <- function(x, MARGIN, STATS, FUN = '-', check.margin = TRUE, ...) {
  if (any(grepl(pattern = 'X', x = names(x = formals(fun = sweep))))) {
    return(sweep(
      X = x,
      MARGIN = MARGIN,
      STATS = STATS,
      FUN = FUN,
      check.margin = check.margin,
      ...
    ))
  } else {
    return(sweep(
      x = x,
      MARGIN = MARGIN,
      STATS = STATS,
      FUN = FUN,
      check.margin = check.margin,
      ...
    ))
  }
}

# Find nearest neighbors
# IMPORTED from Seurat 3.2.1
FindNN.Seurat <- function(
    object,
    cells1 = NULL,
    cells2 = NULL,
    internal.neighbors,
    grouping.var = NULL,
    dims = 1:10,
    reduction = "cca.l2",
    reduction.2 = character(),
    nn.dims = dims,
    nn.reduction = reduction,
    k = 300,
    nn.method = "rann",
    eps = 0,
    integration.name = 'integrated',
    verbose = TRUE
) {
  if (xor(x = is.null(x = cells1), y = is.null(x = cells2))) {
    stop("cells1 and cells2 must both be specified")
  }
  if (!is.null(x = cells1) && !is.null(x = cells2) && !is.null(x = grouping.var)) {
    stop("Specify EITHER grouping.var or cells1/2.")
  }
  if (is.null(x = cells1) && is.null(x = cells2) && is.null(x = grouping.var)) {
    stop("Please set either cells1/2 or grouping.var")
  }
  if (!is.null(x = grouping.var)) {
    if (nrow(x = unique(x = object[[grouping.var]])) != 2) {
      stop("Number of groups in grouping.var not equal to 2.")
    }
    groups <- names(x = sort(x = table(object[[grouping.var]]), decreasing = TRUE))
    cells1 <- colnames(x = object)[object[[grouping.var]] == groups[[1]]]
    cells2 <- colnames(x = object)[object[[grouping.var]] == groups[[2]]]
  }
  if (verbose) {
    message("Finding neighborhoods")
  }
  if (!is.null(x = internal.neighbors[[1]])) {
    nnaa <- internal.neighbors[[1]]
    nnbb <- internal.neighbors[[2]]
  } else {
    dim.data.self <- Embeddings(object = object[[nn.reduction]])[ ,nn.dims]
    dims.cells1.self <- dim.data.self[cells1, ]
    dims.cells2.self <- dim.data.self[cells2, ]
    nnaa <- Seurat:::NNHelper(
      data = dims.cells1.self,
      k = k + 1,
      method = nn.method,
      eps = eps
    )
    nnbb <- Seurat:::NNHelper(
      data = dims.cells2.self,
      k = k + 1,
      method = nn.method,
      eps = eps
    )
  }
  if (length(x = reduction.2) > 0) {
    nnab <- Seurat:::NNHelper(
      data = Embeddings(object = object[[reduction.2]])[cells2, ],
      query = Embeddings(object = object[[reduction.2]])[cells1, ],
      k = k,
      method = nn.method,
      eps = eps
    )
    nnba <- Seurat:::NNHelper(
      data = Embeddings(object = object[[reduction]])[cells1, ],
      query = Embeddings(object = object[[reduction]])[cells2, ],
      k = k,
      method = nn.method,
      eps = eps
    )
  } else {
    dim.data.opposite <- Embeddings(object = object[[reduction]])[ ,dims]
    dims.cells1.opposite <- dim.data.opposite[cells1, ]
    dims.cells2.opposite <- dim.data.opposite[cells2, ]
    nnab <- Seurat:::NNHelper(
      data = dims.cells2.opposite,
      query = dims.cells1.opposite,
      k = k,
      method = nn.method,
      eps = eps
    )
    nnba <- Seurat:::NNHelper(
      data = dims.cells1.opposite,
      query = dims.cells2.opposite,
      k = k,
      method = nn.method,
      eps = eps
    )
  }
  
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'neighbors',
    new.data = list('nnaa' = nnaa, 'nnab' = nnab, 'nnba' = nnba, 'nnbb' = nnbb, 'cells1' = cells1, 'cells2' = cells2)
  )
  return(object)
}


# Find Anchor pairs
# IMPORTED from Seurat 3.2.1
FindAnchorPairs.Seurat <- function(
    object,
    integration.name = 'integrated',
    k.anchor = 5,
    verbose = TRUE
) {
  neighbors <- GetIntegrationData(object = object, integration.name = integration.name, slot = 'neighbors')
  max.nn <- c(ncol(x = neighbors$nnab), ncol(x = neighbors$nnba))
  if (any(k.anchor > max.nn)) {
    message(paste0('warning: requested k.anchor = ', k.anchor, ', only ', min(max.nn), ' in dataset'))
    k.anchor <- min(max.nn)
  }
  if (verbose) {
    message("Finding anchors")
  }
  # convert cell name to neighbor index
  nn.cells1 <- neighbors$cells1
  nn.cells2 <- neighbors$cells2
  cell1.index <-  suppressWarnings(which(colnames(x = object) == nn.cells1, arr.ind = TRUE))
  ncell <- 1:nrow(x = neighbors$nnab)
  ncell <- ncell[ncell %in% cell1.index]
  anchors <- list()
  # pre allocate vector
  anchors$cell1 <- rep(x = 0, length(x = ncell) * 5)
  anchors$cell2 <- anchors$cell1
  anchors$score <- anchors$cell1 + 1
  idx <- 0
  indices.ab <- Indices(object = neighbors$nnab)
  indices.ba <- Indices(object = neighbors$nnba)
  for (cell in ncell) {
    neighbors.ab <- indices.ab[cell, 1:k.anchor]
    mutual.neighbors <- which(
      x = indices.ba[neighbors.ab, 1:k.anchor, drop = FALSE] == cell,
      arr.ind = TRUE
    )[, 1]
    for (i in neighbors.ab[mutual.neighbors]){
      idx <- idx + 1
      anchors$cell1[idx] <- cell
      anchors$cell2[idx] <- i
      anchors$score[idx] <- 1
    }
  }
  anchors$cell1 <- anchors$cell1[1:idx]
  anchors$cell2 <- anchors$cell2[1:idx]
  anchors$score <- anchors$score[1:idx]
  anchors <- t(x = do.call(what = rbind, args = anchors))
  anchors <- as.matrix(x = anchors)
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'anchors',
    new.data = anchors
  )
  if (verbose) {
    message(paste0("\tFound ", nrow(x = anchors), " anchors"))
  }
  return(object)
}

# Top dim features
# IMPORTED from Seurat 3.2.1
TopDimFeatures.Seurat <- function(
    object,
    reduction,
    dims = 1:10,
    features.per.dim = 100,
    max.features = 200,
    projected = FALSE
) {
  dim.reduction <- object[[reduction]]
  max.features <- max(length(x = dims) * 2, max.features)
  num.features <- sapply(X = 1:features.per.dim, FUN = function(y) {
    length(x = unique(x = as.vector(x = sapply(X = dims, FUN = function(x) {
      unlist(x = TopFeatures(object = dim.reduction, dim = x, nfeatures = y, balanced = TRUE, projected = projected))
    }))))
  })
  max.per.pc <- which.max(x = num.features[num.features < max.features])
  features <- unique(x = as.vector(x = sapply(X = dims, FUN = function(x) {
    unlist(x = TopFeatures(object = dim.reduction, dim = x, nfeatures = max.per.pc, balanced = TRUE, projected = projected))
  })))
  features <- unique(x = features)
  return(features)
}

# Filter anchors
# IMPORTED from Seurat 3.2.1
FilterAnchors.Seurat <- function(
    object,
    assay = NULL,
    slot = "data",
    integration.name = 'integrated',
    features = NULL,
    k.filter = 200,
    nn.method = "rann",
    eps = 0,
    verbose = TRUE
) {
  if (verbose) {
    message("Filtering anchors")
  }
  assay <- assay %||% DefaultAssay(object = object)
  features <- features %||% VariableFeatures(object = object)
  if (length(x = features) == 0) {
    stop("No features provided and no VariableFeatures computed.")
  }
  features <- unique(x = features)
  neighbors <- GetIntegrationData(object = object, integration.name = integration.name, slot = 'neighbors')
  nn.cells1 <- neighbors$cells1
  nn.cells2 <- neighbors$cells2
  cn.data1 <- L2Norm(
    mat = as.matrix(x = t(x = GetAssayData(
      object = object[[assay[1]]],
      slot = slot)[features, nn.cells1])),
    MARGIN = 1)
  cn.data2 <- L2Norm(
    mat = as.matrix(x = t(x = GetAssayData(
      object = object[[assay[2]]],
      slot = slot)[features, nn.cells2])),
    MARGIN = 1)
  nn <- NNHelper(
    data = cn.data2[nn.cells2, ],
    query = cn.data1[nn.cells1, ],
    k = k.filter,
    method = nn.method,
    eps = eps
  )
  
  anchors <- GetIntegrationData(object = object, integration.name = integration.name, slot = "anchors")
  position <- sapply(X = 1:nrow(x = anchors), FUN = function(x) {
    which(x = anchors[x, "cell2"] == Indices(object = nn)[anchors[x, "cell1"], ])[1]
  })
  anchors <- anchors[!is.na(x = position), ]
  if (verbose) {
    message("\tRetained ", nrow(x = anchors), " anchors")
  }
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = "anchors",
    new.data = anchors
  )
  return(object)
}

# Score Anchors
# IMPORTED from Seurat 3.2.1
ScoreAnchors.Seurat <- function(
    object,
    assay = NULL,
    integration.name = 'integrated',
    verbose = TRUE,
    k.score = 30,
    do.cpp = TRUE
) {
  assay <- assay %||% DefaultAssay(object = object)
  anchor.df <- as.data.frame(x = GetIntegrationData(object = object, integration.name = integration.name, slot = 'anchors'))
  neighbors <- GetIntegrationData(object = object, integration.name = integration.name, slot = "neighbors")
  offset <- length(x = neighbors$cells1)
  indices.aa <- Indices(object = neighbors$nnaa)
  indices.bb <- Indices(object = neighbors$nnbb)
  indices.ab <- Indices(object = neighbors$nnab)
  indices.ba <- Indices(object = neighbors$nnba)
  nbrsetA <- function(x) c(indices.aa[x, 1:k.score], indices.ab[x, 1:k.score] + offset)
  nbrsetB <- function(x) c(indices.ba[x, 1:k.score], indices.bb[x, 1:k.score] + offset)
  # score = number of shared neighbors
  anchor.new <- data.frame(
    'cell1' = anchor.df[, 1],
    'cell2' = anchor.df[, 2],
    'score' = mapply(
      FUN = function(x, y) {
        length(x = intersect(x = nbrsetA(x = x), nbrsetB(x = y)))},
      anchor.df[, 1],
      anchor.df[, 2]
    )
  )
  # normalize the score
  max.score <- quantile(anchor.new$score, 0.9)
  min.score <- quantile(anchor.new$score, 0.01)
  anchor.new$score <- anchor.new$score - min.score
  anchor.new$score <- anchor.new$score / (max.score - min.score)
  anchor.new$score[anchor.new$score > 1] <-  1
  anchor.new$score[anchor.new$score < 0] <- 0
  anchor.new <- as.matrix(x = anchor.new)
  object <- SetIntegrationData(
    object = object,
    integration.name = integration.name,
    slot = 'anchors',
    new.data = anchor.new
  )
  return(object)
}

###Find Anchors functions, modified from Seurat 3.1 to return the distance in PC space between pairs of datasets
FindAnchors.wdist <- function(
    object.pair,
    assay,
    slot,
    cells1,
    cells2,
    internal.neighbors,
    reduction,
    reduction.2 = character(),
    nn.reduction = reduction,
    dims = 1:10,
    k.anchor = 5,
    k.filter = 200,
    k.score = 30,
    max.features = 200,
    nn.method = "rann",
    eps = 0,
    projected = FALSE,
    verbose = TRUE
) {
  # compute local neighborhoods, use max of k.anchor and k.score if also scoring to avoid
  # recomputing neighborhoods
  k.neighbor <- k.anchor
  if (!is.na(x = k.score)) {
    k.neighbor <- max(k.anchor, k.score)
  }
  object.pair <- FindNN.Seurat(
    object = object.pair,
    cells1 = cells1,
    cells2 = cells2,
    internal.neighbors = internal.neighbors,
    dims = dims,
    reduction = reduction,
    reduction.2 = reduction.2,
    nn.reduction = nn.reduction,
    k = k.neighbor,
    nn.method = nn.method,
    eps = eps,
    verbose = verbose
  )
  
  object.pair <- FindAnchorPairs.Seurat(
    object = object.pair,
    integration.name = "integrated",
    k.anchor = k.anchor,
    verbose = verbose
  )
  
  if (!is.na(x = k.filter)) {
    top.features <- TopDimFeatures.Seurat(
      object = object.pair,
      reduction = reduction,
      dims = dims,
      features.per.dim = 100,
      max.features = max.features,
      projected = projected
    )
    object.pair <- FilterAnchors.Seurat(
      object = object.pair,
      assay = assay,
      slot = slot,
      integration.name = 'integrated',
      features = top.features,
      k.filter = k.filter,
      nn.method = nn.method,
      eps = eps,
      verbose = verbose
    )
  }
  if (!is.na(x = k.score)) {
    object.pair = ScoreAnchors.Seurat(
      object = object.pair,
      assay = DefaultAssay(object = object.pair),
      integration.name = "integrated",
      verbose = verbose,
      k.score = k.score
    )
  }
  
  ###Return distances
  anc.tab <- object.pair@tools$integrated@anchors
  d1.2 <- numeric(length = dim(anc.tab)[1])
  d2.1 <- numeric(length = dim(anc.tab)[1])
  for (r in 1:dim(anc.tab)[1]) {
    c1 <- anc.tab[r,"cell1"]
    c2 <- anc.tab[r,"cell2"]
    d1.2[r] <- object.pair@tools$integrated@neighbors$nnab@nn.dist[c1, which(object.pair@tools$integrated@neighbors$nnab@nn.idx[c1,] == c2 )]
    d2.1[r] <- object.pair@tools$integrated@neighbors$nnba@nn.dist[c2, which(object.pair@tools$integrated@neighbors$nnba@nn.idx[c2,] == c1 )]
  }
  
  object.pair@tools$integrated@anchors <- cbind(object.pair@tools$integrated@anchors, dist1.2=d1.2)
  object.pair@tools$integrated@anchors <- cbind(object.pair@tools$integrated@anchors, dist2.1=d2.1)
  
  anchors <- GetIntegrationData(
    object = object.pair,
    integration.name = 'integrated',
    slot = 'anchors'
  )
  
  return(anchors)
}

AddDatasetID.2 <- function(
    anchor.df,
    offsets,
    obj.lengths
) {
  ndataset <- length(x = offsets)
  total.cells <- sum(obj.lengths)
  offsets <- c(offsets, total.cells)
  row.offset <- rep.int(x = offsets[1:ndataset], times = obj.lengths)
  dataset <- rep.int(x = 1:ndataset, times = obj.lengths)
  
  anchor.df <- data.frame(
    'cell1' = anchor.df[, 1] - row.offset[anchor.df[, 1]],
    'cell2' = anchor.df[, 2] - row.offset[anchor.df[, 2]],
    'score' = anchor.df[, 3],
    'dataset1' = dataset[anchor.df[, 1]],
    'dataset2' = dataset[anchor.df[, 2]],
    'dist1.2' = anchor.df[, 4],
    'dist2.1' = anchor.df[, 5]
  )
  return(anchor.df)
}
