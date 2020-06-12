Run.STACAS <- function (
  object.list = NULL,
  assay = NULL,
  reference = NULL,
  anchor.features = 500,
  dims = 1:10,
  dist.thr = NULL,
  dist.pct = 0.8,
  k.anchor = 5,
  k.score = 30,
  plot.file = "anchor.dist.mean.pairwise.png",
  verbose = TRUE
) {
  ref.anchors <- FindAnchors.STACAS(object.list, dims=dims, k.anchor=k.anchor, k.score=k.score, assay=assay,
                                    reference=reference, anchor.features=anchor.features, verbose=verbose)

  if (is.null(names(object.list))) {
    names <- 1:length(object.list)
  } else {
    names <- names(object.list)
  }

  plots <- PlotAnchors.STACAS(ref.anchors, obj.names=names, dist.thr=dist.thr,  dist.pct = dist.pct)

  g.cols <- 3
  g.rows <- as.integer((length(plots)+2)/g.cols)
  g <- do.call("arrangeGrob", c(plots, ncol=g.cols, nrow=g.rows))
  ggsave(plot.file, plot=g, width = 6*g.cols, height = 6*g.rows)

  plot.new()
  grid.draw(g)

  ref.anchors.filtered <- FilterAnchors.STACAS(ref.anchors, dist.thr=dist.thr, dist.pct = dist.pct)

  return(ref.anchors.filtered)
}

FindAnchors.STACAS <- function (
  object.list = NULL,
  assay = NULL,
  reference = NULL,
  anchor.features = 500,
  dims = 1:10,
  k.anchor = 5,
  k.score = 30,
  verbose = TRUE
) {

  scale = FALSE
  reduction = "rpca"

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

  #anchor features
  if (is.numeric(x = anchor.features)) {
    n.this <- anchor.features
    if (verbose) {
      message("Computing ", anchor.features, " integration features")
    }
    anchor.features <- SelectIntegrationFeatures(
      object.list = object.list,
      nfeatures = n.this*2,
      assay = assay
    )
    #Remove cell cycling genes
    anchor.features <- head(setdiff(anchor.features, cycling.genes.STACAS), n.this)
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

  #perform rPCA to find anchors
  if (is.null(reference)) {
    ref.anchors <- FindIntegrationAnchors.wdist(object.list, dims = dims, k.anchor = k.anchor, anchor.features=anchor.features,
                                                scale=scale, reduction=reduction, assay=assay, k.score=k.score, verbose=verbose)
  } else {
    ref.anchors <- FindIntegrationAnchors.wdist(object.list, reference=reference, dims = dims, k.anchor = k.anchor, anchor.features=anchor.features,
                                                scale=scale, reduction=reduction, assay=assay, k.score=k.score, verbose=verbose)
  }
  for (r in 1:dim(ref.anchors@anchors)[1]) {
    ref.anchors@anchors[r,"dist.mean"] <- mean(c(ref.anchors@anchors[r,"dist1.2"],ref.anchors@anchors[r,"dist2.1"]))
    ref.anchors@anchors[r,"dist.max"] <- max(c(ref.anchors@anchors[r,"dist1.2"],ref.anchors@anchors[r,"dist2.1"]))
    ref.anchors@anchors[r,"dist.min"] <- min(c(ref.anchors@anchors[r,"dist1.2"],ref.anchors@anchors[r,"dist2.1"]))
  }

  return(ref.anchors)
}

FilterAnchors.STACAS <- function(
  ref.anchors = NULL,
  dist.thr = NULL,
  dist.pct = 0.8
) {

  if(dist.pct<0 | dist.pct>1) {
    stop("Variable dist.pct must be a real number between 0 and 1")
  }

  #determine empirically the cutoff threshold, if not provided
  if(is.null(dist.thr)) {
    i=1
    min.pcntile <- vector()
    for (l1 in levels(as.factor(ref.anchors@anchors$dataset1))) {
      j=1
      pcntile <- vector()
      for (l2 in levels(as.factor(ref.anchors@anchors$dataset2))) {
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
  }
  message(sprintf("Filter anchors using distance threshold t=%.3f",dist.thr))
  ref.anchors@anchors <- subset(ref.anchors@anchors, subset=dist.mean<=dist.thr)
  return(ref.anchors)
}

PlotAnchors.STACAS <- function(
  ref.anchors = NULL,
  obj.names = NULL,
  dist.thr = NULL,
  dist.pct = 0.8
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

  #determine empirically the cutoff threshold, if not provided
  if(is.null(dist.thr)) {
    i=1
    min.pcntile <- vector()
    for (l1 in levels(as.factor(ref.anchors@anchors$dataset1))) {
      j=1
      pcntile <- vector()
      for (l2 in levels(as.factor(ref.anchors@anchors$dataset2))) {
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
  }
  ###Make distribution plots
  anchortab.toprint <- anchortab[]

  a.tmp <- anchortab.toprint
  for (r in 1:dim(anchortab.toprint)[1]) {
    anchortab.toprint[r,"dataset1"] <- obj.names[a.tmp[r,"dataset1"]]
    anchortab.toprint[r,"dataset2"] <- obj.names[a.tmp[r,"dataset2"]]
  }
  rm(a.tmp)

  #my.colors=palette(rainbow(length(levs)))
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

SampleTree.STACAS <- function (
  anchorset
) {
  object.list <- slot(object = anchorset, name = "object.list")
  reference.objects <- slot(object = anchorset, name = "reference.objects")
  anchors <- slot(object = anchorset, name = "anchors")
  objects.ncell <- sapply(X = object.list, FUN = ncol)
  offsets <- slot(object = anchorset, name = "offsets")

  similarity.matrix <- Seurat:::CountAnchors(
    anchor.df = anchors,
    offsets = offsets,
    obj.lengths = objects.ncell
  )

  similarity.matrix <- similarity.matrix[reference.objects, reference.objects]
  sample.tree <- Seurat:::BuildSampleTree(similarity.matrix = similarity.matrix)
  sample.tree <- Seurat:::AdjustSampleTree(x = sample.tree, reference.objects = reference.objects)

  #precalculate anchors between sets
  nanch <- list()
  names(x = object.list) <- as.character(-(1:length(x = object.list)))
  for (i in 1:length(object.list)) {
    nanch[[as.character(-i)]] <- sum(anchorset@anchors$dataset1==i)
  }

  for (r in 1:nrow(sample.tree)) {
    pair <- sample.tree[r, ]
    length1 <- nanch[[as.character(pair[1])]]
    length2 <- nanch[[as.character(pair[2])]]

    if (length2 > length1) {
      pair <- rev(pair)
      sample.tree[r, ] <- pair
    }

    nanch[[as.character(r)]] <- length1+length2
  }
  return(sample.tree)
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
  anchor.set <- new(Class = "AnchorSet",
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

##Modified Sweep function
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

###Find Anchors functions, modified to return the distance in PC space between pairs of datasets
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
  object.pair <- Seurat:::FindNN(
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

  object.pair <- Seurat:::FindAnchorPairs(
    object = object.pair,
    integration.name = "integrated",
    k.anchor = k.anchor,
    verbose = verbose
  )

  if (!is.na(x = k.filter)) {
    top.features <- Seurat:::TopDimFeatures(
      object = object.pair,
      reduction = reduction,
      dims = dims,
      features.per.dim = 100,
      max.features = max.features,
      projected = projected
    )
    object.pair <- Seurat:::FilterAnchors(
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
    object.pair = Seurat:::ScoreAnchors(
      object = object.pair,
      assay = DefaultAssay(object = object.pair),
      integration.name = "integrated",
      verbose = verbose,
      k.score = k.score
    )
  }

  ###HERE
  anc.tab <- object.pair@tools$integrated@anchors
  d1.2 <- numeric(length = dim(anc.tab)[1])
  d2.1 <- numeric(length = dim(anc.tab)[1])
  for (r in 1:dim(anc.tab)[1]) {
    c1 <- anc.tab[r,"cell1"]
    c2 <- anc.tab[r,"cell2"]
    d1.2[r] <- object.pair@tools$integrated@neighbors$nnab$nn.dists[c1, which(object.pair@tools$integrated@neighbors$nnab$nn.idx[c1,] == c2 )]
    d2.1[r] <- object.pair@tools$integrated@neighbors$nnba$nn.dists[c2, which(object.pair@tools$integrated@neighbors$nnba$nn.idx[c2,] == c1 )]
  }

  object.pair@tools$integrated@anchors <- cbind(object.pair@tools$integrated@anchors, dist1.2=d1.2)
  object.pair@tools$integrated@anchors <- cbind(object.pair@tools$integrated@anchors, dist2.1=d2.1)

  ##TO HERE

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
