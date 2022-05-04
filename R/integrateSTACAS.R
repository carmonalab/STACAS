Run.STACAS <- function (
  object.list = NULL,
  assay = NULL,
  reference = NULL,
  anchor.features = 500,
  dims = 1:10,
  normalization.method = c("LogNormalize", "SCT"),
  dist.thr = NULL,
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
  normalization.method = c("LogNormalize", "SCT"),
  k.anchor = 5,
  k.score = 30,
  verbose = TRUE
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
                                                normalization.method = normalization.method,
                                                scale=scale, reduction=reduction, assay=assay, k.score=k.score, verbose=verbose)
  } else {
    ref.anchors <- FindIntegrationAnchors.wdist(object.list, reference=reference, dims = dims, k.anchor = k.anchor, anchor.features=anchor.features,
                                                normalization.method = normalization.method,
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
  
  min.anchors.warn = 100
  
  if(dist.pct<0 | dist.pct>1) {
    stop("Variable dist.pct must be a real number between 0 and 1")
  }
  
  lev1 <- levels(as.factor(ref.anchors@anchors$dataset1))
  lev2 <- levels(as.factor(ref.anchors@anchors$dataset2))
  
  #determine empirically the cutoff threshold, if not provided
  if(is.null(dist.thr)) {
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
  }
  message(sprintf("Filter anchors using distance threshold t=%.3f",dist.thr))
  ref.anchors@anchors <- subset(ref.anchors@anchors, subset=dist.mean<=dist.thr)
  
  a1 <- factor(ref.anchors@anchors[,"dataset1"], levels=lev1)
  a2 <- factor(ref.anchors@anchors[,"dataset2"], levels=lev2)
  
  anchor.tab <- table(a1, a2)
  
  max.anchors <- apply(anchor.tab, 1, max)
  
  for (i in seq_along(max.anchors)) {
    if (max.anchors[i] == 0) {
      text <- sprintf("Warning! No anchors left after filtering for dataset %i - integration will fail\n", i)
      text <- paste0(text, "You may want to specify a higher dist.thr parameter")
      warning(text)
    }
    else if (max.anchors[i] < min.anchors.warn) {
      text <- sprintf("Warning! Fewer than %i anchors left after filtering for dataset %i - integration may fail\n", min.anchors.warn, i)
      text <- paste0(text, "You may want to specify a higher dist.thr parameter")
      warning(text)
    }
  }  
  
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

#' Anchor scoring function 
#'
#' Score anchors by mixing MNN based score and PCA based anchor distance
#'
#' @param anchors anchor object obtained from \code{inconsistent_anchors} STACAS function. - Scoring will be applied to this object
#' @param alpha A tabular model with scGate signatures. See Details for this format
#' @param remove_inconsistent_anchors 
#' @param dist.hard.thr hard threshold to emulate the classical STACAS behavior
#' @param dist.mean.center center of logistic function for transforming PCA distance
#' @param dist.mean.scale scale of logistic function for transforming PCA distance
#' @param q.dist.mean.center center of logistic function based on quantile value of PCA  distance distribution
#' @param dist.scale.factor scale factor of logistic function based on standard deviation of PCA distance distribution  
#' @param beta this parameter controls the way in which inconsistent anchors will be penalized. 
#' @param q_boltzmann this parameter defines a a baseline threshold to decide if reject an inconsistent anchor or not.
#' @param plot logical indicating if scoring must be plotted as a function of PCA distance

#' @return A new anchor object with modified anchor's scoring.
#' @details beta = 1: all inconsistent anchors will be removed. 
#' beta = 0 (q_boltzmann != 0)  all inconsistent anchors having a score above the given quantile of non-inconsistent anchor score distribution, will be accepted.
#' In the interval (0,1), this parameter allows inconsistent anchors to be randomly accepted with a boltzman distribution probability which scale is controled by beta. 
#' The higher beta, the lower probability or rejecting an inconsistent anchor.  
#' beta = 0 (q_boltzmann = 0) all anchors will be accepted either inconsistent or not (not recommended)
#' alpha = 1: distance based score (no knn is considered)
#' alpha = 0: knn based score
#' alpha = 0.5 equal weight for distance and knn (dot product) for score computing. 
# 
#' beta = 0.5 allows some anchors randomly be accepted because its score is above a given quantile threshold of "non-inconsistent anchor score distribution" (q_boltzman)
#' beta = 0:  
#' beta = 0, q_boltzman = 0: all anchors will be accepted either inconsistent or not (not recommended)
#' dist.hard.threshold is intended to reproduce similar behaviour to original STACAS algorithm (but avoiding Errors because lack of anchors)  ## ToDo a warning must be set here when few anchors survive

#' @examples
#' library(STACAS)
#' @export
scoring_anchors <- function(anchors, alpha = 0.5, remove_inconsistent_anchors = T, dist.hard.thr = NULL,
                                        dist.mean.center = NULL,dist.mean.scale = NULL, q.dist.mean.center = NULL, dist_scale_factor = 1,
                                        beta = 0.5, q_boltzmann = 0.8, plot =F){
  df <- anchors@anchors
  seurat_knn_score <- df$score
  
  
  # distance based scores
  distance_based_scores = T
  if(distance_based_scores){
    if(is.null(dist.mean.center)){
      dist.mean.center = mean(df$dist.mean,na.rm = T)
    }
    if (!is.null(q.dist.mean.center)){
      dist.mean.center = quantile(df$dist.mean,q.dist.mean.center)  
    }
    
    if(is.null(dist.mean.scale)){
      dist.mean.scale = sd(df$dist.mean,na.rm = T)/dist_scale_factor
    }else if (is.character(dist.mean.scale)){
      dist.mean.scale = as.numeric(dist.mean.scale) * sd(df$dist.mean,na.rm = T)  
    }
  }
  
  if(is.null(dist.hard.thr)){
    stacas_scoring <-  logistic(x = df$dist.mean,invert = T, center = dist.mean.center, scale = dist.mean.scale)**(alpha)  * (seurat_knn_score)**(1-alpha) 
    stacas_scoring[stacas_scoring == 0 ] <- min(stacas_scoring[stacas_scoring>0],na.rm = T)* 0.001  
  }else{
    to.be.killed <- df$dist.mean > dist.hard.thr
    stacas_scoring <- seurat_knn_score  # because STACAS algorithm uses this score as integration score. 
    #stacas_scoring <-  1-df$dist.mean 
    min.to.impute <- min(stacas_scoring[stacas_scoring>0],na.rm = T)* 0.001
    stacas_scoring[to.be.killed] <- min.to.impute  
  }
  
  reject_anchors <- boltzmann_based_rejection(scoring = stacas_scoring,inconsistent_flag = df$flag, beta = beta,q = q_boltzmann)
  if(remove_inconsistent_anchors){
    stacas_scoring[reject_anchors] <- min(stacas_scoring,na.rm = T)* 0.001  
  }
  
  df$score <- stacas_scoring
  anchors@anchors <- df
  if(plot) {
    message(sprintf("dist.mean.center = %s; dist.mean.scale = %s",dist.mean.center,dist.mean.scale))
    oo = order(df$score)
    plot(df$dist.mean[oo],df$score[oo], type = "l")
  }
  return(anchors)
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
#' @details 
#' @examples
#' library(STACAS)
#' @export

SampleTree.weighted <- function (
  anchorset,
  hclust.method = "complete",
  usecol = "score",
  dist.hard.thr = NULL,
  method = "cum.sum",
  plot = T
) {
  object.list <- slot(object = anchorset, name = "object.list")
  reference.objects <- slot(object = anchorset, name = "reference.objects")
  anchors <- slot(object = anchorset, name = "anchors")
  objects.ncell <- sapply(X = object.list, FUN = ncol)
  offsets <- slot(object = anchorset, name = "offsets")
  usecol = match.arg(arg = usecol ,choices = c("dist.mean","score"))
  method = match.arg(arg = method ,choices = c("cum.sum","counts"))
  
  if(!is.null(dist.hard.thr)){
    anchors <- anchors%>%subset(dist.mean < dist.hard.thr)
    usecol = "dist.mean"
    method = "counts"  ## Force count anchors instead of weighted sum  to mimic original STACAS behaviour
  }  
  ## Main change in the way of computing similarity matrix
  similarity.matrix <- weighted.Anchors.STACAS(
    anchor.df = anchors,
    offsets = offsets,
    obj.lengths = objects.ncell,
    usecol = usecol,
    method = method
  )
  
  similarity.matrix <- similarity.matrix[reference.objects, reference.objects]
  similarity.matrix <- similarity.matrix+10^-6 #avoid 1/0
  
  #minor change: use 1- simil instead of 1/simil
  #distance.matrix <- as.dist(m = 1 / similarity.matrix)
  distance.matrix <- as.dist(m = 1 - similarity.matrix)
  if(plot){
    plot(hclust(d = distance.matrix,method = hclust.method))
  }
  sample.tree <- hclust(d = distance.matrix,method = hclust.method)$merge
  sample.tree <- AdjustSampleTree.Seurat(x = sample.tree, reference.objects = reference.objects)
  
  #precalculate anchors between sets
  nanch <- list()
  names(x = object.list) <- as.character(-(1:length(x = object.list)))
  for (i in 1:length(object.list)) {
    #    nanch[[as.character(-i)]] <- sum(anchorset@anchors$dataset1==i)
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
    
    nanch[[as.character(r)]] <- length1+length2
  }
  return(sample.tree)
}


#' Mark anchor as inconsistent based on cell-type labeling 
#'
#' Anchors will be consider as inconsistent if they are linking cells classified as different cell populations. 
#'
#' @param anchors anchor object obtained from STACAS function \code{FindAnchors.STACAS} or the Seurat function \code{FindAnchors}. 
#' @param pred.colname colname of seurat object containing consistent celltype information across all samples/datasets 
#' @param scoring
#' 
#' @return A n integration tree to be passed to the integration function.
#' @details 
#' @examples
#' library(STACAS)
#' @export

inconsistent_anchors <- function(anchors,pred.colname){
  labels <- data.frame()

  for(i in anchors@reference.objects){
    x <-anchors@object.list[[i]]
    labels <- rbind(labels,data.frame(dataset = i, ncell = 1:ncol(x), label = x[[pred.colname]]))
  }
  
  
  not_assigned_flag <- is.na(labels[,pred.colname])| (labels[,pred.colname]== "NaN")
  labels[not_assigned_flag,pred.colname] <- "not.assigned"
  #message(labels%>%head())
  
  df.anch <- anchors@anchors
  merged = merge(df.anch,labels,by.x = c("dataset1","cell1"),by.y = c("dataset","ncell")) 
  merged = merge(merged,labels,by.x = c("dataset2","cell2"),by.y = c("dataset","ncell"), suffixes = c("_1","_2")) 
  
  #not_assigned <- merged[,c(paste0(pred.colname,"_1"),paste0(pred.colname,"_2"))] ==  "not.assigned"
  #no_info = apply(not_assigned,1,any)
  
  match_ok <- (merged[,paste0(pred.colname,"_1")] == merged[,paste0(pred.colname,"_2")])  
  label_asigned <- (merged[,paste0(pred.colname,"_1")] != "not.assigned") & (merged[,paste0(pred.colname,"_2")] != "not.assigned")
  merged[merged == "not.assigned"] <- NaN
  
  merged$flag <- (!match_ok & label_asigned) + 0
  
  col.order = colnames(anchors@anchors)
  extra_cols <- setdiff(colnames(merged), col.order)
  
  anchors@anchors <- merged[,c(col.order,extra_cols)]  # keep only consistent anchors (THIS GET BROKEN BECAUSE FEW ANCHORS AMONG SOME DATASET PAIRS)
  return(anchors)
}



SampleTree.STACAS <- function (
  anchorset
) {
  object.list <- slot(object = anchorset, name = "object.list")
  reference.objects <- slot(object = anchorset, name = "reference.objects")
  anchors <- slot(object = anchorset, name = "anchors")
  objects.ncell <- sapply(X = object.list, FUN = ncol)
  offsets <- slot(object = anchorset, name = "offsets")
  
  similarity.matrix <- CountAnchors.Seurat(
    anchor.df = anchors,
    offsets = offsets,
    obj.lengths = objects.ncell
  )
  
  similarity.matrix <- similarity.matrix[reference.objects, reference.objects]
  similarity.matrix <- similarity.matrix+10^-6 #avoid 1/0
  
  distance.matrix <- as.dist(m = 1 / similarity.matrix)
  sample.tree <- hclust(d = distance.matrix)$merge
  sample.tree <- AdjustSampleTree.Seurat(x = sample.tree, reference.objects = reference.objects)
  
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
