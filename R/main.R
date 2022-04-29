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



