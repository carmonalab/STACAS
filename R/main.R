#alpha = 1: distance based score (no knn is considered)
#alpha = 0: knn based score
#alpha = 0.5 equal weight for distance and knn (dot product) for score computing. 
# beta = 1: all inconsistent anchors will be removed
# beta = 0.5 allows some anchors randomly be accepted because its score is above a given quantile threshold of "non-inconsistent anchor score distribution" (q_boltzman)
# beta = 0: all inconsistent anchors having a score above the given quantile of non-inconsistent anchor score distribution, will be accepted. 
# beta = 0, q_boltzman = 0: all anchors will be accepted either inconsistent or not (not recommended)
# dist.hard.threshold is intended to reproduce similar behaviour to original STACAS algorithm (but avoiding Errors because lack of anchors)  ## ToDo a warning must be set here when few anchors survive
scoring_anchors_by_distance <- function(anchors, alpha = 0.5, remove_inconsistent_anchors = T, dist.hard.thr = NULL,
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


## This is actually not recovering the STACAS precedent behaviour when dist.hard.thr is set (REVISE AND MODIFY)
SampleTree.weighted <- function (
  anchorset,
  hclust.method = "complete",
  usecol = "score",
  dist.hard.thr = NULL,
  method = "cum.sum"
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
  plot(hclust(d = distance.matrix,method = hclust.method))
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


# Count anchors between data sets
# IMPORTED from Seurat 3.2.1
weighted.Anchors.STACAS <- function(
  anchor.df,
  offsets,
  obj.lengths,
  usecol = "score",
  method = "cum.sum"
) {
  usecol = match.arg(arg = usecol ,choices = c("dist.mean","score"))
  method = match.arg(arg = method ,choices = c("cum.sum","counts"))
  
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



