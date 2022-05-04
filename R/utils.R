get_all_genes <- function(obj.list){
  all.genes <- row.names(obj.list[[1]])
  for (i in 2:length(obj.list)) {
    all.genes <- intersect(all.genes, row.names(obj.list[[i]]))
  }
  return(all.genes)
}

# Logistic transformation
logistic <- function(x,scale = NULL, center = NULL, invert = F){
  a <- ifelse(invert,1, -1)
  if(is.null(scale))  scale = sd(x)
  if(is.null(center)) center = mean(x)
  sigm <- 1/(1+ exp(a*(x-center)/scale))
  return(sigm)    
}

#boltzmann_based_rejection of inconsistent anchors
# q: quantile score of non-inconsistent anchors
# Beta = 1: q is completely ignored. This represent a maximum confidence in our labeling. None inconsistent anchor will be accepted.
# Beta = 0: deterministic theta quantile behavior:  All inconsistent anchors with scoring higher than "q", will be accepted.
# Beta in [0-1] range: Boltzmann based acceptation behavior. Only anchors with higher score than "q" could be accepted with a probability that increases with their score.  
# Beta = 0, q = 0: Completely ignoring of labeling: All anchors will be accepted, either being inconsistent or not
# Default: theta = 0.8; beta = 0.5
boltzmann_based_rejection <- function(scoring,inconsistent_flag,beta = 0.5, q = 0.8){
  e0 <- quantile(scoring[!as.logical(inconsistent_flag)],q)  # scale is given by q75 score in NON inconsistent anchors
  DE <- scoring[as.logical(inconsistent_flag)] - e0             #Delta Energy is given by the actual score of inconsistent anchors.   
  kT = quantile(DE[DE>0],beta)
  rejection_proba <- sapply(exp(-DE/kT),function(x){min(1,x)})
  #plot(E+e0,rejection_proba,xlab = "score")
  rejection_inconsist <- runif(length(DE),0,1) < rejection_proba**(1-beta)
  rejection <- rep(F,length(scoring))
  
  # this if is necesary because b = 0 and q = 0 could yet reject a few anchors (those inconsistent ones having higher score than the minumum of the non-inconsistent distribution)
  if(beta != 0 & q != 0){
    rejection[as.logical(inconsistent_flag)] <- rejection_inconsist
  }
  #table(rejection)  
  return(rejection)
}


# to be used in the distance matrix among datasets (when building SampleTree intregration)
strength_function <- function(anchor.df,
                              method = "cum.sum", 
                              usecol = "dist.mean"
                              ){
  method = match.arg(arg = method ,choices = c("counts","mean.score","cum.sum"))
  usecol = match.arg(arg = usecol ,choices = c("dist.mean","score"))
  
  if(method== "counts") {
    strength <- anchor.df%>%nrow()
  }
  
  if(method== "mean.score") {
    if(usecol == "dist.mean")    strength <- 1- anchor.df[,usecol]%>%mean()
    if(usecol == "score")    strength <- anchor.df[,usecol]%>%mean()
    
  }

  if(method== "cum.sum") {
    if(usecol == "dist.mean")  strength <- (1-anchor.df[,usecol])%>%sum()
    if(usecol == "score")  strength <- (anchor.df[,usecol])%>%sum()
    
  }
  
  return(strength)    
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
