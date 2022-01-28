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



## Radar plot for comparing 2 or more clusters in the same dataset. 
RadarPlot <- function(obj,meta.data.column, genes4radar = NULL, assay = "RNA",slot = "data"){
  
  if(is.null(genes4radar)){  
    genes4radar <- c("Cd4", "Cd8a", "Tcf7", "Ccr7", 
                     "Gzmb", "Gzmk", "Pdcd1", "Tox", "Mki67")
  }
  useData <- GetAssayData(obj, assay = assay, slot = slot)
  genes4radar <- intersect(genes4radar, row.names(useData))
  genes4radar <- sort(genes4radar)
  
  df <- data.frame(list("Gene"= genes4radar))
  rownames(df) <- df$Gene
  df2 <- merge(df, useData, by = 0);
  df2 <- df2[,-1]
  df3 <- df2%>%melt(id.vars = "Gene")
  colnames(df3)[2:3] <- c("cell","expression")
  
  aux = obj@meta.data[,meta.data.column,drop =F]
  aux$cell <- aux%>%rownames()
  
  ncolors = length(unique(aux[,meta.data.column]))
  radar.colors <- c("black", hue_pal()(ncolors - 1))
  
  df4 <- merge(df3,aux,by = "cell")
  df4$cluster <- df4[,meta.data.column]
  
  df5 <- df4%>%group_by(Gene,cluster) %>% summarise(Expression = mean(expression))
  
  # Now, radar plot
  ymin = min(df5$Expression)
  ymax = max(df5$Expression)
  #                 ,"Expression","cluster")
  plt <- ggplot(data = df5, aes(x = Gene, y = Expression, 
                                group = cluster, colour = cluster, fill = cluster)) + 
    geom_point(size = 2) + geom_polygon(size = 0.75, 
                                        alpha = 0.1) + ylim(ymin, ymax) + scale_x_discrete()  +
    scale_fill_manual(values = radar.colors) + scale_colour_manual(values = radar.colors) +
    theme_light() + theme(axis.text.x = element_blank()) + 
    annotate(geom = "text", x = seq(1, length(genes4radar))*0.98, 
             y = ymax - 0.05  , label = genes4radar, 
             size = 3) + 
    coord_polar()
  return(plt)
}


# Deprecated function
inconsistent_anchors <- function(anchors,pred.colname,scoring = F){
  labels <- data.frame()
  if(scoring){
    score.colname <- sprintf("score_%s",pred.colname)
    for(i in anchors@reference.objects){
      x <-anchors@object.list[[i]]
      labels <- rbind(labels,data.frame(dataset = i, ncell = 1:ncol(x), label = x[[pred.colname]], scoring = x[[score.colname]]))
    }
  }else{
    for(i in anchors@reference.objects){
      x <-anchors@object.list[[i]]
      labels <- rbind(labels,data.frame(dataset = i, ncell = 1:ncol(x), label = x[[pred.colname]]))
    }
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


