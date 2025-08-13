run_prerds<-function(data,minFeature=0,
                     maxFeature=NULL,
                     sample="sample",
                     minCount=0,
                     maxCount=NULL,
                     tissue=NULL){

data$groups<- "0"
obj<-data

#' filter object based on nCount and nFeature
if (is.null(maxCount)){
  maxCount <- max(obj$nCount_Spatial)
}
if (is.null(maxFeature)){
  maxFeature <- max(obj$nFeature_Spatial)
}

maxva<-max(obj$nFeature_Spatial)
print(maxva)
median_value <- median(obj$nFeature_Spatial)
print(median_value)

obj <- subset(obj, subset = nCount_Spatial >= minCount & nCount_Spatial <= maxCount &
                nFeature_Spatial >= minFeature & nFeature_Spatial <= maxFeature)

#' lasso object based on tissue parameter, need further test
if (!is.null(tissue)){
  obj <- Lasso(obj, tissue, sample)
}

obj$cell <- paste0("sample:", obj$x, '_', obj$y)
obj$x_y=paste0("x",obj$x,"_y",obj$y)

# finaldata<-list(rds=obj,unique=unique)
return(obj)

}


runrds<-function(rds,samplelist){

  ##rds
  rds$x_y=paste0("x",rds$x,"_y",rds$y)
  rds<-subset(x=rds,x_y %in% samplelist$x_y)

return(rds)
}

RUNSCT<- function(obj){
  obj[["SCT"]] <-CreateAssayObject(counts = as.matrix(obj@assays$Spatial$counts))
  DefaultAssay(obj) <- "SCT"
  obj <-NormalizeData(obj,
                      normalization.method = "LogNormalize",
                      scale.factor = 1e4)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  return(obj)
}
