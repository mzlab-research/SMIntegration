colname_change<-function(data){

colnames(data)[1]<-"geneID"     
colnames(data)[2]<-"x"     
colnames(data)[3]<-"y"     
colnames(data)[4]<-"MIDCount"
data$geneID<-as.character(data$geneID)
data$x<-as.numeric(data$x)
data$y<-as.numeric(data$y)
data$MIDCount<-as.numeric(data$MIDCount)
return(data)
}
###########################################################runxy
runxy<-function(input,binsize=1,minFeature=0,
                maxFeature=NULL,
                sample="sample",
                minCount=0,
                maxCount=NULL){
data <- input
binsize<-as.numeric(binsize)
data %<>% arrange(y) %>% arrange(x)
if(length(grep("group",colnames(data)))==0){
  data$group<-"0"
}
data$x <- trunc(data$x / binsize) * binsize
data$y <- trunc(data$y / binsize) * binsize

samplelist<-data.frame(cell=paste0(sample, ':', data$x, '_', data$y),group=data$group)
  data<-data %>%
    dplyr::select(geneID,x,y,MIDCount)
  data$MIDCount <- as.numeric(data$MIDCount) 
  data <- data[, .(counts=sum(MIDCount)), by = .(geneID, x, y)]

#' create sparse matrix from stereo
data$cell <- paste0(sample, ':', data$x, '_', data$y)
data$geneIdx <- match(data$geneID, unique(data$geneID))
data$cellIdx <- match(data$cell, unique(data$cell))

mat <- sparseMatrix(i = data$geneIdx, j = data$cellIdx, x = data$counts,
                    dimnames = list(unique(data$geneID), unique(data$cell)))


  unique_rows <- unique(data[, c('cell', 'x', 'y')])
  rm(data);gc()

  groups <-samplelist[match(unique_rows$cell, samplelist$cell), "group"]
  cell_coords <- cbind(unique_rows, groups)
rownames(cell_coords) <- cell_coords$cell

seurat_spatialObj <- CreateSeuratObject(counts = mat, project = 'Stereo', assay = 'Spatial',
                                        names.delim = ':', meta.data = cell_coords)

rm(mat);gc()
#' create pseudo image
cell_coords$x <- cell_coords$x - min(cell_coords$x) + 1
cell_coords$y <- cell_coords$y - min(cell_coords$y) + 1

tissue_lowres_image <- matrix(1, max(cell_coords$y), max(cell_coords$x))

tissue_positions_list <- data.frame(row.names = cell_coords$cell,
                                    tissue = 1,
                                    row = cell_coords$y, col = cell_coords$x,
                                    imagerow = cell_coords$y, imagecol = cell_coords$x)


scalefactors_json <- toJSON(list(fiducial_diameter_fullres = binsize,
                                 tissue_hires_scalef = 1,
                                 tissue_lowres_scalef = 1))

#' function to create image object
generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE){
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
  }

  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef

  spot.radius <- unnormalized.radius / max(dim(x = image))

  return(new(Class = 'VisiumV1',
             image = image,
             scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef,
                                          fiducial = scale.factors$fiducial_diameter_fullres,
                                          hires = scale.factors$tissue_hires_scalef,
                                          lowres = scale.factors$tissue_lowres_scalef),
             coordinates = tissue.positions,
             spot.radius = spot.radius))
}

spatialObj <- generate_spatialObj(image = tissue_lowres_image,
                                  scale.factors = fromJSON(scalefactors_json),
                                  tissue.positions = tissue_positions_list)

#' import image into seurat object
spatialObj <- spatialObj[Cells(x = seurat_spatialObj)]
DefaultAssay(spatialObj) <- 'Spatial'

seurat_spatialObj[['slice1']] <- spatialObj

#' filter out empty cell
seurat_spatialObj <- subset(seurat_spatialObj, subset = nCount_Spatial > 0)

############ start seurat processing

obj <- seurat_spatialObj

#' filter object based on nCount and nFeature
if (is.null(maxCount)){
  maxCount <- max(obj$nCount_Spatial)
}
if (is.null(maxFeature)){
  maxFeature <- max(obj$nFeature_Spatial)
}


obj <- subset(obj, subset = nCount_Spatial >= minCount & nCount_Spatial <= maxCount &
                nFeature_Spatial >= minFeature & nFeature_Spatial <= maxFeature)

obj$cell <- paste0("sample:", obj$x, '_', obj$y)
obj$x_y=paste0("x",obj$x,"_y",obj$y)

return(obj)
}





generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE){
  if (filter.matrix) {
    tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
  }
  
  unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef
  spot.radius <- unnormalized.radius / max(dim(x = image))
  
  return(new(Class = 'VisiumV1', 
             image = image, 
             scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, 
                                          fiducial = scale.factors$fiducial_diameter_fullres, 
                                          hires = scale.factors$tissue_hires_scalef, 
                                          lowres = scale.factors$tissue_lowres_scalef), 
             coordinates = tissue.positions, 
             spot.radius = spot.radius))
}

create_seurat <- function(mat,cell_coords,binsize){
  rownames(cell_coords) <- colnames(mat)
  colnames(cell_coords) <- c("x","y")
  seurat_spatialObj <- CreateSeuratObject(counts = mat, project = 'Stereo', assay = 'Spatial',names.delim = ':', meta.data = cell_coords)
  
  cell_coords$x <- cell_coords$x - min(cell_coords$x) + 1
  cell_coords$y <- cell_coords$y - min(cell_coords$y) + 1
  tissue_positions_list <- data.frame(row.names = rownames(cell_coords),
                                      tissue = 1,
                                      row = cell_coords$y, col = cell_coords$x,
                                      imagerow = cell_coords$y, imagecol = cell_coords$x)
  
  tissue_lowres_image <- matrix(1, max(cell_coords$y), max(cell_coords$x))
  scalefactors_json <- toJSON(list(fiducial_diameter_fullres = binsize,
                                   tissue_hires_scalef = binsize,
                                   tissue_lowres_scalef = binsize))
  
  spatialObj <- generate_spatialObj(image = tissue_lowres_image, 
                                    scale.factors = fromJSON(scalefactors_json), 
                                    tissue.positions = tissue_positions_list)
  
  spatialObj <- spatialObj[Cells(x = seurat_spatialObj)]
  DefaultAssay(spatialObj) <- 'Spatial'
  
  seurat_spatialObj[['slice1']] <- spatialObj
  return(seurat_spatialObj)
}