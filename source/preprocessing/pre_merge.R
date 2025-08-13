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

create_seurat <- function(mat,cell_coords,binsize=1){
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
run_merge_rds<-function(decon_mtrx,decon_ttrx){
  stopifnot(identical(colnames(decon_mtrx), colnames(decon_ttrx)))
  mz_decon_mtrx = decon_mtrx
  st_decon_mtrx = decon_ttrx

  mat <- rbind(mz_decon_mtrx,st_decon_mtrx)
  
   mat <- rbind2(decon_mtrx, decon_ttrx)
  cell_coords <- data.frame(spot=colnames(mat)) %>%
    tidyr::separate(col = spot, 
             into = c("prefix", "x", "y"), 
             sep = "[:\\.]|_") %>%
    dplyr::select(-prefix) %>%
    dplyr::mutate(x = as.numeric(x), y = as.numeric(y))
  merge_rds<-create_seurat(mat,cell_coords)
  return(merge_rds)
}

