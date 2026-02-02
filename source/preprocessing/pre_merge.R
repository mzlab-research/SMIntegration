#' @title Helper Function to Create SpatialImage Object
#' @description Internal helper to create a VisiumV1 object for Seurat.
#' @param image Matrix representing the tissue image.
#' @param scale.factors List of scale factors.
#' @param tissue.positions Data frame of tissue positions.
#' @param filter.matrix Logical. Whether to filter by tissue coverage.
#' @return A VisiumV1 object.
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

#' @title Create Seurat Object from Matrix and Coordinates
#' @description Wraps sparse matrix and coordinates into a spatial Seurat object.
#' @param mat Sparse matrix (Features x Spots).
#' @param cell_coords Data frame of coordinates (x, y).
#' @param binsize Numeric. Bin size.
#' @return A Seurat object.
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

#' @title Merge Spatial Datasets
#' @description Combines two spatial matrices (e.g., Metabolomics and Transcriptomics) into a single Seurat object
#' by feature concatenation (rbind). Requires identical column names (spots).
#' @param decon_mtrx Matrix 1 (e.g., Metabolomics counts).
#' @param decon_ttrx Matrix 2 (e.g., Transcriptomics counts).
#' @return A merged Seurat object containing features from both modalities.
run_merge_rds<-function(decon_mtrx,decon_ttrx){
  stopifnot(identical(colnames(decon_mtrx), colnames(decon_ttrx)))
  mz_decon_mtrx = decon_mtrx
  st_decon_mtrx = decon_ttrx

  # Combine features
  mat <- rbind(mz_decon_mtrx,st_decon_mtrx)
  
   mat <- rbind2(decon_mtrx, decon_ttrx)
  # Extract coordinates from spot names (assumed format: prefix:x_y or similar)
  cell_coords <- data.frame(spot=colnames(mat)) %>%
    tidyr::separate(col = spot, 
             into = c("prefix", "x", "y"), 
             sep = "[:\\.]|_") %>%
    dplyr::select(-prefix) %>%
    dplyr::mutate(x = as.numeric(x), y = as.numeric(y))
  merge_rds<-create_seurat(mat,cell_coords)
  return(merge_rds)
}

