slim_seurat <- function(seurat_obj, keep_counts = FALSE, keep_scale = FALSE, keep_reductions = FALSE){

  if (!"SCT" %in% names(seurat_obj@assays)){
    return(seurat_obj)
    showNotification("SCT assay not found in Seurat object")
  }
  

  if (!keep_counts && length(seurat_obj@assays$SCT@counts) > 0) {
    seurat_obj@assays$SCT@counts <- new("dgCMatrix")
  }
  

  if (!keep_scale && length(seurat_obj@assays$SCT@scale.data) > 0) {
    seurat_obj@assays$SCT@scale.data <- matrix()
  }

  if (!keep_reductions) {
    for (red in names(seurat_obj@reductions)) {
      seurat_obj@reductions[[red]] <- NULL
    }
  }
  

  seurat_obj@commands <- list()
  

  return(seurat_obj)
}

