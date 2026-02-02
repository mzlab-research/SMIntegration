#' @title Initial Seurat Object Preprocessing and Filtering
#' @description Filters a Seurat object based on minimum/maximum counts and features.
#' Also adds cell identifiers based on spatial coordinates.
#' @param data Seurat object to be processed.
#' @param minFeature Numeric. Minimum number of features required per spot.
#' @param maxFeature Numeric. Maximum number of features allowed per spot. Defaults to max in dataset.
#' @param sample Character. Sample name prefix.
#' @param minCount Numeric. Minimum total counts required per spot.
#' @param maxCount Numeric. Maximum total counts allowed per spot. Defaults to max in dataset.
#' @param tissue Optional parameter for Lasso filtering (experimental).
#' @return A filtered Seurat object with updated metadata.
run_prerds<-function(data,minFeature=0,
                     maxFeature=NULL,
                     sample="sample",
                     minCount=0,
                     maxCount=NULL,
                     tissue=NULL){

  data$groups<- "0"
  obj<-data

  # Determine max thresholds if not provided
  if (is.null(maxCount)){
    maxCount <- max(obj$nCount_Spatial)
  }
  if (is.null(maxFeature)){
    maxFeature <- max(obj$nFeature_Spatial)
  }

  # Debug output
  maxva<-max(obj$nFeature_Spatial)
  print(maxva)
  median_value <- median(obj$nFeature_Spatial)
  print(median_value)

  # Apply filtering based on count/feature thresholds
  obj <- subset(obj, subset = nCount_Spatial >= minCount & nCount_Spatial <= maxCount &
                  nFeature_Spatial >= minFeature & nFeature_Spatial <= maxFeature)

  # Lasso object based on tissue parameter, need further test
  if (!is.null(tissue)){
    obj <- Lasso(obj, tissue, sample)
  }

  # Create unique cell identifiers
  obj$cell <- paste0("sample:", obj$x, '_', obj$y)
  obj$x_y=paste0("x",obj$x,"_y",obj$y)

  return(obj)
}


#' @title Subset Seurat Object by Coordinates
#' @description Subsets a Seurat object to keep only spots present in a sample list.
#' @param rds Seurat object.
#' @param samplelist Data frame containing a column 'x_y' with coordinates to keep.
#' @return A subsetted Seurat object.
runrds<-function(rds,samplelist){

  ## Create matching key
  rds$x_y=paste0("x",rds$x,"_y",rds$y)
  
  ## Filter object
  rds<-subset(x=rds,x_y %in% samplelist$x_y)

  return(rds)
}

#' @title Flexible Preprocessing Pipeline (SCT-like)
#' @description Performs normalization, transformation, feature selection, and scaling on a Seurat object.
#' Supports specific workflows for metabolomics (TIC/RMS) and standard transcriptomics processing.
#' @param obj Seurat object.
#' @param norm_method Character. 'TIC' (Total Ion Current), 'RMS' (Root Mean Squared), or 'None'.
#' @param transform_method Character. 'LogNormalize' or 'None'.
#' @param scale_data Logical. Whether to perform Z-score scaling.
#' @return Processed Seurat object with a new 'SCT' assay containing normalized/scaled data.
RUNSCT<- function(obj, norm_method="TIC", transform_method="LogNormalize", scale_data=TRUE){
  # Create SCT assay if not present (copy of Spatial)
  if(!"SCT" %in% names(obj@assays)){
      obj[["SCT"]] <- CreateAssayObject(counts = as.matrix(obj@assays$Spatial$counts))
  }
  DefaultAssay(obj) <- "SCT"
  
  # Get raw counts
  counts_mat <- GetAssayData(obj, slot = "counts", assay = "SCT")
  norm_mat <- counts_mat

  # 1. Normalization (Pixel-wise intensity correction)
  if(norm_method == "TIC") {
      # Total Ion Current Normalization (Relative Counts * 1e4)
      col_sums <- colSums(counts_mat)
      col_sums[col_sums == 0] <- 1
      norm_mat <- t(t(counts_mat) / col_sums) * 1e4
  } else if (norm_method == "RMS") {
      # Root Mean Squared Normalization
      # RMS = sqrt(sum(x^2) / n) where n is number of features (rows)
      n_features <- nrow(counts_mat)
      col_sq_sums <- colSums(counts_mat^2)
      rms_vals <- sqrt(col_sq_sums / n_features)
      rms_vals[rms_vals == 0] <- 1
      norm_mat <- t(t(counts_mat) / rms_vals)
  }
  # If "None", norm_mat remains counts_mat

  # Update 'data' slot with normalized values
  obj <- SetAssayData(obj, slot = "data", new.data = as.matrix(norm_mat), assay = "SCT")

  # 2. Transformation (Variance stabilization)
  if(transform_method == "LogNormalize") {
      # If normalization was performed (TIC/RMS), just log transform the normalized values
      if(norm_method %in% c("TIC", "RMS")) {
          curr_data <- GetAssayData(obj, slot = "data", assay = "SCT")
          log_data <- log1p(curr_data)
          obj <- SetAssayData(obj, slot = "data", new.data = log_data, assay = "SCT")
      } else {
          # If No Normalization selected, use standard Seurat LogNormalize
          # This implicitly performs library size normalization followed by log transformation
          obj <- NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 1e4)
      }
  }

  # 3. Feature Selection (Top variable features)
  obj <- FindVariableFeatures(obj)

  # 4. Scaling (Z-score for PCA)
  if(scale_data) {
      obj <- ScaleData(obj)
  }
  
  return(obj)
}
