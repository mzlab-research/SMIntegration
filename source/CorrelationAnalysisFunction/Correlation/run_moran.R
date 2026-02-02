#' @title Compute Pairwise Bivariate Moran's I
#' @description Calculates the bivariate Moran's I correlation matrix for all pairs of variables starting with "Metabolite" or "Gene".
#' Uses a k-nearest neighbor (k=5) spatial weight matrix.
#' Computations are parallelized.
#' @param df Data frame containing 'x', 'y' coordinates and feature columns.
#' @return A list containing the correlation matrix (moran_matrix) and the p-value matrix (p_matrix).
run_moran <- function(df) {
  # Create spatial weight matrix using k-nearest neighbors (k=5)
  coords <- as.matrix(df[, c("x", "y")])
  knn <- spdep::knearneigh(coords, k = 5)
  nb <- spdep::knn2nb(knn)
  weights <- spdep::nb2listw(nb, style = "W")
  
  # Identify all variables starting with "Metabolite" or "Gene"
  all_vars <- grep("^Metabolite|^Gene", names(df), value = TRUE)
  n <- length(all_vars)
  var_idx <- seq_along(all_vars)
  
  # Generate all unique pairs of variables
  pairs <- combn(var_idx, 2, simplify = FALSE)
  
  # Function to calculate Moran's I for a single pair
  calc_moran_pair <- function(pair) {
    i <- pair[1]
    j <- pair[2]
    var1 <- all_vars[i]
    var2 <- all_vars[j]
    
    # Calculate bivariate Moran's I with permutation test
    res <- spdep::moran_bv(
      x = df[[var1]], 
      y = df[[var2]], 
      listw = weights, 
      nsim = 99  # Number of permutations
    )
    
    # Return indices, Moran's I statistic, and p-value
    c(i, j, round(res$t0, 3), round(mean(abs(res$t) >= abs(res$t0)), 3))
  }
  
  # Set up parallel computation
  nCores <- min(4, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(nCores)
  
  # Export necessary objects to cluster
  parallel::clusterExport(cl, c("df", "weights", "all_vars"), envir = environment())
  parallel::clusterEvalQ(cl, library(spdep))
  
  # Run computations in parallel
  results <- parallel::parLapplyLB(cl, pairs, calc_moran_pair)
  parallel::stopCluster(cl)
  
  # Initialize result matrices
  moran_matrix <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(all_vars, all_vars))
  p_matrix <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(all_vars, all_vars))
  
  # Fill matrices with results
  for (res in results) {
    i <- res[1]
    j <- res[2]
    moran_ij <- res[3]
    p_ij <- res[4]
    
    moran_matrix[i, j] <- moran_ij
    moran_matrix[j, i] <- moran_ij
    p_matrix[i, j] <- p_ij
    p_matrix[j, i] <- p_ij
  }
  
  # Set diagonal values
  diag(moran_matrix) <- 1
  diag(p_matrix) <- 0
  
  return(list(moran_matrix, p_matrix))
}

################################
# Main analysis function
#' @title Analyze Specific Spatial Correlations
#' @description Computes bivariate Moran's I between a specific target feature (singlefeature) and a set of other features (otherfeature).
#' Useful for finding top spatially correlated neighbors for a selected molecule.
#' @param df Data frame containing 'x', 'y', and features.
#' @param singlefeature Name of the target feature.
#' @param otherfeature Vector of names of other features to correlate with.
#' @return A list of results, each containing the partner feature name, Moran's I value, and p-value.
analyze_spatial_correlations <- function(df, singlefeature, otherfeature) {
  # Create spatial weight matrix
  coords <- as.matrix(df[, c("x", "y")])
  knn <- spdep::knearneigh(coords, k = 5)
  nb <- spdep::knn2nb(knn)
  weights <- spdep::nb2listw(nb, style = "W")
  
  # Generate all combinations of singlefeature with each otherfeature
  pairs <- expand.grid(k = singlefeature, a = otherfeature, stringsAsFactors = FALSE)
  pairs <- split(pairs, seq(nrow(pairs)))
  
  # Function: Compute bivariate Moran's I and its p-value
  compute_bivariate_moran <- function(pair) {
    pair <- as.character(pair)
    i <- pair[1]
    j <- pair[2]
    
    # Compute bivariate Moran's I
    res <- spdep::moran_bv(
      x = df[[i]], 
      y = df[[j]], 
      listw = weights, 
      nsim = 99  # Reduced number of permutations for speed
    )
    
    # Return results: feature name, Moran's I, and p-value
    c(
      j,
      round(res$t0, 3),
      round(mean(abs(res$t) >= abs(res$t0)), 3)
    )
  }
  
  # Set up optimized parallel computation
  nCores <- min(6, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(nCores)
  
  # Preload necessary packages and data to cluster
  parallel::clusterExport(cl, c("df", "weights"), envir = environment())
  parallel::clusterEvalQ(cl, library(spdep))
  
  # Use parLapplyLB for load-balanced parallel computation
  results <- parallel::parLapplyLB(cl, pairs, compute_bivariate_moran)
  parallel::stopCluster(cl)
  
  return(results)
}