run_moran <- function(df) {

  coords <- as.matrix(df[, c("x", "y")])
  knn <- spdep::knearneigh(coords, k = 5)
  nb <- spdep::knn2nb(knn)
  weights <- spdep::nb2listw(nb, style = "W")
  

  all_vars <- grep("^Metabolite|^Gene", names(df), value = TRUE)
  n <- length(all_vars)
  var_idx <- seq_along(all_vars)
  

  pairs <- combn(var_idx, 2, simplify = FALSE)
  

  calc_moran_pair <- function(pair) {
    i <- pair[1]
    j <- pair[2]
    var1 <- all_vars[i]
    var2 <- all_vars[j]
    

    res <- spdep::moran_bv(
      x = df[[var1]], 
      y = df[[var2]], 
      listw = weights, 
      nsim = 99  
    )
    
    c(i, j, round(res$t0,3), round(mean(abs(res$t) >= abs(res$t0)),3))
  }
  

  nCores <- min(4, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(nCores)
  

  parallel::clusterExport(cl, c("df", "weights", "all_vars"), envir = environment())
  parallel::clusterEvalQ(cl, library(spdep))
  

  results <- parallel::parLapplyLB(cl, pairs, calc_moran_pair)
  parallel::stopCluster(cl)
  

  moran_matrix <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(all_vars, all_vars))
  p_matrix <- matrix(NA_real_, nrow = n, ncol = n, dimnames = list(all_vars, all_vars))
  
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
  

  diag(moran_matrix) <- 1
  diag(p_matrix) <- 0
  
  return(list(moran_matrix, p_matrix))
}
################################
# 主分析函数
analyze_spatial_correlations <- function(df, singlefeature, otherfeature) {
  coords <- as.matrix(df[, c("x", "y")])
  knn <- spdep::knearneigh(coords, k = 5)
  nb <- spdep::knn2nb(knn)
  weights <- spdep::nb2listw(nb, style = "W")
  #所有组合
  pairs <- expand.grid(k = singlefeature, a = otherfeature, stringsAsFactors = FALSE)
  pairs <- split(pairs, seq(nrow(pairs)))
  # 函数：计算双变量莫兰指数及其 p 值
  compute_bivariate_moran <- function(pair) {
    pair<-as.character(pair)
    i <- pair[1]
    j <- pair[2]
    # 计算双变量莫兰指数
    res <- spdep::moran_bv(
      x = df[[i]], 
      y = df[[j]], 
      listw = weights, 
      nsim = 99  # 减少置换次数
    )
    
    # 返回结果
    c(
      j,
      round(res$t0,3),
      round(mean(abs(res$t) >= abs(res$t0)),3)
    )
  }
  # 5. 优化的并行计算
  nCores <- min(6, parallel::detectCores() - 1)
  cl <- parallel::makeCluster(nCores)
  
  # 预加载必要包和数据
  parallel::clusterExport(cl, c("df", "weights"), envir = environment())
  parallel::clusterEvalQ(cl, library(spdep))
  
  # 使用parLapplyLB实现负载均衡
  results <- parallel::parLapplyLB(cl, pairs, compute_bivariate_moran)
  parallel::stopCluster(cl)
  return(results)
}
