
# -------------------------------------------------------------------------
# Unit Test 02: Clustering Analysis Logic
# -------------------------------------------------------------------------
# Purpose:
#   Verify that the clustering function ('run_clusterplot') correctly assigns
#   clusters to the Seurat objects.
#
# Input:
#   - Processed Seurat objects (generated on-the-fly from example data)
#
# Output:
#   - Seurat objects with 'seurat_clusters' metadata.
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(testthat)
  library(Seurat)
  library(dplyr)
  library(jsonlite)
  library(RColorBrewer)
  library(ggplot2)
})

context("Module 2: Clustering Analysis")

# Source required functions
source("../source/preprocessing/run_prerds.R")
source("../source/preprocessing/pre_merge.R")
source("../source/OverallAnalysisFunction/Clustering/clusterplot.R")

# Helper to prepare data (simplified from Module 1)
prepare_data <- function() {
  m <- readRDS("../example_data/pre_metab1.rds")
  t <- readRDS("../example_data/pre_trans1.rds")
  m <- RUNSCT(m, "None", "LogNormalize", TRUE)
  t <- RUNSCT(t, "None", "LogNormalize", TRUE)
  c <- run_merge_rds(m@assays$Spatial$counts, t@assays$Spatial$counts)
  c <- run_prerds(c)
  c <- RUNSCT(c, "None", "LogNormalize", TRUE)
  # S4 logic: bind data
  c@assays$SCT$data=rbind(m@assays$SCT$data, t@assays$SCT$data)
  c@assays$SCT$scale.data=rbind(m@assays$SCT$scale.data, t@assays$SCT$scale.data)
  list(m=m, t=t, c=c)
}

test_that("Clustering Functionality", {
  
  data_list <- prepare_data()
  data_mrds <- data_list$m
  
  # Test Parameters
  resolution <- 0.5
  clustertype <- 1 # Louvain
  
  # Run Clustering on Metabolomics
  # run_clusterplot returns a list: [[1]]=plot_obj, [[2]]=csv_data, [[3]]=seurat_obj
  # We care about [[3]] for logic validation
  
  result <- run_clusterplot(data_mrds, type="Metabolite", pointSize=1, breakseq=50, resolution=resolution, clustertype=clustertype)
  
  expect_type(result, "list")
  expect_length(result, 3)
  
  obj_clustered <- result[[3]]
  
  expect_s4_class(obj_clustered, "Seurat")
  expect_true("seurat_clusters" %in% colnames(obj_clustered@meta.data), label = "Clusters assigned")
  
  n_clusters <- length(unique(obj_clustered$seurat_clusters))
  expect_gt(n_clusters, 1, label = "Found at least 2 clusters")
  
  message(paste("Clustering found", n_clusters, "clusters using Louvain algorithm."))
  
  # Optional: Test K-means (clustertype=0 or 4) if needed, but Louvain is primary.
})
