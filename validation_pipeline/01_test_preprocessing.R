
# -------------------------------------------------------------------------
# Unit Test 01: Preprocessing and Data Integration Logic
# -------------------------------------------------------------------------
# Purpose:
#   Verify that raw Seurat objects (Metabolomics & Transcriptomics) can be
#   correctly normalized, scaled, and merged into a combined object.
#   This mimics the logic in 'server/S2_upload.R' and 'server/S4_clustering.R'.
#
# Input:
#   - example_data/pre_metab1.rds
#   - example_data/pre_trans1.rds
#
# Output:
#   - Validated Seurat objects (Metab, Trans, Merge) ready for clustering.
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(testthat)
  library(Seurat)
  library(dplyr)
  library(jsonlite)
})

context("Module 1: Preprocessing & Integration")

# 1. Load Source Functions
# We need RUNSCT (likely in preprocessing/run_prerds.R or similar, but checking S2 code)
# S2 uses RUNSCT. Let's find where it is defined.
# Based on file structure, likely source/preprocessing/run_prerds.R
if(file.exists("../source/preprocessing/run_prerds.R")) {
  source("../source/preprocessing/run_prerds.R")
} else {
  stop("Critical dependency missing: ../source/preprocessing/run_prerds.R")
}

if(file.exists("../source/preprocessing/pre_merge.R")) {
  source("../source/preprocessing/pre_merge.R")
} else {
  stop("Critical dependency missing: ../source/preprocessing/pre_merge.R")
}

test_that("Data Loading and Normalization", {
  
  # 2. Load Demo Data
  # These paths are relative to the 'validation_pipeline' directory
  metab_path <- "../example_data/pre_metab1.rds"
  trans_path <- "../example_data/pre_trans1.rds"
  
  expect_true(file.exists(metab_path), label = "Example Metabolomics data exists")
  expect_true(file.exists(trans_path), label = "Example Transcriptomics data exists")
  
  data_mrds <- readRDS(metab_path)
  data_trds <- readRDS(trans_path)
  
  expect_s4_class(data_mrds, "Seurat")
  expect_s4_class(data_trds, "Seurat")
  
  # 3. Test Metabolomics Preprocessing (RUNSCT)
  # Mimic S4_clustering.R logic: 
  # norm_method="None", transform_method="LogNormalize", scale_data=TRUE
  # Note: S2 allows user choice, we test default "LogNormalize" path.
  
  data_mrds <- RUNSCT(data_mrds, norm_method = "None", transform_method = "LogNormalize", scale_data = TRUE)
  
  expect_true("SCT" %in% names(data_mrds@assays), label = "SCT assay created for Metabolomics")
  expect_true(!is.null(data_mrds@assays$SCT@scale.data), label = "Scale data exists for Metabolomics")
  
  # 4. Test Transcriptomics Preprocessing
  data_trds <- RUNSCT(data_trds, norm_method = "None", transform_method = "LogNormalize", scale_data = TRUE)
  
  expect_true("SCT" %in% names(data_trds@assays), label = "SCT assay created for Transcriptomics")
  
  # 5. Test Integration (Merge)
  decon_mtrx = data_mrds@assays$Spatial$counts
  decon_ttrx = data_trds@assays$Spatial$counts
  
  # run_merge_rds logic
  data_crds <- run_merge_rds(decon_mtrx, decon_ttrx)
  
  expect_s4_class(data_crds, "Seurat")
  # Check if dimensions make sense (cols should match intersection or union depending on logic)
  # Usually strict intersection of spots
  expect_true(ncol(data_crds) <= ncol(data_mrds)) 
  
  # 6. Process Combined Object
  data_crds <- run_prerds(data = data_crds) # Basic setup
  data_crds <- RUNSCT(data_crds, norm_method = "None", transform_method = "LogNormalize", scale_data = TRUE)
  
  # 7. Mock Integration Step (rbind data/scale.data from parts)
  # This is specific logic found in S4_clustering.R
  # For unit test, we just verify the objects are structurally sound
  
  expect_true("SCT" %in% names(data_crds@assays), label = "SCT assay created for Merge")
  
  # Save intermediate objects for next test step?
  # For true unit tests, usually we don't save state between files, 
  # but here we are building a pipeline validation.
  # We will NOT save, but allow next script to re-run or save to a temp RDS if needed.
  # To keep tests independent, next script will re-run this fast logic.
  
  message("Preprocessing module validated successfully.")
})
