context("Manuscript Benchmark Data Integrity")

# ==============================================================================
# Unit Test: Benchmark Data Integrity (Meta-Test)
# ==============================================================================
# Purpose:
#   This script performs a 'meta-test' on the benchmark dataset itself.
#   It verifies that the RDS files (archived from the original manuscript analysis)
#   are present, readable, and contain data structures consistent with expectations
#   (e.g., Seurat objects with specific metadata columns).
#
#   This ensures that any failure in figure reproduction is due to code/environment
#   issues, not corrupted benchmark data.
# ==============================================================================

# 1. Locate Benchmark Directory
# Try to find the 'benchmarks/rds' directory relative to the execution context.
benchmark_dir <- file.path("benchmarks", "rds")

if (!dir.exists(benchmark_dir)) {
  # Fallback logic for different running contexts (root vs subdir)
  if (dir.exists("benchmarks/rds")) {
    benchmark_dir <- "benchmarks/rds"
  } else {
    benchmark_dir <- file.path("validation_figure", "benchmarks", "rds")
  }
}

test_that("Benchmark directory exists", {
  expect_true(dir.exists(benchmark_dir), 
              info = paste("Benchmark directory not found at:", benchmark_dir))
})

# 2. Validation of Figure 1 Data (cluster.rds)
test_that("cluster.rds contains valid Seurat objects for Figure 1", {
  skip_if_not(file.exists(file.path(benchmark_dir, "cluster.rds")))
  
  data <- readRDS(file.path(benchmark_dir, "cluster.rds"))
  
  # Check structure: List of 3 elements (Metab, Trans, Merge)
  expect_type(data, "list")
  expect_length(data, 3)
  
  m <- data[[1]]
  t <- data[[2]]
  c <- data[[3]]
  
  # Helper function to validate Seurat-like objects
  # Checks for essential metadata columns used in plotting: 'seurat_clusters', 'x', 'y'
  check_obj <- function(obj, name) {
    if (inherits(obj, "Seurat")) {
      expect_true("seurat_clusters" %in% colnames(obj@meta.data), info = paste(name, "missing seurat_clusters"))
      expect_true("x" %in% colnames(obj@meta.data), info = paste(name, "missing x coord"))
      expect_true("y" %in% colnames(obj@meta.data), info = paste(name, "missing y coord"))
    } else {
      # Assume dataframe structure if not Seurat
      expect_true("seurat_clusters" %in% colnames(obj), info = paste(name, "missing seurat_clusters"))
      expect_true("x" %in% colnames(obj), info = paste(name, "missing x coord"))
      expect_true("y" %in% colnames(obj), info = paste(name, "missing y coord"))
    }
  }
  
  check_obj(m, "Metabolite Object")
  check_obj(t, "Transcript Object")
  check_obj(c, "Merge Object")
})

# 3. Validation of Figure 1b Data (find_spatial_pattern.rds)
test_that("find_spatial_pattern.rds contains pattern data for Figure 1b", {
  skip_if_not(file.exists(file.path(benchmark_dir, "find_spatial_pattern.rds")))
  
  pattern_data <- readRDS(file.path(benchmark_dir, "find_spatial_pattern.rds"))
  
  expect_type(pattern_data, "list")
  expect_length(pattern_data, 2) # m and t patterns
  
  m_pattern <- pattern_data[[1]]
  
  # Validate internal structure: Pattern matrix and Location data
  expect_true(is.list(m_pattern))
  expect_true(!is.null(m_pattern[[1]])) # Pattern matrix
  expect_true(!is.null(m_pattern[[2]])) # Location
})