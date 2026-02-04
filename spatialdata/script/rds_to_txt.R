# SMIntegration: Spatial Multi-omics Integration Platform
# ==============================================================================
# 
# Purpose:
#   Converts Seurat RDS objects containing spatial transcriptomics and 
#   metabolomics data into text-based formats (long format) for compatibility 
#   with downstream tools or visualization.
#
# Input Requirements:
#   - 'metab.rds': Seurat object for metabolomics data.
#   - 'trans.rds': Seurat object for transcriptomics data.
#   - Input files should be located in the specified project directory.
#
# Mathematical/Analytical Logic:
#   - Loads Seurat objects.
#   - Extracts count matrices and spatial coordinates.
#   - Transforms data from wide matrix format to long format (gene/mz, x, y, count).
#   - Handles string manipulation to clean coordinate identifiers.
#
# Output:
#   - 'trans_bin_filter.txt': Processed transcriptomics data in tab-delimited format.
#   - 'metab_bin_filter.txt': Processed metabolomics data in tab-delimited format.
#
# ==============================================================================

# Load required R packages
# tibble: Provides flexible data frame structures for data manipulation
# Seurat, SeuratObject, SeuratDisk: For single-cell and spatial transcriptomics data analysis, object handling, and format conversion
# ggplot2, patchwork: For data visualization and multi-panel figure layout
# dplyr: Provides syntax for data manipulation and transformation
# data.table: Efficient handling of large datasets
# Matrix: Sparse matrix operations
# rjson: JSON data parsing and generation
# RColorBrewer: Color palette schemes
# magrittr: Pipe operator for improved code readability
# stringr: String manipulation
# tidyr: Data tidying and reshaping
# readr: Efficient reading/writing of tabular data
library(tibble)
library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)
library(Matrix)
library(rjson)
library(RColorBrewer)
library(magrittr)
library(stringr)
library(tidyr)
library(readr)

# Get project directory path from command line arguments
# commandArgs(T) returns command line arguments, args[1] should be the project directory
args <- commandArgs(T)
projectdir <- args[1]
outdir <- projectdir  # Output directory is the same as project directory

# Construct input file paths and read data
# Read spatial metabolomics data
input_m <- file.path(projectdir, "metab.rds")
data_m <- readRDS(file = input_m)

# Read spatial transcriptomics data
input_t <- file.path(projectdir, "trans.rds")
data_t <- readRDS(file = input_t)

# Convert Seurat objects for spatial transcriptomics and metabolomics to data frames
# Note: The code directly uses ft and fm objects, but based on the above code these should be data_t and data_m
# This may be an inconsistency in the code, but we are following the instruction not to modify the code
t = as.data.frame(ft@assays$Spatial$counts)  # Transcriptomics expression matrix
m = as.data.frame(fm@assays$Spatial$counts)  # Metabolomics expression matrix

# Process transcriptomics data
# Convert transcriptomics data from wide to long format
t_long <- t %>%
  rownames_to_column(var = "geneID") %>%  # Convert row names to geneID column
  pivot_longer(
    cols = -geneID,  # Transform all columns except geneID
    names_to = "x_y",  # Store column names in x_y column
    values_to = "value"  # Store expression values in value column
  )

# Clean prefix from coordinate column names
t_long$x_y <- gsub("^sample:", "", t_long$x_y)

# Split coordinate column into separate x and y columns
t_long <- t_long %>%
  separate(col = x_y, into = c("x", "y"), sep = "_")

# Rename columns for clarity
colnames(t_long) = c("geneID", "x", "y", "MIDCount")

# Ensure correct data types for each column
t_long$geneID <- as.character(t_long$geneID)
t_long$x <- as.numeric(t_long$x)
t_long$y <- as.numeric(t_long$y)
t_long$MIDCount <- as.numeric(t_long$MIDCount)

# Output processed transcriptomics data
write.table(t_long, file = file.path(outdir, "trans_bin_filter.txt"), sep = "\t", row.names = F, quote = F)

# Process metabolomics data (same workflow as transcriptomics data)
m_long <- m %>%
  rownames_to_column(var = "geneID") %>%  # Note: This should be metabolite identifiers, not necessarily genes
  pivot_longer(
    cols = -geneID,
    names_to = "x_y",
    values_to = "value"
  )

m_long$x_y <- gsub("^sample:", "", m_long$x_y)

m_long <- m_long %>%
  separate(col = x_y, into = c("x", "y"), sep = "_")

colnames(m_long) = c("mz", "x", "y", "Intensity")  # Rename geneID to mz, representing mass-to-charge ratio of metabolites

# Ensure correct data types for each column
m_long$mz <- as.character(m_long$mz)
m_long$x <- as.numeric(m_long$x)
m_long$y <- as.numeric(m_long$y)
m_long$Intensity <- as.numeric(m_long$Intensity)

# Output processed metabolomics data
write.table(m_long, file = file.path(outdir, "metab_bin_filter.txt"), sep = "\t", row.names = F, quote = F)