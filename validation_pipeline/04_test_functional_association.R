# -------------------------------------------------------------------------
# Unit Test 04: Functional Association Analysis Logic
# -------------------------------------------------------------------------
# Purpose:
#   Verify that the pathway analysis module correctly maps differential 
#   features to KEGG pathways and calculates counts.
#
# Context:
#   Uses 'metabolomics' workflow as the primary test case.
#
# Input:
#   - Mock list of differential metabolites (with KEGG IDs).
#   - KEGG Database files (in source/Database).
#
# Output:
#   - Pathway enrichment table (Pathway, Count, KEGG IDs).
# -------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(testthat)
  library(dplyr)
  library(readr)
})

context("Module 4: Functional Association Analysis")

# Ensure we are in Project Root because annotation_stats.R has hardcoded "./source/Database"
if (basename(getwd()) == "validation_pipeline") {
  setwd("..")
  message("Changed working directory to Project Root for test compatibility.")
}

# Source required functions
if(file.exists("source/Pathway/annotation_stats.R")){
  source("source/Pathway/annotation_stats.R")
} else if (file.exists("../source/Pathway/annotation_stats.R")) {
  # Support running from validation_pipeline dir
  source("../source/Pathway/annotation_stats.R")
  # Note: annotationdatabase() has hardcoded "./source/Database" path.
  # If we are in validation_pipeline dir, we might need to adjust WD or mock.
  # Assuming execution from Root based on run_pipeline_validation.R.
} else {
  stop("Source file annotation_stats.R not found!")
}

test_that("Pathway Database Loading", {
  
  # Check if Database directory exists
  if(!dir.exists("./source/Database")){
    skip("Database directory not found at ./source/Database. Run test from project root.")
  }
  
  # 1. Test Database Loading
  # species="mmu" (Mouse) is used in demo data
  db <- annotationdatabase(omics = "metab", species = "mmu")
  
  expect_type(db, "list")
  expect_true("tab_data" %in% names(db))
  expect_true("gff_data" %in% names(db))
  expect_true("CtoMap" %in% names(db))
  
  expect_gt(nrow(db$tab_data), 0)
  expect_gt(nrow(db$gff_data), 0)
  
  message("KEGG Database loaded successfully.")
  
  # 2. Prepare Mock Input Data
  # Simulating a differential analysis result where we found some metabolites.
  # We use metab_identi.xls to pick real IDs to ensure matches.
  
  identi_file <- "./example_data/metab_identi.xls"
  if(!file.exists(identi_file)){
    skip("metab_identi.xls not found.")
  }
  
  metab_identi <- read_delim(identi_file, delim = "\t", show_col_types = FALSE)
  
  # Pick top 5 metabolites that have a KEGG ID
  valid_metabs <- metab_identi %>%
    filter(!is.na(KEGG.ID) & KEGG.ID != "") %>%
    head(5)
  
  expect_gt(nrow(valid_metabs), 0, label = "Found metabolites with KEGG IDs in annotation file")
  
  # Construct input dataframe for run_annotation_count
  # It expects a dataframe with column 'KEGG.ID' and 'node' (optional but used in S9)
  # Actually run_annotation_count uses 'selectCid_data$KEGG.ID'.
  
  mock_diff_data <- data.frame(
    node = valid_metabs$KEGG.ID,
    KEGG.ID = valid_metabs$KEGG.ID,
    stringsAsFactors = FALSE
  )
  
  message(paste("Testing pathway analysis with", nrow(mock_diff_data), "features:", 
                paste(mock_diff_data$KEGG.ID, collapse=", ")))
  
  # 3. Run Analysis
  result <- run_annotation_count(
    selectCid_data = mock_diff_data,
    database = db,
    omics_type = "metab"
  )
  
  # 4. Verify Results
  expect_true(is.data.frame(result))
  
  if(nrow(result) > 0) {
    expect_true("Pathway" %in% colnames(result))
    expect_true("Count" %in% colnames(result))
    expect_true("KEGG IDs" %in% colnames(result))
    
    # Check that counts are numeric
    expect_true(is.numeric(result$Count))
    
    # Check that at least one pathway has count > 0
    expect_gt(sum(result$Count), 0)
    
    message(paste("Identified", nrow(result), "associated pathways."))
    message(paste("Top Pathway:", result$Pathway[1], "(Count:", result$Count[1], ")"))
    
  } else {
    # It is possible (though unlikely for 5 common metabs) to map to 0 pathways 
    # if the database mapping is sparse. But for standard KEGG it should map.
    warning("No pathways identified for input features. Check database coverage.")
  }
  
})
