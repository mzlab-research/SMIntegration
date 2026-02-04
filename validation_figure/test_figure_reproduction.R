context("Manuscript Figure Reproduction")

# ==============================================================================
# Helper Function: Locate Benchmark Directory
# ==============================================================================
# This function attempts to locate the 'benchmarks' directory containing RDS files
# and scripts. It checks common locations relative to the execution context.
find_benchmark_dir <- function() {
  candidates <- c(
    file.path("validation_figure", "benchmarks"), # If running from root
    file.path("benchmarks")                       # If running from validation_figure dir
  )
  for (d in candidates) {
    if (dir.exists(d)) return(d)
  }
  return(NULL)
}

test_that("Manuscript Figure Reproduction (Figs 1-4)", {
  
  # 1. Locate Benchmarks
  benchmark_dir <- find_benchmark_dir()
  skip_if(is.null(benchmark_dir), "Benchmark directory not found. Cannot proceed.")
  
  # Check for figure1.R as a proxy for valid benchmark content
  script_path <- file.path(benchmark_dir, "figure1.R")
  skip_if_not(file.exists(script_path), "figure1.R not found in benchmark directory.")

  # 2. Setup Output Directory
  # Define where the generated figures will be saved for manual inspection.
  # We use the absolute path to ensure it persists after working directory changes.
  project_root <- getwd()
  output_dir <- file.path(project_root, "validation_figure", "test_output")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # Optimization: Convert to short path (8.3 format) on Windows to avoid encoding issues
  # with Chinese characters in the project path (e.g., "分析调研").
  try({
    output_dir <- utils::shortPathName(output_dir)
    message(paste("Debug: Using short path for output:", output_dir))
  }, silent = TRUE)

  # 3. Setup Temporary Execution Environment
  # We use a temporary directory to run the scripts to ensure a clean state
  # and prevent pollution of the source directory.
  temp_dir <- tempdir()
  work_dir <- file.path(temp_dir, "SMIntegration_Validation_Figures")
  
  # Clean up previous runs if any
  if (dir.exists(work_dir)) unlink(work_dir, recursive = TRUE)
  dir.create(work_dir)
  
  # Create subdirectory structure required by the figure scripts
  # The scripts write to "./Figures", so we must create it.
  dir.create(file.path(work_dir, "Figures"))
  
  # 4. Copy Assets to Temporary Directory
  
  # A. Copy RDS Data Files
  # The scripts expect data in "./Figures/rds/".
  # We copy the entire 'rds' folder from benchmarks to work_dir/Figures.
  rds_source <- file.path(benchmark_dir, "rds")
  if (dir.exists(rds_source)) {
    file.copy(rds_source, file.path(work_dir, "Figures"), recursive = TRUE)
  } else {
    warning("RDS directory not found in benchmarks!")
  }
  
  # B. Copy Source Code (Library Functions)
  # Figure scripts (especially Fig 4) depend on 'source/Pathway/...' files.
  if (dir.exists("source")) {
    file.copy("source", work_dir, recursive = TRUE)
  } else if (dir.exists("../source")) {
    file.copy("../source", work_dir, recursive = TRUE)
  }

  # C. Copy Example Data
  # Figure 4 depends on files in 'example_data/'.
  if (dir.exists("example_data")) {
    file.copy("example_data", work_dir, recursive = TRUE)
  } else if (dir.exists("../example_data")) {
    file.copy("../example_data", work_dir, recursive = TRUE)
  }
  
  # 5. Execute Figure Generation Scripts
  # We loop through figure1.R to figure4.R and execute them sequentially.
  scripts <- c("figure1.R", "figure2.R", "figure3.R", "figure4.R")
  
  for (script_name in scripts) {
      source_script_path <- file.path(benchmark_dir, script_name)
      
      # Verify script exists
      if (!file.exists(source_script_path)) {
          warning(paste("Script not found:", script_name))
          next
      }
      
      # Copy script to working directory
      file.copy(source_script_path, file.path(work_dir, script_name), overwrite = TRUE)
      
      message(paste("Running:", script_name, "..."))
      
      # Define log file for capturing stdout/stderr
      log_file <- file.path(work_dir, paste0(script_name, ".log"))
      
      # Run the script using Rscript in the temporary directory
      # We switch working directory temporarily to work_dir
      old_wd_loop <- setwd(work_dir)
      result <- system2("Rscript", args = c(script_name), stdout = log_file, stderr = log_file)
      setwd(old_wd_loop) # Restore WD
      
      # Check execution result
      if (result == 0) {
        message(paste("SUCCESS:", script_name))
      } else {
        # If failure, read and print the last 20 lines of the log for debugging
        log_content <- readLines(log_file, warn = FALSE)
        warning(paste("FAILURE:", script_name, "\nLog (last 20 lines):\n", 
                      paste(tail(log_content, 20), collapse = "\n")))
      }
  }

  # 6. Verify and Collect Output
  # Check if PNG or PDF files were generated in the 'Figures' folder
  generated_figures <- list.files(file.path(work_dir, "Figures"), pattern = "\\.(png|pdf)$", full.names = TRUE)
  
  expect_true(length(generated_figures) > 0, info = "No figures (png/pdf) were generated.")
  
  # Copy the generated figures from the temp dir to the user-visible output directory
  file.copy(generated_figures, output_dir, overwrite = TRUE)
  
  message(paste("Generated", length(generated_figures), "figures. Copies saved to:", output_dir))
  
  # 7. Cleanup
  # Remove the temporary directory to save space
  unlink(work_dir, recursive = TRUE)
})