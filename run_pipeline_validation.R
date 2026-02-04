# Master Validation Script for Analytical Pipeline (Scheme 3)
# ==============================================================================
# This script executes the modular unit tests for the core analytical pipeline.
# It uses the downsampled example data to verify:
# 1. Preprocessing & Integration
# 2. Clustering Logic
# 3. Differential Analysis Logic
# ==============================================================================

pipeline_dir <- "./validation_pipeline"
scripts <- list.files(pipeline_dir, pattern = "\\.R$", full.names = TRUE)

if(length(scripts) == 0) {
  stop("No pipeline validation scripts found!")
}

message("Starting SMIntegration Analytical Pipeline Validation...")
message("=======================================================")

pass_count <- 0
fail_count <- 0

for (script in scripts) {
  message(sprintf(">> Executing Module: %s", basename(script)))
  tryCatch({
    # Use 'list' reporter to show detailed output for every test case
    # This provides visual confirmation of what is actually being tested
    testthat::test_file(script, reporter = "list")
    pass_count <- pass_count + 1
  }, error = function(e) {
    message(sprintf("!! FAILED: %s", basename(script)))
    message(e$message)
    fail_count <- fail_count + 1
  })
  message("-------------------------------------------------------")
}

message(sprintf("Pipeline Validation Complete. Modules Passed: %d/%d", pass_count, length(scripts)))

if (fail_count == 0) {
  message("SUCCESS: All core analytical modules are functioning correctly.")
} else {
  stop("FAILURE: Issues detected in analytical modules.")
}

