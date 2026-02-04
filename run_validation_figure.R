# SMIntegration: Master Validation Script
# ==============================================================================
#
# Purpose:
#   This script serves as the master controller for the SMIntegration validation suite.
#   It automatically discovers and executes all unit tests located in the '/validation_figure'
#   directory to verify the correctness of key analytical modules.
#
# Input Requirements:
#   - Directory: '/validation_figure' containing R scripts named with 'test_' prefix (e.g., 'test_clustering.R').
#   - Libraries: 'testthat' package must be installed.
#
# Logic:
#   1. Scans the './validation_figure' directory for test scripts.
#   2. Iteratively executes each script using 'testthat::test_file()'.
#   3. Captures success/failure status and reports a summary.
#
# Output:
#   - Console output detailing the execution of each test file.
#   - Final summary of passed/failed tests.
#   - Stops execution with error if any test fails (useful for CI/CD pipelines).
#
# ==============================================================================

# Master Validation Script for SMIntegration
# Runs all unit tests in the /validation_figure directory

validation_dir <- "./validation_figure"
scripts <- list.files(validation_dir, pattern = "^test_.*\\.R$", full.names = TRUE)

if(length(scripts) == 0) {
  stop("No validation scripts found!")
}

message("Starting SMIntegration Validation Suite...")
message("-----------------------------------------")

pass_count <- 0
fail_count <- 0

for (script in scripts) {
  message(sprintf("Running: %s", basename(script)))
  tryCatch({
    testthat::test_file(script)
    pass_count <- pass_count + 1
  }, error = function(e) {
    message(sprintf("FAILED: %s", basename(script)))
    message(e$message)
    fail_count <- fail_count + 1
  })
  message("-----------------------------------------")
}

message(sprintf("Validation Complete. Passed: %d, Failed: %d", pass_count, fail_count))

if (fail_count == 0) {
  message("All systems operational.")
} else {
  stop("Validation failures detected.")
}

