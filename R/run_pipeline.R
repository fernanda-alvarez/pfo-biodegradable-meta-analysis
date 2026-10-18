# ==============================================================================
# Master Replication Pipeline
# Efficacy and Safety of Biodegradable versus Metallic PFO Occluders
# Systematic Review and Meta-Analysis (v5)
# ==============================================================================

cat("\n========================================================================\n")
cat(" Starting Meta-Analysis Replication Pipeline\n")
cat("========================================================================\n\n")

# 1. Check required R packages
required_pkgs <- c("metafor", "readr", "dplyr", "jsonlite")
missing_pkgs <- required_pkgs[!vapply(required_pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]

if (length(missing_pkgs) > 0) {
  cat("Installing missing packages:", paste(missing_pkgs, collapse = ", "), "\n")
  install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(metafor)
  library(readr)
  library(dplyr)
  library(jsonlite)
})

# Establish project root
PROJ <- normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(PROJ)

run_stage <- function(name, script_path) {
  cat(sprintf("[STAGE] %-30s ... ", name))
  t0 <- proc.time()
  tryCatch({
    source(script_path, local = new.env(parent = globalenv()))
    elapsed <- round((proc.time() - t0)[3], 2)
    cat(sprintf("OK (%s s)\n", elapsed))
    return(TRUE)
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", conditionMessage(e)))
    return(FALSE)
  })
}

# 2. Execution order
stages <- list(
  list(name = "01: Risk of Bias Ledger", script = "R/05b_rob_join.R"),
  list(name = "02: Primary & Core Models", script = "R/02_run_models.R"),
  list(name = "03: Forest Plots & Report", script = "R/03_render_figures_report.R"),
  list(name = "04: Sensitivity Analyses", script = "R/04_sensitivity_analyses.R"),
  list(name = "05: Adjunct Analyses", script = "R/06_adjunct_analyses.R"),
  list(name = "06: Bayesian Random-Effects", script = "R/07_bayesian.R"),
  list(name = "07: Clinical & Person-Time", script = "R/08_extraction_analyses.R"),
  list(name = "08: Network Meta-Analysis", script = "R/09_nma.R"),
  list(name = "09: GRADE & TSA", script = "R/10_grade_tsa.R")
)

stage_results <- logical(length(stages))
for (i in seq_along(stages)) {
  stage_results[i] <- run_stage(stages[[i]]$name, stages[[i]]$script)
}

# 3. Automated Unit Tests
cat("\n========================================================================\n")
cat(" Running Verification Test Suite\n")
cat("========================================================================\n\n")

test_scripts <- c(
  "tests/test_manifest.R",
  "tests/test_models.R",
  "tests/test_sensitivity.R",
  "tests/test_adjunct.R",
  "tests/test_bayesian.R",
  "tests/test_nma.R",
  "tests/test_grade_tsa.R",
  "tests/test_extraction.R"
)

test_results <- logical(length(test_scripts))
for (i in seq_along(test_scripts)) {
  ts <- test_scripts[i]
  cat(sprintf("[TEST]  %-30s ... ", basename(ts)))
  tryCatch({
    source(ts, local = new.env(parent = globalenv()))
    test_results[i] <- TRUE
  }, error = function(e) {
    cat(sprintf("FAILED: %s\n", conditionMessage(e)))
    test_results[i] <- FALSE
  })
}

# 4. Summary
cat("\n========================================================================\n")
cat(sprintf(" Replication Summary: %d/%d stages OK | %d/%d tests PASSED\n",
            sum(stage_results), length(stage_results),
            sum(test_results), length(test_results)))
cat(" Outputs written to 'outputs/tables/' and 'outputs/figures/'.\n")
cat("========================================================================\n\n")

if (!all(stage_results) || !all(test_results)) {
  stop("Replication failed: one or more stages or tests failed.")
}

cat("SUCCESS: All replication stages and tests completed successfully!\n")
