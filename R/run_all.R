#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_argument)) sub("^--file=", "", script_argument[[1]]) else "R/run_all.R"
root <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = TRUE)
Sys.setenv(LOWVOLTAGE_PROJECT_ROOT = root)

scripts <- c(
  "01_prepare_and_impute.R",
  "02_primary_and_sensitivity_models.R",
  "03_primary_diagnostics.R",
  "04_mice_diagnostics.R",
  "05_broad_and_survival_models.R",
  "06_tables_and_figures.R"
)
for (script in scripts) {
  message("Running ", script)
  source(file.path(root, "R", script), chdir = FALSE)
}
message("All statistical scripts completed")
