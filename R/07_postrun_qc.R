#!/usr/bin/env Rscript

source(file.path(Sys.getenv("LOWVOLTAGE_PROJECT_ROOT", unset = "."), "R", "00_functions.R"))

required_files <- c(
  file.path(results_dir, "final_landmark_model_results.csv"),
  file.path(results_dir, "malignancy_inclusive_marginal_standardization.csv"),
  file.path(results_dir, "final_flexible_model_diagnostics_summary.csv"),
  file.path(results_dir, "primary_landmark_mice_convergence_summary.csv"),
  file.path(figures_dir, "Figure_2_landmark_Kaplan_Meier.tif"),
  file.path(figures_dir, "Figure_3_adjusted_associations.tif"),
  file.path(figures_dir, "figure1_bundle", "Figure_1_cohort_selection_600dpi.png")
)
missing_files <- required_files[!file.exists(required_files)]
if (length(missing_files)) {
  stop("Required outputs are missing: ", paste(missing_files, collapse = ", "))
}

primary_results <- read_csv(file.path(results_dir, "final_landmark_model_results.csv"), show_col_types = FALSE)
primary_row <- primary_results %>%
  filter(analysis == "Malignancy-inclusive flexible Model 3")
stopifnot(
  nrow(primary_row) == 1L,
  primary_row$n == 14129L,
  primary_row$exposed == 3213L,
  primary_row$deaths == 2346L
)

marginal <- read_csv(
  file.path(results_dir, "malignancy_inclusive_marginal_standardization.csv"),
  show_col_types = FALSE
)
stopifnot(
  nrow(marginal) == 1L,
  all(c("risk_unexposed", "risk_exposed", "risk_difference", "risk_ratio") %in% names(marginal))
)

diagnostics <- read_csv(
  file.path(results_dir, "final_flexible_model_diagnostics_summary.csv"),
  show_col_types = FALSE
)
stopifnot(
  diagnostics$minimum[diagnostics$diagnostic == "Model rank"] == 23,
  diagnostics$maximum[diagnostics$diagnostic == "Model rank"] == 23,
  diagnostics$minimum[diagnostics$diagnostic == "Events per fitted parameter"] == 102,
  diagnostics$maximum[diagnostics$diagnostic == "Events per fitted parameter"] == 102
)

figure1_png <- file.path(figures_dir, "figure1_bundle", "Figure_1_cohort_selection_600dpi.png")
figure1_tiff <- file.path(figures_dir, "Figure_1_cohort_selection.tif")
if (requireNamespace("magick", quietly = TRUE)) {
  figure <- magick::image_read(figure1_png) %>%
    magick::image_background("white", flatten = TRUE) %>%
    magick::image_convert(type = "TrueColor", colorspace = "sRGB", matte = FALSE)
  magick::image_write(
    figure,
    figure1_tiff,
    format = "tiff",
    density = "600x600",
    compression = "lzw"
  )
} else {
  warning("Package 'magick' is unavailable; the 600-dpi PNG remains the Figure 1 raster output.")
}

qc_lines <- c(
  "# Post-run QC report",
  "",
  "- PASS: primary risk set = 14,129 ICU stays.",
  "- PASS: primary exposure-positive group = 3,213 ICU stays.",
  "- PASS: deaths after the 24-hour landmark through day 28 = 2,346.",
  "- PASS: the final flexible model has rank 23 and 102 events per fitted parameter.",
  "- PASS: marginal standardized risks, risk difference, and risk ratio were generated.",
  "- PASS: Figure 1-3 source files were generated.",
  if (file.exists(figure1_tiff)) "- PASS: Figure 1 was exported as a white-background RGB, LZW-compressed TIFF." else "- NOTE: Figure 1 TIFF conversion was unavailable; use the 600-dpi PNG.",
  "- PASS: no patient-level input files are required in the submission ZIP; data/input and local_cache must remain excluded."
)
writeLines(qc_lines, file.path(project_root, "outputs", "postrun_QC_report.md"), useBytes = TRUE)
message_time("Post-run QC completed")
