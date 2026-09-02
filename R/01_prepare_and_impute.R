#!/usr/bin/env Rscript

source(file.path(Sys.getenv("LOWVOLTAGE_PROJECT_ROOT", unset = "."), "R", "00_functions.R"))

nononc <- prepare_landmark_data(
  file.path(input_dir, "nononc_post_icu_24h.csv"),
  file.path(input_dir, "landmark_extension.csv"),
  include_malignancy = FALSE
)
inclusive <- prepare_landmark_data(
  file.path(input_dir, "malignancy_inclusive_post_icu_24h.csv"),
  file.path(input_dir, "landmark_extension_inclusive.csv"),
  include_malignancy = TRUE
)
inclusive_broad <- prepare_broad_data(
  file.path(input_dir, "malignancy_inclusive_broad.csv"),
  include_malignancy = TRUE
)

stopifnot(
  nrow(nononc) == 12285L,
  n_distinct(nononc$subject_id) == 11152L,
  sum(nononc$low_qrs_landmark) == 2695L,
  sum(nononc$landmark_day28_outcome) == 1762L,
  nrow(inclusive) == 14129L,
  n_distinct(inclusive$subject_id) == 12812L,
  sum(inclusive$low_qrs_landmark) == 3213L,
  sum(inclusive$landmark_day28_outcome) == 2346L,
  nrow(inclusive_broad) == 22569L,
  n_distinct(inclusive_broad$subject_id) == 18930L,
  sum(inclusive_broad$low_qrs_landmark) == 4960L,
  sum(inclusive_broad$landmark_day28_outcome) == 4315L
)

nononc_spec <- fixed_spline_specification(nononc)
inclusive_spec <- fixed_spline_specification(inclusive)
saveRDS(nononc_spec, file.path(cache_dir, "nononc_spline_specification.rds"))
saveRDS(inclusive_spec, file.path(cache_dir, "inclusive_spline_specification.rds"))

specification_table <- bind_rows(lapply(names(nononc_spec), function(variable) {
  tibble(
    cohort = "nononcologic",
    variable = variable,
    knot_1 = nononc_spec[[variable]]$knots[1],
    knot_2 = nononc_spec[[variable]]$knots[2],
    boundary_low = nononc_spec[[variable]]$boundary[1],
    boundary_high = nononc_spec[[variable]]$boundary[2]
  )
}), lapply(names(inclusive_spec), function(variable) {
  tibble(
    cohort = "malignancy_inclusive",
    variable = variable,
    knot_1 = inclusive_spec[[variable]]$knots[1],
    knot_2 = inclusive_spec[[variable]]$knots[2],
    boundary_low = inclusive_spec[[variable]]$boundary[1],
    boundary_high = inclusive_spec[[variable]]$boundary[2]
  )
}))
write_csv(specification_table, file.path(diagnostics_dir, "fixed_spline_specification.csv"))

nononc_mids <- make_landmark_mids(
  nononc,
  file.path(cache_dir, "nononc_landmark_m40_maxit10.rds"),
  seed = 20260804L,
  include_malignancy = FALSE
)
inclusive_mids <- make_landmark_mids(
  inclusive,
  file.path(cache_dir, "malignancy_inclusive_landmark_with_extension_m40_maxit10.rds"),
  seed = 20260805L,
  include_malignancy = TRUE
)
inclusive_broad_spec <- fixed_spline_specification(inclusive_broad)
saveRDS(inclusive_broad_spec, file.path(cache_dir, "inclusive_broad_spline_specification.rds"))
inclusive_broad_mids <- make_landmark_mids(
  inclusive_broad,
  file.path(cache_dir, "malignancy_inclusive_broad_m40_maxit10.rds"),
  seed = 20260806L,
  include_malignancy = TRUE
)

write_csv(
  tibble(
    cohort = c("nononcologic_landmark", "malignancy_inclusive_landmark", "malignancy_inclusive_broad"),
    stays = c(nrow(nononc), nrow(inclusive), nrow(inclusive_broad)),
    patients = c(n_distinct(nononc$subject_id), n_distinct(inclusive$subject_id), n_distinct(inclusive_broad$subject_id)),
    exposed = c(sum(nononc$low_qrs_landmark), sum(inclusive$low_qrs_landmark), sum(inclusive_broad$low_qrs_landmark)),
    deaths = c(sum(nononc$landmark_day28_outcome), sum(inclusive$landmark_day28_outcome), sum(inclusive_broad$landmark_day28_outcome)),
    imputations = c(nononc_mids$m, inclusive_mids$m, inclusive_broad_mids$m),
    iterations = c(nononc_mids$iteration, inclusive_mids$iteration, inclusive_broad_mids$iteration),
    logged_events = c(
      if (is.null(nononc_mids$loggedEvents)) 0L else nrow(nononc_mids$loggedEvents),
      if (is.null(inclusive_mids$loggedEvents)) 0L else nrow(inclusive_mids$loggedEvents),
      if (is.null(inclusive_broad_mids$loggedEvents)) 0L else nrow(inclusive_broad_mids$loggedEvents)
    )
  ),
  file.path(results_dir, "landmark_cohort_and_imputation_summary.csv")
)

write_session_information(file.path(results_dir, "R_sessionInfo.txt"))
message_time("Preparation and imputation completed")
