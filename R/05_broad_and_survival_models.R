#!/usr/bin/env Rscript

source(file.path(Sys.getenv("LOWVOLTAGE_PROJECT_ROOT", unset = "."), "R", "00_functions.R"))

primary_mids <- readRDS(file.path(cache_dir, "malignancy_inclusive_landmark_with_extension_m40_maxit10.rds"))
primary_spec <- readRDS(file.path(cache_dir, "inclusive_spline_specification.rds"))
primary_formulas <- model_formulas(primary_spec, include_malignancy = TRUE)

nononc_mids <- readRDS(file.path(cache_dir, "nononc_landmark_m40_maxit10.rds"))
nononc_spec <- readRDS(file.path(cache_dir, "nononc_spline_specification.rds"))
nononc_formulas <- model_formulas(nononc_spec, include_malignancy = FALSE)

broad_mids <- readRDS(file.path(cache_dir, "malignancy_inclusive_broad_m40_maxit10.rds"))
broad_spec <- readRDS(file.path(cache_dir, "inclusive_broad_spline_specification.rds"))
broad_formulas <- model_formulas(broad_spec, include_malignancy = TRUE)

message_time("Fitting malignancy-inclusive broad-window logistic model")
broad_logistic <- fit_mi_logistic(
  broad_mids,
  broad_formulas$model3,
  "low_qrs_landmark",
  "Malignancy-inclusive broad-window flexible Model 3"
)
broad_marginal <- pool_marginal_standardization(broad_mids, broad_formulas$model3)

message_time("Fitting secondary Cox models")
primary_cox <- fit_mi_cox(
  primary_mids,
  primary_formulas$model3,
  "low_qrs_landmark",
  "Malignancy-inclusive landmark flexible Model 3 Cox"
)
nononc_cox <- fit_mi_cox(
  nononc_mids,
  nononc_formulas$model3,
  "low_qrs_landmark",
  "Non-oncologic landmark flexible Model 3 Cox"
)
broad_cox <- fit_mi_cox(
  broad_mids,
  broad_formulas$model3,
  "low_qrs_landmark",
  "Malignancy-inclusive broad-window flexible Model 3 Cox"
)

cox_results <- bind_rows(
  primary_cox$effect,
  nononc_cox$effect,
  broad_cox$effect
) %>%
  mutate(
    estimate_ci = sprintf("%.3f (%.3f-%.3f)", estimate, conf_low, conf_high),
    across(c(estimate, conf_low, conf_high), ~ round(.x, 6)),
    p_value = signif(p_value, 6),
    df = round(df, 2),
    fraction_missing_information = signif(fraction_missing_information, 5)
  )

summarize_ph <- function(data, label) {
  data %>%
    filter(term %in% c("low_qrs_landmark", "GLOBAL")) %>%
    group_by(term) %>%
    summarise(
      analysis = label,
      median_chisq = median(chisq, na.rm = TRUE),
      median_df = median(df, na.rm = TRUE),
      median_p = median(p, na.rm = TRUE),
      p_q1 = quantile(p, 0.25, na.rm = TRUE),
      p_q3 = quantile(p, 0.75, na.rm = TRUE),
      imputations_p_below_0_05 = sum(p < 0.05, na.rm = TRUE),
      imputations = n(),
      .groups = "drop"
    ) %>%
    select(analysis, everything())
}

ph_summary <- bind_rows(
  summarize_ph(primary_cox$ph, "Malignancy-inclusive landmark Cox"),
  summarize_ph(nononc_cox$ph, "Non-oncologic landmark Cox"),
  summarize_ph(broad_cox$ph, "Malignancy-inclusive broad-window Cox")
)

write_csv(broad_logistic, file.path(results_dir, "malignancy_inclusive_broad_logistic.csv"))
write_csv(broad_marginal$pooled, file.path(results_dir, "malignancy_inclusive_broad_marginal_standardization.csv"))
write_csv(cox_results, file.path(results_dir, "secondary_cox_model_results.csv"))
write_csv(ph_summary, file.path(results_dir, "cox_proportional_hazards_summary.csv"))
write_csv(primary_cox$ph, file.path(diagnostics_dir, "primary_landmark_cox_ph_by_imputation.csv"))
write_csv(broad_cox$ph, file.path(diagnostics_dir, "broad_window_cox_ph_by_imputation.csv"))

print(broad_logistic)
print(broad_marginal$pooled)
print(cox_results)
print(ph_summary)
message_time("Broad-window and survival models completed")
