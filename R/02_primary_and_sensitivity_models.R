#!/usr/bin/env Rscript

source(file.path(Sys.getenv("LOWVOLTAGE_PROJECT_ROOT", unset = "."), "R", "00_functions.R"))

nononc_mids <- readRDS(file.path(cache_dir, "nononc_landmark_m40_maxit10.rds"))
inclusive_mids <- readRDS(file.path(cache_dir, "malignancy_inclusive_landmark_with_extension_m40_maxit10.rds"))
nononc_spec <- readRDS(file.path(cache_dir, "nononc_spline_specification.rds"))
inclusive_spec <- readRDS(file.path(cache_dir, "inclusive_spline_specification.rds"))

nononc_formulas <- model_formulas(nononc_spec, include_malignancy = FALSE)
inclusive_formulas <- model_formulas(inclusive_spec, include_malignancy = TRUE)

observed_lactate <- inclusive_mids$data$lab_24h_lactate_first
observed_lactate <- observed_lactate[is.finite(observed_lactate)]
lactate_spline_spec <- list(
  knots = unname(quantile(observed_lactate, probs = c(1 / 3, 2 / 3), type = 7)),
  boundary = range(observed_lactate)
)
inclusive_model3_terms <- attr(terms(inclusive_formulas$model3), "term.labels")
inclusive_log_lactate_formula <- reformulate(
  ifelse(
    inclusive_model3_terms == "lab_24h_lactate_first",
    "log1p(lab_24h_lactate_first)",
    inclusive_model3_terms
  ),
  response = "landmark_day28_outcome"
)
inclusive_lactate_spline_formula <- reformulate(
  ifelse(
    inclusive_model3_terms == "lab_24h_lactate_first",
    ns_term("lab_24h_lactate_first", lactate_spline_spec),
    inclusive_model3_terms
  ),
  response = "landmark_day28_outcome"
)

message_time("Fitting sequential non-oncologic landmark models")
sequential_results <- bind_rows(
  fit_mi_logistic(
    nononc_mids, nononc_formulas$model1, "low_qrs_landmark",
    "Model 1: unadjusted"
  ),
  fit_mi_logistic(
    nononc_mids, nononc_formulas$model2, "low_qrs_landmark",
    "Model 2: demographics and comorbidities"
  ),
  fit_mi_logistic(
    nononc_mids, nononc_formulas$model3, "low_qrs_landmark",
    "Model 3: flexible primary model with categorical ECG count"
  )
)

message_time("Fitting prespecified landmark sensitivity models")
first_ecg_formulas <- model_formulas(nononc_spec, exposure = "first_ecg_low_qrs", include_malignancy = FALSE)
sensitivity_results <- bind_rows(
  fit_mi_logistic(
    nononc_mids, nononc_formulas$model3_continuous_count, "low_qrs_landmark",
    "Continuous ECG-count adjustment"
  ),
  fit_mi_logistic(
    nononc_mids, first_ecg_formulas$model3, "first_ecg_low_qrs",
    "Exposure defined by first post-ICU ECG"
  ),
  fit_mi_logistic(
    nononc_mids, nononc_formulas$model3_no_ecg_count, "low_qrs_landmark",
    "Exactly one ECG in first 24 hours",
    subset_rows = function(data) filter(data, n_ecg_pm24 == 1L)
  ),
  fit_mi_logistic(
    nononc_mids, nononc_formulas$model3, "low_qrs_landmark",
    "First eligible landmark stay per patient",
    subset_rows = function(data) filter(data, first_patient_flag == 1L)
  ),
  fit_mi_logistic(
    nononc_mids, nononc_formulas$model3_no_lactate, "low_qrs_landmark",
    "Flexible Model 3 without lactate"
  ),
  fit_mi_logistic(
    nononc_mids, nononc_formulas$model3, "low_qrs_landmark",
    "Sepsis-3 suspected-infection time within 24 hours of ICU admission",
    subset_rows = function(data) filter(
      data,
      !is.na(sepsis_hours_from_icu),
      abs(sepsis_hours_from_icu) <= 24
    )
  ),
  fit_mi_logistic(
    nononc_mids, nononc_formulas$expanded, "low_qrs_landmark",
    "Expanded acuity model"
  )
)

message_time("Fitting complete-case primary model")
complete_data <- nononc_mids$data %>%
  drop_na(
    landmark_day28_outcome, low_qrs_landmark, age, gender, weight,
    all_of(comorbidities), sofa, lab_24h_lactate_first, n_ecg_cat
  )
complete_fit <- fit_cluster_logistic(complete_data, nononc_formulas$model3)
complete_beta <- coef(complete_fit$model)["low_qrs_landmark"]
complete_variance <- complete_fit$covariance["low_qrs_landmark", "low_qrs_landmark"]
complete_df <- n_distinct(complete_data$subject_id) - complete_fit$model$rank
complete_critical <- qt(0.975, complete_df)
complete_case_result <- tibble(
  analysis = "Complete-case flexible Model 3",
  n = nrow(complete_data),
  exposed = sum(complete_data$low_qrs_landmark),
  deaths = sum(complete_data$landmark_day28_outcome),
  patients = n_distinct(complete_data$subject_id),
  parameters = complete_fit$model$rank,
  m = 0L,
  estimate = exp(complete_beta),
  conf_low = exp(complete_beta - complete_critical * sqrt(complete_variance)),
  conf_high = exp(complete_beta + complete_critical * sqrt(complete_variance)),
  p_value = 2 * pt(abs(complete_beta / sqrt(complete_variance)), complete_df, lower.tail = FALSE),
  df = complete_df,
  within_variance = complete_variance,
  between_variance = NA_real_,
  total_variance = complete_variance,
  relative_increase_variance = NA_real_,
  fraction_missing_information = NA_real_,
  monte_carlo_error_fraction_of_se = NA_real_
)

message_time("Fitting malignancy-inclusive landmark model with the same flexible specification")
inclusive_results <- bind_rows(
  fit_mi_logistic(
    inclusive_mids, inclusive_formulas$model1, "low_qrs_landmark",
    "Malignancy-inclusive Model 1"
  ),
  fit_mi_logistic(
    inclusive_mids, inclusive_formulas$model2, "low_qrs_landmark",
    "Malignancy-inclusive Model 2"
  ),
  fit_mi_logistic(
    inclusive_mids, inclusive_formulas$model3, "low_qrs_landmark",
    "Malignancy-inclusive flexible Model 3"
  )
)

inclusive_first_ecg_formulas <- model_formulas(
  inclusive_spec,
  exposure = "first_ecg_low_qrs",
  include_malignancy = TRUE
)
inclusive_sensitivity_results <- bind_rows(
  fit_mi_logistic(
    inclusive_mids, inclusive_formulas$model3_continuous_count, "low_qrs_landmark",
    "Malignancy-inclusive continuous ECG-count adjustment"
  ),
  fit_mi_logistic(
    inclusive_mids, inclusive_first_ecg_formulas$model3, "first_ecg_low_qrs",
    "Malignancy-inclusive exposure defined by first post-ICU ECG"
  ),
  fit_mi_logistic(
    inclusive_mids, inclusive_formulas$model3_no_ecg_count, "low_qrs_landmark",
    "Malignancy-inclusive exactly one ECG in first 24 hours",
    subset_rows = function(data) filter(data, n_ecg_pm24 == 1L)
  ),
  fit_mi_logistic(
    inclusive_mids, inclusive_formulas$model3, "low_qrs_landmark",
    "Malignancy-inclusive first eligible landmark stay per patient",
    subset_rows = function(data) filter(data, first_patient_flag == 1L)
  ),
  fit_mi_logistic(
    inclusive_mids, inclusive_formulas$model3_no_lactate, "low_qrs_landmark",
    "Malignancy-inclusive flexible Model 3 without lactate"
  ),
  fit_mi_logistic(
    inclusive_mids, inclusive_log_lactate_formula, "low_qrs_landmark",
    "Malignancy-inclusive log-transformed lactate model"
  ),
  fit_mi_logistic(
    inclusive_mids, inclusive_lactate_spline_formula, "low_qrs_landmark",
    "Malignancy-inclusive lactate-spline model"
  ),
  fit_mi_logistic(
    inclusive_mids, inclusive_formulas$model3, "low_qrs_landmark",
    "Malignancy-inclusive Sepsis-3 timing within 24 hours of ICU admission",
    subset_rows = function(data) filter(
      data,
      !is.na(sepsis_hours_from_icu),
      abs(sepsis_hours_from_icu) <= 24
    )
  ),
  fit_mi_logistic(
    inclusive_mids, inclusive_formulas$expanded, "low_qrs_landmark",
    "Malignancy-inclusive expanded acuity model"
  )
)

inclusive_complete_data <- inclusive_mids$data %>%
  drop_na(
    landmark_day28_outcome, low_qrs_landmark, age, gender, weight,
    all_of(comorbidities), icd_malignancy, sofa, lab_24h_lactate_first,
    n_ecg_cat
  )
inclusive_complete_fit <- fit_cluster_logistic(inclusive_complete_data, inclusive_formulas$model3)
inclusive_complete_beta <- coef(inclusive_complete_fit$model)["low_qrs_landmark"]
inclusive_complete_variance <- inclusive_complete_fit$covariance["low_qrs_landmark", "low_qrs_landmark"]
inclusive_complete_df <- n_distinct(inclusive_complete_data$subject_id) - inclusive_complete_fit$model$rank
inclusive_complete_critical <- qt(0.975, inclusive_complete_df)
inclusive_complete_case_result <- tibble(
  analysis = "Malignancy-inclusive complete-case flexible Model 3",
  n = nrow(inclusive_complete_data),
  exposed = sum(inclusive_complete_data$low_qrs_landmark),
  deaths = sum(inclusive_complete_data$landmark_day28_outcome),
  patients = n_distinct(inclusive_complete_data$subject_id),
  parameters = inclusive_complete_fit$model$rank,
  m = 0L,
  estimate = exp(inclusive_complete_beta),
  conf_low = exp(inclusive_complete_beta - inclusive_complete_critical * sqrt(inclusive_complete_variance)),
  conf_high = exp(inclusive_complete_beta + inclusive_complete_critical * sqrt(inclusive_complete_variance)),
  p_value = 2 * pt(
    abs(inclusive_complete_beta / sqrt(inclusive_complete_variance)),
    inclusive_complete_df,
    lower.tail = FALSE
  ),
  df = inclusive_complete_df,
  within_variance = inclusive_complete_variance,
  between_variance = NA_real_,
  total_variance = inclusive_complete_variance,
  relative_increase_variance = NA_real_,
  fraction_missing_information = NA_real_,
  monte_carlo_error_fraction_of_se = NA_real_
)

message_time("Estimating marginal standardized risks")
nononc_marginal <- pool_marginal_standardization(
  nononc_mids,
  nononc_formulas$model3
)
inclusive_marginal <- pool_marginal_standardization(
  inclusive_mids,
  inclusive_formulas$model3
)

message_time("Performing pooled exposure-by-ECG-count interaction test")
interaction_formula <- update(
  nononc_formulas$model3,
  . ~ . + low_qrs_landmark:n_ecg_cat
)
interaction_models <- vector("list", nononc_mids$m)
reduced_models <- vector("list", nononc_mids$m)
for (index in seq_len(nononc_mids$m)) {
  completed <- complete(nononc_mids, index)
  interaction_models[[index]] <- glm(interaction_formula, data = completed, family = binomial())
  reduced_models[[index]] <- glm(nononc_formulas$model3, data = completed, family = binomial())
}
interaction_d1 <- mice::D1(
  as.mira(interaction_models),
  as.mira(reduced_models),
  dfcom = n_distinct(nononc_mids$data$subject_id) - max(vapply(interaction_models, function(x) x$rank, integer(1)))
)
interaction_test <- as.data.frame(interaction_d1$result) %>%
  as_tibble() %>%
  mutate(test = "D1 pooled exposure-by-ECG-count interaction", .before = 1)

inclusive_interaction_formula <- update(
  inclusive_formulas$model3,
  . ~ . + low_qrs_landmark:n_ecg_cat
)
inclusive_interaction_models <- vector("list", inclusive_mids$m)
inclusive_reduced_models <- vector("list", inclusive_mids$m)
for (index in seq_len(inclusive_mids$m)) {
  completed <- complete(inclusive_mids, index)
  inclusive_interaction_models[[index]] <- glm(
    inclusive_interaction_formula,
    data = completed,
    family = binomial()
  )
  inclusive_reduced_models[[index]] <- glm(
    inclusive_formulas$model3,
    data = completed,
    family = binomial()
  )
}
inclusive_interaction_d1 <- mice::D1(
  as.mira(inclusive_interaction_models),
  as.mira(inclusive_reduced_models),
  dfcom = n_distinct(inclusive_mids$data$subject_id) -
    max(vapply(inclusive_interaction_models, function(x) x$rank, integer(1)))
)
inclusive_interaction_test <- as.data.frame(inclusive_interaction_d1$result) %>%
  as_tibble() %>%
  mutate(test = "Malignancy-inclusive D1 pooled exposure-by-ECG-count interaction", .before = 1)

message_time("Performing MI-pooled D2 test of lactate nonlinearity")
inclusive_lactate_spline_models <- vector("list", inclusive_mids$m)
inclusive_linear_lactate_models <- vector("list", inclusive_mids$m)
for (index in seq_len(inclusive_mids$m)) {
  completed <- complete(inclusive_mids, index)
  inclusive_lactate_spline_models[[index]] <- glm(
    inclusive_lactate_spline_formula,
    data = completed,
    family = binomial()
  )
  inclusive_linear_lactate_models[[index]] <- glm(
    inclusive_formulas$model3,
    data = completed,
    family = binomial()
  )
}
inclusive_lactate_d2 <- mice::D2(
  as.mira(inclusive_lactate_spline_models),
  as.mira(inclusive_linear_lactate_models),
  use = "likelihood"
)
inclusive_lactate_nonlinearity_test <- as.data.frame(inclusive_lactate_d2$result) %>%
  as_tibble() %>%
  mutate(test = "MI-pooled D2 likelihood test of lactate nonlinearity", .before = 1)

all_results <- bind_rows(
  sequential_results,
  sensitivity_results,
  complete_case_result,
  inclusive_results,
  inclusive_sensitivity_results,
  inclusive_complete_case_result
) %>%
  mutate(
    estimate_ci = sprintf("%.3f (%.3f-%.3f)", estimate, conf_low, conf_high),
    across(c(estimate, conf_low, conf_high), ~ round(.x, 6)),
    p_value = signif(p_value, 6),
    df = round(df, 2),
    fraction_missing_information = signif(fraction_missing_information, 5),
    monte_carlo_error_fraction_of_se = signif(monte_carlo_error_fraction_of_se, 5)
  )

write_csv(all_results, file.path(results_dir, "final_landmark_model_results.csv"))
write_csv(
  bind_rows(interaction_test, inclusive_interaction_test),
  file.path(results_dir, "pooled_ecg_count_interaction_D1.csv")
)
write_csv(
  inclusive_lactate_nonlinearity_test,
  file.path(results_dir, "pooled_lactate_nonlinearity_D2.csv")
)
write_csv(
  tibble(
    variable = "First lactate, mmol/L",
    internal_knot_1 = lactate_spline_spec$knots[1],
    internal_knot_2 = lactate_spline_spec$knots[2],
    boundary_knot_low = lactate_spline_spec$boundary[1],
    boundary_knot_high = lactate_spline_spec$boundary[2]
  ),
  file.path(diagnostics_dir, "lactate_spline_specification.csv")
)
write_csv(nononc_marginal$by_imputation, file.path(diagnostics_dir, "nononc_marginal_by_imputation.csv"))
write_csv(nononc_marginal$pooled, file.path(results_dir, "nononc_primary_marginal_standardization.csv"))
write_csv(inclusive_marginal$by_imputation, file.path(diagnostics_dir, "inclusive_marginal_by_imputation.csv"))
write_csv(inclusive_marginal$pooled, file.path(results_dir, "malignancy_inclusive_marginal_standardization.csv"))

print(all_results, n = Inf)
print(interaction_test)
print(inclusive_lactate_nonlinearity_test)
print(nononc_marginal$pooled)
print(inclusive_marginal$pooled)
message_time("Primary and sensitivity models completed")
