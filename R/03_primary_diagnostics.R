#!/usr/bin/env Rscript

source(file.path(Sys.getenv("LOWVOLTAGE_PROJECT_ROOT", unset = "."), "R", "00_functions.R"))

imputation <- readRDS(file.path(cache_dir, "malignancy_inclusive_landmark_with_extension_m40_maxit10.rds"))
spline_spec <- readRDS(file.path(cache_dir, "inclusive_spline_specification.rds"))
formulas <- model_formulas(spline_spec, include_malignancy = TRUE)

safe_condition_index <- function(design) {
  no_intercept <- design[, colnames(design) != "(Intercept)", drop = FALSE]
  keep <- apply(no_intercept, 2, sd) > 0
  scaled <- scale(no_intercept[, keep, drop = FALSE])
  singular_values <- svd(scaled, nu = 0, nv = 0)$d
  max(singular_values) / min(singular_values[singular_values > sqrt(.Machine$double.eps)])
}

diagnostic_rows <- vector("list", imputation$m)
vif_rows <- vector("list", imputation$m)
flexible_models <- vector("list", imputation$m)
all_linear_models <- vector("list", imputation$m)

for (index in seq_len(imputation$m)) {
  completed <- complete(imputation, index)
  flexible_fit <- fit_cluster_logistic(completed, formulas$model3)
  linear_fit <- glm(formulas$model3_all_linear, data = completed, family = binomial())
  model <- flexible_fit$model
  covariance <- flexible_fit$covariance
  coefficient <- "low_qrs_landmark"
  conventional_se <- sqrt(vcov(model)[coefficient, coefficient])
  robust_se <- sqrt(covariance[coefficient, coefficient])
  design <- model.matrix(model)
  dfbeta <- dfbetas(model)[, coefficient]
  cooks <- cooks.distance(model)
  leverage <- hatvalues(model)

  diagnostic_rows[[index]] <- tibble(
    imputation = index,
    converged = isTRUE(model$converged),
    iterations = model$iter,
    rank = model$rank,
    parameters = length(coef(model)),
    events = sum(completed$landmark_day28_outcome),
    events_per_parameter = sum(completed$landmark_day28_outcome) / model$rank,
    exposure_or = exp(coef(model)[coefficient]),
    conventional_se = conventional_se,
    cluster_robust_se = robust_se,
    robust_to_conventional_se_ratio = robust_se / conventional_se,
    condition_index_flexible_design = safe_condition_index(design),
    max_abs_exposure_dfbeta = max(abs(dfbeta), na.rm = TRUE),
    dfbeta_above_2_over_sqrt_n = sum(abs(dfbeta) > 2 / sqrt(nrow(completed)), na.rm = TRUE),
    max_cooks_distance = max(cooks, na.rm = TRUE),
    cooks_above_4_over_n = sum(cooks > 4 / nrow(completed), na.rm = TRUE),
    max_leverage = max(leverage, na.rm = TRUE),
    leverage_above_2p_over_n = sum(leverage > 2 * model$rank / nrow(completed), na.rm = TRUE)
  )

  vif_value <- car::vif(linear_fit)
  if (is.matrix(vif_value)) {
    vif_frame <- as.data.frame(vif_value) %>%
      rownames_to_column("term") %>%
      transmute(
        imputation = index,
        term,
        gvif = GVIF,
        degrees_freedom = Df,
        adjusted_gvif = `GVIF^(1/(2*Df))`
      )
  } else {
    vif_frame <- tibble(
      imputation = index,
      term = names(vif_value),
      gvif = as.numeric(vif_value),
      degrees_freedom = 1,
      adjusted_gvif = sqrt(as.numeric(vif_value))
    )
  }
  vif_rows[[index]] <- vif_frame
  flexible_models[[index]] <- model
  all_linear_models[[index]] <- linear_fit
}

diagnostics <- bind_rows(diagnostic_rows)
vif_values <- bind_rows(vif_rows)

summary_table <- tibble(
  diagnostic = c(
    "Model convergence proportion",
    "Fisher-scoring iterations",
    "Model rank",
    "Events per fitted parameter",
    "Maximum adjusted GVIF from all-linear surrogate",
    "Condition index of final flexible design",
    "Cluster-robust/conventional SE ratio",
    "Maximum absolute DFBETA for LQRSV",
    "Maximum Cook's distance",
    "Maximum leverage"
  ),
  minimum = c(
    mean(diagnostics$converged),
    min(diagnostics$iterations),
    min(diagnostics$rank),
    min(diagnostics$events_per_parameter),
    min(tapply(vif_values$adjusted_gvif, vif_values$imputation, max)),
    min(diagnostics$condition_index_flexible_design),
    min(diagnostics$robust_to_conventional_se_ratio),
    min(diagnostics$max_abs_exposure_dfbeta),
    min(diagnostics$max_cooks_distance),
    min(diagnostics$max_leverage)
  ),
  median = c(
    mean(diagnostics$converged),
    median(diagnostics$iterations),
    median(diagnostics$rank),
    median(diagnostics$events_per_parameter),
    median(tapply(vif_values$adjusted_gvif, vif_values$imputation, max)),
    median(diagnostics$condition_index_flexible_design),
    median(diagnostics$robust_to_conventional_se_ratio),
    median(diagnostics$max_abs_exposure_dfbeta),
    median(diagnostics$max_cooks_distance),
    median(diagnostics$max_leverage)
  ),
  maximum = c(
    mean(diagnostics$converged),
    max(diagnostics$iterations),
    max(diagnostics$rank),
    max(diagnostics$events_per_parameter),
    max(tapply(vif_values$adjusted_gvif, vif_values$imputation, max)),
    max(diagnostics$condition_index_flexible_design),
    max(diagnostics$robust_to_conventional_se_ratio),
    max(diagnostics$max_abs_exposure_dfbeta),
    max(diagnostics$max_cooks_distance),
    max(diagnostics$max_leverage)
  )
) %>% mutate(across(where(is.numeric), ~ signif(.x, 6)))

vif_summary <- vif_values %>%
  group_by(term) %>%
  summarise(
    minimum_adjusted_gvif = min(adjusted_gvif),
    median_adjusted_gvif = median(adjusted_gvif),
    maximum_adjusted_gvif = max(adjusted_gvif),
    .groups = "drop"
  ) %>%
  arrange(desc(maximum_adjusted_gvif))

message_time("Performing MI-pooled likelihood comparisons for functional form")
disease_terms <- paste(c(comorbidities, "icd_malignancy"), collapse = " + ")
base_tail <- paste0("gender + ", disease_terms, " + lab_24h_lactate_first + n_ecg_cat")
reduced_formulas <- list(
  global_all_linear = formulas$model3_all_linear,
  age_linear = as.formula(paste0(
    "landmark_day28_outcome ~ low_qrs_landmark + age + gender + ",
    ns_term("weight", spline_spec$weight), " + ", disease_terms, " + ",
    ns_term("sofa", spline_spec$sofa), " + lab_24h_lactate_first + n_ecg_cat"
  )),
  weight_linear = as.formula(paste0(
    "landmark_day28_outcome ~ low_qrs_landmark + ", ns_term("age", spline_spec$age),
    " + gender + weight + ", disease_terms, " + ", ns_term("sofa", spline_spec$sofa),
    " + lab_24h_lactate_first + n_ecg_cat"
  )),
  sofa_linear = as.formula(paste0(
    "landmark_day28_outcome ~ low_qrs_landmark + ", ns_term("age", spline_spec$age),
    " + gender + ", ns_term("weight", spline_spec$weight), " + ", disease_terms,
    " + sofa + lab_24h_lactate_first + n_ecg_cat"
  ))
)

functional_form_tests <- bind_rows(lapply(names(reduced_formulas), function(label) {
  reduced_models <- lapply(seq_len(imputation$m), function(index) {
    glm(reduced_formulas[[label]], data = complete(imputation, index), family = binomial())
  })
  result <- mice::D2(
    as.mira(flexible_models),
    as.mira(reduced_models),
    use = "likelihood"
  )
  as.data.frame(result$result) %>%
    as_tibble() %>%
    mutate(comparison = label, .before = 1)
}))

write_csv(diagnostics, file.path(diagnostics_dir, "final_flexible_model_diagnostics_by_imputation.csv"))
write_csv(summary_table, file.path(results_dir, "final_flexible_model_diagnostics_summary.csv"))
write_csv(vif_values, file.path(diagnostics_dir, "all_linear_surrogate_vif_by_imputation.csv"))
write_csv(vif_summary, file.path(results_dir, "all_linear_surrogate_vif_summary.csv"))
write_csv(functional_form_tests, file.path(results_dir, "pooled_functional_form_D2_tests.csv"))

print(summary_table)
print(vif_summary)
print(functional_form_tests)
message_time("Primary model diagnostics completed")
