#!/usr/bin/env Rscript

source(file.path(Sys.getenv("LOWVOLTAGE_PROJECT_ROOT", unset = "."), "R", "00_functions.R"))

primary <- prepare_landmark_data(
  file.path(input_dir, "malignancy_inclusive_post_icu_24h.csv"),
  file.path(input_dir, "landmark_extension_inclusive.csv"),
  include_malignancy = TRUE
)
primary <- primary %>% mutate(exposure_group = factor(low_qrs_landmark, levels = c(0, 1), labels = c("No machine-reported LQRSV", "Machine-reported LQRSV")))

standardized_mean_difference_continuous <- function(x, group) {
  x0 <- x[group == 0 & is.finite(x)]
  x1 <- x[group == 1 & is.finite(x)]
  (mean(x1) - mean(x0)) / sqrt((var(x1) + var(x0)) / 2)
}

standardized_mean_difference_binary <- function(x, group) {
  value <- as.integer(as.character(x) %in% c("1", "Yes", "F", "TRUE"))
  p0 <- mean(value[group == 0], na.rm = TRUE)
  p1 <- mean(value[group == 1], na.rm = TRUE)
  (p1 - p0) / sqrt((p1 * (1 - p1) + p0 * (1 - p0)) / 2)
}

continuous_row <- function(data, variable, label, digits = 1, group_variable = "low_qrs_landmark") {
  group <- data[[group_variable]]
  display <- function(level) {
    values <- data[[variable]][group == level]
    values <- values[is.finite(values)]
    paste0(
      formatC(median(values), format = "f", digits = digits), " [",
      formatC(quantile(values, 0.25), format = "f", digits = digits), "-",
      formatC(quantile(values, 0.75), format = "f", digits = digits), "]"
    )
  }
  tibble(
    characteristic = label,
    group_0 = display(0),
    group_1 = display(1),
    smd = standardized_mean_difference_continuous(data[[variable]], group),
    missing_n = sum(is.na(data[[variable]])),
    missing_percent = mean(is.na(data[[variable]])) * 100
  )
}

binary_row <- function(data, variable, label, positive = "Yes", group_variable = "low_qrs_landmark") {
  group <- data[[group_variable]]
  positive_flag <- if (variable == "gender") {
    as.integer(as.character(data[[variable]]) == positive)
  } else {
    as.integer(as.character(data[[variable]]) %in% c(positive, "1", "TRUE"))
  }
  display <- function(level) {
    denominator <- sum(group == level)
    numerator <- sum(positive_flag[group == level] == 1L, na.rm = TRUE)
    sprintf("%s (%.1f%%)", format(numerator, big.mark = ","), numerator / denominator * 100)
  }
  tibble(
    characteristic = label,
    group_0 = display(0),
    group_1 = display(1),
    smd = standardized_mean_difference_binary(
      if (variable == "gender") factor(if_else(as.character(data[[variable]]) == positive, "Yes", "No")) else data[[variable]],
      group
    ),
    missing_n = sum(is.na(data[[variable]])),
    missing_percent = mean(is.na(data[[variable]])) * 100
  )
}

table1 <- bind_rows(
  continuous_row(primary, "age", "Age, years", 1),
  binary_row(primary, "gender", "Female sex", positive = "F"),
  continuous_row(primary, "weight", "Weight, kg", 1),
  continuous_row(primary, "sofa", "SOFA score", 0),
  continuous_row(primary, "sapsii", "SAPS II", 0),
  continuous_row(primary, "lab_24h_lactate_first", "First lactate, mmol/L", 1),
  binary_row(primary, "icd_hf", "Heart failure"),
  binary_row(primary, "icd_afib", "Atrial fibrillation"),
  binary_row(primary, "icd_renal", "Renal disease"),
  binary_row(primary, "icd_liver", "Liver disease"),
  binary_row(primary, "icd_copd", "Chronic pulmonary disease"),
  binary_row(primary, "icd_cad", "Coronary artery disease"),
  binary_row(primary, "icd_stroke", "Cerebrovascular disease"),
  binary_row(primary, "icd_malignancy", "Malignancy"),
  binary_row(primary, "itvtn_24h_vent_tag", "Mechanical ventilation"),
  binary_row(primary, "drug_24h_vaso_tag", "Vasopressor use"),
  binary_row(primary, "itvtn_24h_rrt_tag", "Renal replacement therapy")
) %>% mutate(smd = round(smd, 3), missing_percent = round(missing_percent, 2))

write_csv(table1, file.path(results_dir, "Table_1_primary_baseline_by_LQRSV.csv"))

selection <- read_csv(file.path(input_dir, "landmark_extension_inclusive.csv"), show_col_types = FALSE) %>%
  filter(landmark_selection_group %in% c("included_landmark", "alive_24h_no_post_ecg")) %>%
  mutate(
    comparison_group = as.integer(landmark_selection_group == "included_landmark"),
    gender = factor(gender),
    weight = as.numeric(weight),
    age = as.numeric(age),
    sofa = as.numeric(sofa),
    sapsii = as.numeric(sapsii),
    lab_24h_lactate_first = as.numeric(lab_24h_lactate_first)
  )
for (variable in intersect(binary_covariates, names(selection))) selection[[variable]] <- as_binary_factor(selection[[variable]])

table_selection <- bind_rows(
  continuous_row(selection, "age", "Age, years", 1, "comparison_group"),
  binary_row(selection, "gender", "Female sex", "F", "comparison_group"),
  continuous_row(selection, "sofa", "SOFA score", 0, "comparison_group"),
  continuous_row(selection, "sapsii", "SAPS II", 0, "comparison_group"),
  continuous_row(selection, "lab_24h_lactate_first", "First lactate, mmol/L", 1, "comparison_group"),
  binary_row(selection, "icd_hf", "Heart failure", "Yes", "comparison_group"),
  binary_row(selection, "icd_afib", "Atrial fibrillation", "Yes", "comparison_group"),
  binary_row(selection, "icd_malignancy", "Malignancy", "Yes", "comparison_group"),
  binary_row(selection, "itvtn_24h_vent_tag", "Mechanical ventilation", "Yes", "comparison_group"),
  binary_row(selection, "drug_24h_vaso_tag", "Vasopressor use", "Yes", "comparison_group"),
  binary_row(selection, "day28_outcome", "Death through day 28", "1", "comparison_group")
) %>%
  rename(alive_24h_no_post_ecg = group_0, included_landmark = group_1) %>%
  mutate(smd = round(smd, 3), missing_percent = round(missing_percent, 2))
write_csv(table_selection, file.path(results_dir, "Table_S_selection_included_vs_no_ECG.csv"))

missing_variables <- c(
  "age", "gender", "weight", comorbidities, "icd_malignancy", "sofa",
  "lab_24h_lactate_first", "sapsii", "itvtn_24h_vent_tag",
  "drug_24h_vaso_tag", "itvtn_24h_rrt_tag", "lab_24h_creatinine_first",
  "lab_24h_ph_first"
)
missingness <- bind_rows(lapply(missing_variables, function(variable) {
  tibble(
    variable = variable,
    missing_n = sum(is.na(primary[[variable]])),
    missing_percent = mean(is.na(primary[[variable]])) * 100
  )
})) %>% mutate(missing_percent = round(missing_percent, 2))
write_csv(missingness, file.path(results_dir, "Table_S_primary_covariate_missingness.csv"))

model_results <- read_csv(file.path(results_dir, "final_landmark_model_results.csv"), show_col_types = FALSE)
broad_result <- read_csv(file.path(results_dir, "malignancy_inclusive_broad_logistic.csv"), show_col_types = FALSE)
cox_results <- read_csv(file.path(results_dir, "secondary_cox_model_results.csv"), show_col_types = FALSE)

primary_model_table <- model_results %>%
  filter(analysis %in% c(
    "Malignancy-inclusive Model 1",
    "Malignancy-inclusive Model 2",
    "Malignancy-inclusive flexible Model 3"
  )) %>%
  select(analysis, n, exposed, deaths, estimate, conf_low, conf_high, p_value, df, fraction_missing_information)
write_csv(primary_model_table, file.path(results_dir, "Table_3_primary_sequential_models.csv"))

sensitivity_labels <- c(
  "Malignancy-inclusive flexible Model 3",
  "Malignancy-inclusive continuous ECG-count adjustment",
  "Malignancy-inclusive exposure defined by first post-ICU ECG",
  "Malignancy-inclusive exactly one ECG in first 24 hours",
  "Malignancy-inclusive first eligible landmark stay per patient",
  "Malignancy-inclusive flexible Model 3 without lactate",
  "Malignancy-inclusive log-transformed lactate model",
  "Malignancy-inclusive Sepsis-3 timing within 24 hours of ICU admission",
  "Malignancy-inclusive expanded acuity model",
  "Malignancy-inclusive complete-case flexible Model 3",
  "Model 3: flexible primary model with categorical ECG count"
)
sensitivity_table <- bind_rows(
  model_results %>% filter(analysis %in% sensitivity_labels),
  broad_result
) %>%
  select(analysis, n, exposed, deaths, estimate, conf_low, conf_high, p_value, df, fraction_missing_information)
write_csv(sensitivity_table, file.path(results_dir, "Table_4_primary_sensitivity_models.csv"))
write_csv(cox_results, file.path(results_dir, "Table_S_secondary_Cox_models.csv"))

forest_data <- sensitivity_table %>%
  mutate(
    label = recode(
      analysis,
      "Malignancy-inclusive flexible Model 3" = "Primary landmark analysis",
      "Malignancy-inclusive continuous ECG-count adjustment" = "Continuous ECG count",
      "Malignancy-inclusive exposure defined by first post-ICU ECG" = "First post-ICU ECG",
      "Malignancy-inclusive exactly one ECG in first 24 hours" = "Exactly one ECG",
      "Malignancy-inclusive first eligible landmark stay per patient" = "First eligible stay per patient",
      "Malignancy-inclusive flexible Model 3 without lactate" = "Model without lactate",
      "Malignancy-inclusive log-transformed lactate model" = "Log-transformed lactate",
      "Malignancy-inclusive Sepsis-3 timing within 24 hours of ICU admission" = "Sepsis timing within +/-24 h",
      "Malignancy-inclusive expanded acuity model" = "Expanded acuity model",
      "Malignancy-inclusive complete-case flexible Model 3" = "Complete-case analysis",
      "Model 3: flexible primary model with categorical ECG count" = "Non-oncologic cohort",
      "Malignancy-inclusive broad-window flexible Model 3" = "Broad ECG window"
    ),
    label = factor(label, levels = rev(c(
      "Primary landmark analysis", "Continuous ECG count", "First post-ICU ECG",
      "Exactly one ECG", "First eligible stay per patient", "Model without lactate",
      "Log-transformed lactate", "Sepsis timing within +/-24 h", "Expanded acuity model",
      "Complete-case analysis",
      "Non-oncologic cohort", "Broad ECG window"
    ))),
    estimate_text = sprintf("%.3f (%.3f-%.3f)", estimate, conf_low, conf_high)
  )
write_csv(forest_data, file.path(results_dir, "Figure_3_forest_data.csv"))

forest_plot <- ggplot(forest_data, aes(estimate, label)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.45, colour = "#666666") +
  geom_errorbar(
    aes(xmin = conf_low, xmax = conf_high),
    width = 0.18,
    linewidth = 0.55,
    orientation = "y"
  ) +
  geom_point(shape = 15, size = 2.5) +
  geom_text(aes(x = 1.55, label = estimate_text), hjust = 0, size = 3.2) +
  coord_cartesian(xlim = c(0.90, 1.90), clip = "off") +
  scale_x_continuous(breaks = c(1.0, 1.2, 1.4, 1.6, 1.8)) +
  labs(x = "Adjusted odds ratio (95% CI)", y = NULL) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.y = element_text(colour = "#111111"),
    plot.margin = margin(8, 120, 8, 8)
  )
ggsave(
  file.path(figures_dir, "Figure_3_adjusted_associations.tif"),
  forest_plot,
  width = 8.2,
  height = 7.4,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

survival_fit <- survfit(
  Surv(landmark_los_days, landmark_day28_outcome) ~ exposure_group,
  data = primary
)
group_levels <- levels(primary$exposure_group)
km_curve_data <- survminer::surv_summary(survival_fit, data = primary) %>%
  mutate(
    exposure_group = factor(sub("^exposure_group=", "", strata), levels = group_levels),
    icu_day = time + 1
  )
km_plot <- ggplot(
  km_curve_data,
  aes(icu_day, surv, colour = exposure_group, linetype = exposure_group)
) +
  geom_step(linewidth = 0.65) +
  geom_point(
    data = filter(km_curve_data, n.censor > 0),
    shape = 3,
    size = 1.0,
    stroke = 0.45
  ) +
  scale_colour_manual(values = c("#111111", "#777777")) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  scale_x_continuous(breaks = c(1, 7, 14, 21, 28), limits = c(1, 28)) +
  coord_cartesian(ylim = c(0.70, 1.00)) +
  labs(x = NULL, y = "Survival probability", colour = NULL, linetype = NULL) +
  theme_classic(base_size = 10) +
  theme(
    legend.position = "top",
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  )

risk_times <- c(0, 6, 13, 20, 27)
risk_days <- c(1, 7, 14, 21, 28)
risk_summary <- summary(survival_fit, times = risk_times, extend = TRUE)
risk_data <- tibble(
  time = risk_summary$time,
  n_risk = risk_summary$n.risk,
  strata = as.character(risk_summary$strata)
) %>%
  mutate(
    exposure_group = factor(sub("^exposure_group=", "", strata), levels = rev(group_levels)),
    icu_day = risk_days[match(time, risk_times)]
  )
risk_plot <- ggplot(risk_data, aes(icu_day, exposure_group, label = format(n_risk, big.mark = ","))) +
  geom_text(size = 3.0, colour = "#111111") +
  scale_x_continuous(breaks = risk_days, limits = c(1, 28)) +
  labs(x = "Day after ICU admission", y = "Number at risk") +
  theme_classic(base_size = 9) +
  theme(
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_text(margin = margin(r = 8))
  )
km_combined <- patchwork::wrap_plots(km_plot, risk_plot, ncol = 1, heights = c(3.2, 1.15))
ggsave(
  file.path(figures_dir, "Figure_2_landmark_Kaplan_Meier.tif"),
  km_combined,
  width = 7.4,
  height = 6.3,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

write_csv(
  primary %>%
    count(n_ecg_cat, low_qrs_landmark, name = "stays") %>%
    group_by(n_ecg_cat) %>%
    mutate(percent_within_ecg_count = stays / sum(stays) * 100) %>%
    ungroup(),
  file.path(results_dir, "Table_S_ECG_count_distribution.csv")
)

message_time("Tables and figures completed")
