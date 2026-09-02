options(stringsAsFactors = FALSE, scipen = 999)

required_packages <- c(
  "broom", "car", "dplyr", "ggplot2", "lattice", "MASS", "mice",
  "mitml", "patchwork", "readr", "sandwich", "splines", "survival",
  "survminer", "tibble", "tidyr"
)
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(mice)
  library(readr)
  library(sandwich)
  library(splines)
  library(survival)
  library(tibble)
  library(tidyr)
})

project_root <- normalizePath(
  Sys.getenv("LOWVOLTAGE_PROJECT_ROOT", unset = "."),
  winslash = "/",
  mustWork = TRUE
)
input_dir <- normalizePath(
  Sys.getenv("LOWVOLTAGE_INPUT_DIR", unset = file.path(project_root, "data", "input")),
  winslash = "/",
  mustWork = TRUE
)
cache_dir <- file.path(project_root, "local_cache")
results_dir <- file.path(project_root, "outputs", "results")
diagnostics_dir <- file.path(project_root, "outputs", "diagnostics")
figures_dir <- file.path(project_root, "outputs", "figures")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diagnostics_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

M <- 40L
MAXIT <- 10L
PMM_DONORS <- 5L

comorbidities <- c(
  "icd_hf", "icd_afib", "icd_renal", "icd_liver", "icd_copd",
  "icd_cad", "icd_stroke"
)
binary_covariates <- c(
  comorbidities, "icd_malignancy", "itvtn_24h_vent_tag",
  "drug_24h_vaso_tag", "itvtn_24h_rrt_tag"
)

message_time <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
}

add_nelson_aalen <- function(data, time_var, event_var, output_var) {
  fit <- coxph(reformulate("1", response = sprintf("Surv(%s, %s)", time_var, event_var)), data = data)
  hazard <- basehaz(fit)
  index <- findInterval(data[[time_var]], hazard$time)
  data[[output_var]] <- ifelse(index == 0L, 0, hazard$hazard[index])
  data
}

as_binary_factor <- function(x) {
  factor(as.integer(x), levels = c(0L, 1L), labels = c("No", "Yes"))
}

prepare_landmark_data <- function(post_path, extension_path = NULL, include_malignancy = FALSE) {
  post <- read_csv(post_path, show_col_types = FALSE) %>%
    mutate(
      stay_id = as.integer(stay_id),
      subject_id = as.integer(subject_id),
      hadm_id = as.integer(hadm_id),
      landmark_los_days = (as.numeric(day28_los) - 24) / 24,
      landmark_day28_outcome = as.integer(day28_outcome),
      low_qrs_landmark = as.integer(low_qrs),
      n_ecg_pm24 = as.integer(n_ecg_in_window),
      n_ecg_cat = factor(
        case_when(n_ecg_pm24 == 1L ~ "1", n_ecg_pm24 == 2L ~ "2", TRUE ~ ">=3"),
        levels = c("1", "2", ">=3")
      ),
      gender = factor(gender),
      weight = if_else(as.numeric(weight) >= 20 & as.numeric(weight) <= 300, as.numeric(weight), NA_real_),
      age = as.numeric(age),
      sofa = as.numeric(sofa),
      sapsii = as.numeric(sapsii),
      lab_24h_lactate_first = as.numeric(lab_24h_lactate_first),
      lab_24h_creatinine_first = as.numeric(lab_24h_creatinine_first),
      lab_24h_ph_first = as.numeric(lab_24h_ph_first),
      anchor_year_group = factor(anchor_year_group),
      first_careunit = factor(first_careunit)
    ) %>%
    filter(landmark_los_days > 0, !is.na(landmark_day28_outcome), n_ecg_pm24 >= 1L)

  for (variable in intersect(binary_covariates, names(post))) {
    post[[variable]] <- as_binary_factor(post[[variable]])
  }

  if (!is.null(extension_path)) {
    extension <- read_csv(extension_path, show_col_types = FALSE) %>%
      filter(landmark_selection_group == "included_landmark") %>%
      transmute(
        stay_id = as.integer(stay_id),
        first_ecg_low_qrs = as.integer(first_post_ecg_low_qrs),
        sepsis_hours_from_icu = as.numeric(crtr_sepsis3_suspected_infection_hours_from_icu)
      )
    post <- post %>% left_join(extension, by = "stay_id")
    stopifnot(!anyNA(post$first_ecg_low_qrs))
  }

  if (!include_malignancy) {
    stopifnot(all(as.character(post$icd_malignancy) == "No"))
  }

  post <- post %>%
    arrange(subject_id, intime, stay_id) %>%
    group_by(subject_id) %>%
    mutate(first_patient_flag = as.integer(row_number() == 1L)) %>%
    ungroup()

  post <- add_nelson_aalen(
    post,
    "landmark_los_days",
    "landmark_day28_outcome",
    "nelson_aalen_landmark"
  )
  post %>% arrange(stay_id)
}

prepare_broad_data <- function(broad_path, include_malignancy = FALSE) {
  data <- read_csv(broad_path, show_col_types = FALSE) %>%
    mutate(
      stay_id = as.integer(stay_id),
      subject_id = as.integer(subject_id),
      hadm_id = as.integer(hadm_id),
      landmark_los_days = as.numeric(day28_los) / 24,
      landmark_day28_outcome = as.integer(day28_outcome),
      low_qrs_landmark = as.integer(low_qrs),
      n_ecg_pm24 = as.integer(n_ecg_in_window),
      n_ecg_cat = factor(
        case_when(n_ecg_pm24 == 1L ~ "1", n_ecg_pm24 == 2L ~ "2", TRUE ~ ">=3"),
        levels = c("1", "2", ">=3")
      ),
      gender = factor(gender),
      weight = if_else(as.numeric(weight) >= 20 & as.numeric(weight) <= 300, as.numeric(weight), NA_real_),
      age = as.numeric(age),
      sofa = as.numeric(sofa),
      sapsii = as.numeric(sapsii),
      lab_24h_lactate_first = as.numeric(lab_24h_lactate_first),
      lab_24h_creatinine_first = as.numeric(lab_24h_creatinine_first),
      lab_24h_ph_first = as.numeric(lab_24h_ph_first),
      anchor_year_group = factor(anchor_year_group),
      first_careunit = factor(first_careunit)
    ) %>%
    filter(landmark_los_days > 0, !is.na(landmark_day28_outcome), n_ecg_pm24 >= 1L)

  for (variable in intersect(binary_covariates, names(data))) {
    data[[variable]] <- as_binary_factor(data[[variable]])
  }
  if (!include_malignancy) stopifnot(all(as.character(data$icd_malignancy) == "No"))
  data <- data %>%
    arrange(subject_id, intime, stay_id) %>%
    group_by(subject_id) %>%
    mutate(first_patient_flag = as.integer(row_number() == 1L)) %>%
    ungroup()
  data <- add_nelson_aalen(data, "landmark_los_days", "landmark_day28_outcome", "nelson_aalen_landmark")
  data %>% arrange(stay_id)
}

fixed_spline_specification <- function(data) {
  make_spec <- function(x) {
    observed <- x[is.finite(x)]
    list(
      knots = unname(quantile(observed, probs = c(1 / 3, 2 / 3), type = 7)),
      boundary = range(observed)
    )
  }
  list(
    age = make_spec(data$age),
    weight = make_spec(data$weight),
    sofa = make_spec(data$sofa)
  )
}

format_num <- function(x) formatC(x, digits = 12, format = "fg", flag = "#")

ns_term <- function(variable, specification) {
  paste0(
    "ns(", variable,
    ", knots=c(", paste(format_num(specification$knots), collapse = ","), ")",
    ", Boundary.knots=c(", paste(format_num(specification$boundary), collapse = ","), "))"
  )
}

model_formulas <- function(spline_spec, exposure = "low_qrs_landmark", include_malignancy = FALSE) {
  disease_terms <- c(comorbidities, if (include_malignancy) "icd_malignancy")
  demographic_terms <- c(
    ns_term("age", spline_spec$age), "gender", ns_term("weight", spline_spec$weight), disease_terms
  )
  list(
    model1 = reformulate(exposure, response = "landmark_day28_outcome"),
    model2 = reformulate(c(exposure, demographic_terms), response = "landmark_day28_outcome"),
    model3 = reformulate(
      c(exposure, demographic_terms, ns_term("sofa", spline_spec$sofa), "lab_24h_lactate_first", "n_ecg_cat"),
      response = "landmark_day28_outcome"
    ),
    model3_continuous_count = reformulate(
      c(exposure, demographic_terms, ns_term("sofa", spline_spec$sofa), "lab_24h_lactate_first", "n_ecg_pm24"),
      response = "landmark_day28_outcome"
    ),
    model3_no_ecg_count = reformulate(
      c(exposure, demographic_terms, ns_term("sofa", spline_spec$sofa), "lab_24h_lactate_first"),
      response = "landmark_day28_outcome"
    ),
    model3_no_lactate = reformulate(
      c(exposure, demographic_terms, ns_term("sofa", spline_spec$sofa), "n_ecg_cat"),
      response = "landmark_day28_outcome"
    ),
    model3_all_linear = reformulate(
      c(exposure, "age", "gender", "weight", disease_terms, "sofa", "lab_24h_lactate_first", "n_ecg_cat"),
      response = "landmark_day28_outcome"
    ),
    expanded = reformulate(
      c(
        exposure, ns_term("age", spline_spec$age), "gender", ns_term("weight", spline_spec$weight),
        ns_term("sofa", spline_spec$sofa), "sapsii", "itvtn_24h_vent_tag",
        "drug_24h_vaso_tag", "itvtn_24h_rrt_tag", "lab_24h_lactate_first",
        "lab_24h_creatinine_first", "lab_24h_ph_first", "n_ecg_cat"
      ),
      response = "landmark_day28_outcome"
    )
  )
}

make_landmark_mids <- function(data, cache_path, seed, include_malignancy = FALSE) {
  fixed_variables <- c(
    "stay_id", "subject_id", "hadm_id", "landmark_los_days",
    "landmark_day28_outcome", "nelson_aalen_landmark", "low_qrs_landmark",
    "n_ecg_pm24", "n_ecg_cat", "anchor_year_group", "first_careunit",
    "first_patient_flag",
    if ("first_ecg_low_qrs" %in% names(data)) "first_ecg_low_qrs",
    if ("sepsis_hours_from_icu" %in% names(data)) "sepsis_hours_from_icu",
    if (include_malignancy) "icd_malignancy"
  )
  analysis_variables <- c(
    "age", "gender", "weight", comorbidities,
    if (include_malignancy) "icd_malignancy",
    "sofa", "lab_24h_lactate_first", "sapsii",
    "itvtn_24h_vent_tag", "drug_24h_vaso_tag", "itvtn_24h_rrt_tag",
    "lab_24h_creatinine_first", "lab_24h_ph_first"
  )
  variables <- unique(c(fixed_variables, analysis_variables))
  imputation_data <- data[, variables, drop = FALSE]

  method <- make.method(imputation_data)
  method[] <- ""
  imputed_variables <- intersect(
    c("weight", "lab_24h_lactate_first", "lab_24h_creatinine_first", "lab_24h_ph_first"),
    names(imputation_data)
  )
  method[imputed_variables] <- "pmm"

  predictors <- make.predictorMatrix(imputation_data)
  predictors[,] <- 0L
  predictor_candidates <- setdiff(
    names(imputation_data),
    c("stay_id", "subject_id", "hadm_id", "first_ecg_low_qrs", "sepsis_hours_from_icu")
  )
  for (target in imputed_variables) {
    predictors[target, setdiff(predictor_candidates, target)] <- 1L
  }
  diag(predictors) <- 0L

  write_csv(
    as.data.frame(predictors) %>% rownames_to_column("target_variable"),
    file.path(diagnostics_dir, paste0(tools::file_path_sans_ext(basename(cache_path)), "_predictor_matrix.csv"))
  )
  write_csv(
    tibble(variable = names(method), method = unname(method)),
    file.path(diagnostics_dir, paste0(tools::file_path_sans_ext(basename(cache_path)), "_methods.csv"))
  )

  if (file.exists(cache_path) && Sys.getenv("LOWVOLTAGE_FORCE_IMPUTE", "0") != "1") {
    message_time("Reading cached imputation object: ", cache_path)
    return(readRDS(cache_path))
  }

  message_time("Running MICE: n=", nrow(imputation_data), ", m=", M, ", maxit=", MAXIT)
  imputation <- mice(
    imputation_data,
    m = M,
    maxit = MAXIT,
    method = method,
    predictorMatrix = predictors,
    donors = PMM_DONORS,
    seed = seed,
    printFlag = FALSE
  )
  saveRDS(imputation, cache_path)
  imputation
}

fit_cluster_logistic <- function(data, formula) {
  model <- glm(formula, data = data, family = binomial())
  covariance <- vcovCL(model, cluster = data$subject_id, type = "HC0", fix = TRUE)
  list(model = model, covariance = covariance)
}

pool_rubin_scalar <- function(q, u, n_clusters, parameter_count, transform = c("identity", "exp")) {
  transform <- match.arg(transform)
  pooled <- mice::pool.scalar(
    Q = q,
    U = u,
    n = n_clusters,
    k = parameter_count,
    rule = "rubin1987"
  )
  critical <- qt(0.975, df = pooled$df)
  estimate <- pooled$qbar
  lower <- pooled$qbar - critical * sqrt(pooled$t)
  upper <- pooled$qbar + critical * sqrt(pooled$t)
  if (transform == "exp") {
    estimate <- exp(estimate)
    lower <- exp(lower)
    upper <- exp(upper)
  }
  tibble(
    estimate = estimate,
    conf_low = lower,
    conf_high = upper,
    p_value = 2 * pt(abs(pooled$qbar / sqrt(pooled$t)), df = pooled$df, lower.tail = FALSE),
    df = pooled$df,
    within_variance = pooled$ubar,
    between_variance = pooled$b,
    total_variance = pooled$t,
    relative_increase_variance = pooled$r,
    fraction_missing_information = pooled$fmi,
    monte_carlo_error_fraction_of_se = sqrt(pooled$b / length(q)) / sqrt(pooled$t)
  )
}

fit_mi_logistic <- function(imputation, formula, exposure, label, subset_rows = NULL, transform_data = identity) {
  q <- u <- numeric(imputation$m)
  n_values <- events <- exposed <- clusters <- ranks <- integer(imputation$m)
  models <- vector("list", imputation$m)
  for (index in seq_len(imputation$m)) {
    completed <- transform_data(complete(imputation, index))
    if (!is.null(subset_rows)) completed <- subset_rows(completed)
    fitted <- fit_cluster_logistic(completed, formula)
    q[index] <- coef(fitted$model)[exposure]
    u[index] <- fitted$covariance[exposure, exposure]
    n_values[index] <- nrow(completed)
    events[index] <- sum(completed$landmark_day28_outcome)
    exposed[index] <- sum(as.integer(completed[[exposure]]) == 1L)
    clusters[index] <- n_distinct(completed$subject_id)
    ranks[index] <- fitted$model$rank
    models[[index]] <- fitted$model
  }
  stopifnot(length(unique(n_values)) == 1L, length(unique(events)) == 1L)
  pooled <- pool_rubin_scalar(q, u, min(clusters), max(ranks), "exp")
  pooled %>%
    mutate(
      analysis = label,
      n = unique(n_values),
      exposed = unique(exposed),
      deaths = unique(events),
      patients = min(clusters),
      parameters = max(ranks),
      m = imputation$m,
      .before = 1
    )
}

cox_formula_from_logistic <- function(logistic_formula) {
  rhs <- attr(terms(logistic_formula), "term.labels")
  as.formula(paste(
    "Surv(landmark_los_days, landmark_day28_outcome) ~",
    paste(rhs, collapse = " + ")
  ))
}

fit_mi_cox <- function(imputation, logistic_formula, exposure, label) {
  formula <- cox_formula_from_logistic(logistic_formula)
  q <- u <- numeric(imputation$m)
  clusters <- ranks <- integer(imputation$m)
  ph_rows <- vector("list", imputation$m)
  for (index in seq_len(imputation$m)) {
    completed <- complete(imputation, index)
    model <- coxph(
      formula,
      data = completed,
      cluster = subject_id,
      robust = TRUE,
      x = TRUE,
      model = TRUE
    )
    q[index] <- coef(model)[exposure]
    u[index] <- vcov(model)[exposure, exposure]
    clusters[index] <- n_distinct(completed$subject_id)
    ranks[index] <- length(coef(model))
    ph <- as.data.frame(cox.zph(model, transform = "km")$table) %>%
      rownames_to_column("term") %>%
      as_tibble()
    ph_rows[[index]] <- ph %>% mutate(imputation = index, .before = 1)
  }
  pooled <- pool_rubin_scalar(q, u, min(clusters), max(ranks), "exp") %>%
    mutate(
      analysis = label,
      n = nrow(imputation$data),
      exposed = sum(imputation$data[[exposure]]),
      deaths = sum(imputation$data$landmark_day28_outcome),
      patients = min(clusters),
      parameters = max(ranks),
      m = imputation$m,
      .before = 1
    )
  list(effect = pooled, ph = bind_rows(ph_rows))
}

pool_marginal_standardization <- function(imputation, formula, exposure = "low_qrs_landmark") {
  rows <- vector("list", imputation$m)
  clusters <- ranks <- integer(imputation$m)
  for (index in seq_len(imputation$m)) {
    completed <- complete(imputation, index)
    fitted <- fit_cluster_logistic(completed, formula)
    model <- fitted$model
    covariance <- fitted$covariance
    data0 <- completed
    data1 <- completed
    data0[[exposure]] <- 0L
    data1[[exposure]] <- 1L
    terms_without_outcome <- delete.response(terms(model))
    x0 <- model.matrix(terms_without_outcome, data0)
    x1 <- model.matrix(terms_without_outcome, data1)
    p0 <- plogis(drop(x0 %*% coef(model)))
    p1 <- plogis(drop(x1 %*% coef(model)))
    risk0 <- mean(p0)
    risk1 <- mean(p1)
    gradient0 <- colMeans(x0 * as.numeric(p0 * (1 - p0)))
    gradient1 <- colMeans(x1 * as.numeric(p1 * (1 - p1)))
    variance0 <- drop(t(gradient0) %*% covariance %*% gradient0)
    variance1 <- drop(t(gradient1) %*% covariance %*% gradient1)
    covariance01 <- drop(t(gradient0) %*% covariance %*% gradient1)
    rows[[index]] <- tibble(
      risk0 = risk0,
      risk1 = risk1,
      variance0 = variance0,
      variance1 = variance1,
      covariance01 = covariance01,
      risk_difference = risk1 - risk0,
      variance_risk_difference = variance1 + variance0 - 2 * covariance01,
      log_risk_ratio = log(risk1 / risk0),
      variance_log_risk_ratio = variance1 / risk1^2 + variance0 / risk0^2 -
        2 * covariance01 / (risk1 * risk0)
    )
    clusters[index] <- n_distinct(completed$subject_id)
    ranks[index] <- model$rank
  }
  values <- bind_rows(rows)
  n_clusters <- min(clusters)
  parameter_count <- max(ranks)

  logit_pool <- function(estimates, variances) {
    transformed <- qlogis(estimates)
    transformed_variance <- variances / (estimates^2 * (1 - estimates)^2)
    pooled <- pool_rubin_scalar(transformed, transformed_variance, n_clusters, parameter_count)
    pooled %>% mutate(across(c(estimate, conf_low, conf_high), plogis))
  }
  risk0 <- logit_pool(values$risk0, values$variance0)
  risk1 <- logit_pool(values$risk1, values$variance1)
  rd <- pool_rubin_scalar(
    values$risk_difference,
    values$variance_risk_difference,
    n_clusters,
    parameter_count
  )
  rr <- pool_rubin_scalar(
    values$log_risk_ratio,
    values$variance_log_risk_ratio,
    n_clusters,
    parameter_count,
    "exp"
  )
  list(
    by_imputation = values,
    pooled = tibble(
      risk_unexposed = risk0$estimate,
      risk_unexposed_low = risk0$conf_low,
      risk_unexposed_high = risk0$conf_high,
      risk_exposed = risk1$estimate,
      risk_exposed_low = risk1$conf_low,
      risk_exposed_high = risk1$conf_high,
      risk_difference = rd$estimate,
      risk_difference_low = rd$conf_low,
      risk_difference_high = rd$conf_high,
      risk_ratio = rr$estimate,
      risk_ratio_low = rr$conf_low,
      risk_ratio_high = rr$conf_high
    )
  )
}

write_session_information <- function(path) {
  writeLines(capture.output(sessionInfo()), path, useBytes = TRUE)
}
