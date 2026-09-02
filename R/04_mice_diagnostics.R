#!/usr/bin/env Rscript

source(file.path(Sys.getenv("LOWVOLTAGE_PROJECT_ROOT", unset = "."), "R", "00_functions.R"))

imputation <- readRDS(file.path(cache_dir, "malignancy_inclusive_landmark_with_extension_m40_maxit10.rds"))
imputed_variables <- names(imputation$method)[imputation$method != ""]

convergence_values <- mice::convergence(imputation) %>%
  filter(vrb %in% imputed_variables)
final_convergence <- convergence_values %>%
  group_by(vrb) %>%
  filter(.it == max(.it)) %>%
  ungroup() %>%
  transmute(
    variable = vrb,
    iteration = .it,
    autocorrelation = ac,
    potential_scale_reduction_factor = psrf
  )

write_csv(convergence_values, file.path(diagnostics_dir, "primary_landmark_mice_convergence_by_iteration.csv"))
write_csv(final_convergence, file.path(results_dir, "primary_landmark_mice_convergence_summary.csv"))

if (is.null(imputation$loggedEvents) || nrow(imputation$loggedEvents) == 0L) {
  write_csv(
    tibble(status = "No logged MICE events"),
    file.path(diagnostics_dir, "primary_landmark_mice_logged_events.csv")
  )
} else {
  write_csv(
    as_tibble(imputation$loggedEvents),
    file.path(diagnostics_dir, "primary_landmark_mice_logged_events.csv")
  )
}

chain_rows <- list()
row_index <- 1L
for (variable in imputed_variables) {
  means <- imputation$chainMean[variable, , , drop = FALSE]
  for (chain in seq_len(dim(means)[3])) {
    chain_values <- as.numeric(means[1, seq_len(dim(means)[2]), chain, drop = TRUE])
    chain_rows[[row_index]] <- tibble(
      variable = variable,
      iteration = seq_along(chain_values),
      chain = factor(rep(chain, length(chain_values))),
      chain_mean = chain_values
    )
    row_index <- row_index + 1L
  }
}
chain_data <- bind_rows(chain_rows)
write_csv(chain_data, file.path(diagnostics_dir, "primary_landmark_mice_chain_means.csv"))

variable_labels <- c(
  weight = "Weight",
  lab_24h_lactate_first = "First lactate",
  lab_24h_creatinine_first = "First creatinine",
  lab_24h_ph_first = "First pH"
)
chain_plot_data <- chain_data %>%
  mutate(variable_label = recode(variable, !!!variable_labels))

trace_plot <- ggplot(chain_plot_data, aes(iteration, chain_mean, group = chain)) +
  geom_line(linewidth = 0.25, alpha = 0.35, colour = "#303030") +
  facet_wrap(~ variable_label, scales = "free_y", ncol = 2) +
  scale_x_continuous(breaks = seq_len(MAXIT)) +
  labs(x = "Iteration", y = "Chain mean") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
ggsave(
  file.path(figures_dir, "Supplementary_Figure_S1_primary_landmark_MICE_trace.tif"),
  trace_plot,
  width = 7.2,
  height = 5.8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

convergence_plot <- convergence_values %>%
  select(.it, vrb, ac, psrf) %>%
  pivot_longer(c(ac, psrf), names_to = "diagnostic", values_to = "value") %>%
  filter(is.finite(value)) %>%
  mutate(
    variable_label = recode(vrb, !!!variable_labels),
    diagnostic = recode(
      diagnostic,
      ac = "Lag-one autocorrelation",
      psrf = "Potential scale reduction factor"
    )
  ) %>%
  ggplot(aes(.it, value)) +
  geom_hline(
    data = tibble(
      diagnostic = c("Lag-one autocorrelation", "Potential scale reduction factor"),
      reference = c(0, 1)
    ),
    aes(yintercept = reference),
    linetype = 2,
    colour = "#777777",
    inherit.aes = FALSE
  ) +
  geom_line(linewidth = 0.5, colour = "#202020") +
  facet_grid(diagnostic ~ variable_label, scales = "free_y") +
  labs(x = "Iteration", y = NULL) +
  theme_classic(base_size = 9) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
ggsave(
  file.path(figures_dir, "Supplementary_Figure_S2_primary_landmark_MICE_convergence.tif"),
  convergence_plot,
  width = 8.4,
  height = 4.8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

distribution_rows <- list()
row_index <- 1L
for (variable in imputed_variables) {
  observed <- imputation$data[[variable]]
  observed <- observed[!is.na(observed)]
  distribution_rows[[row_index]] <- tibble(
    variable = variable,
    source = "Observed",
    imputation = 0L,
    n = length(observed),
    mean = mean(observed),
    sd = sd(observed),
    q25 = unname(quantile(observed, 0.25)),
    median = median(observed),
    q75 = unname(quantile(observed, 0.75))
  )
  row_index <- row_index + 1L
  for (index in seq_len(imputation$m)) {
    values <- complete(imputation, index)[[variable]]
    distribution_rows[[row_index]] <- tibble(
      variable = variable,
      source = "Completed",
      imputation = index,
      n = length(values),
      mean = mean(values),
      sd = sd(values),
      q25 = unname(quantile(values, 0.25)),
      median = median(values),
      q75 = unname(quantile(values, 0.75))
    )
    row_index <- row_index + 1L
  }
}
distribution_summary <- bind_rows(distribution_rows)
write_csv(distribution_summary, file.path(diagnostics_dir, "primary_landmark_observed_vs_imputed_distribution.csv"))

set.seed(20260804)
density_rows <- list()
row_index <- 1L
for (variable in imputed_variables) {
  observed <- imputation$data[[variable]]
  observed <- observed[!is.na(observed)]
  density_rows[[row_index]] <- tibble(
    variable = variable,
    source = "Observed",
    value = sample(observed, min(length(observed), 30000L), replace = FALSE)
  )
  row_index <- row_index + 1L
  imputed_values <- unlist(imputation$imp[[variable]], use.names = FALSE)
  density_rows[[row_index]] <- tibble(
    variable = variable,
    source = "Imputed",
    value = sample(imputed_values, min(length(imputed_values), 30000L), replace = FALSE)
  )
  row_index <- row_index + 1L
}
density_data <- bind_rows(density_rows)
density_plot_data <- density_data %>%
  mutate(variable_label = recode(variable, !!!variable_labels))
density_plot <- ggplot(density_plot_data, aes(value, colour = source, linetype = source)) +
  geom_density(linewidth = 0.6, adjust = 1.1) +
  facet_wrap(~ variable_label, scales = "free", ncol = 2) +
  scale_colour_manual(values = c("Observed" = "#111111", "Imputed" = "#777777")) +
  labs(x = NULL, y = "Density", colour = NULL, linetype = NULL) +
  theme_classic(base_size = 10) +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold"), legend.position = "bottom")
ggsave(
  file.path(figures_dir, "Supplementary_Figure_S3_observed_vs_imputed_distributions.tif"),
  density_plot,
  width = 7.2,
  height = 5.8,
  units = "in",
  dpi = 600,
  compression = "lzw",
  bg = "white"
)

print(final_convergence)
message_time("Primary landmark MICE diagnostics completed")
