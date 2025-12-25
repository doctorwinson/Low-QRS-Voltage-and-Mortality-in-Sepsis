## ## Multivariable Cox subgroup analysis (jstable) + forest plot ## 
rm(list = ls())

library(dplyr)
library(tidyr)
library(purrr)
library(survival)
library(jstable)
library(forestploter)
library(grid)
library(tibble)

load("./load/matched_cohort.Rdata", .GlobalEnv)

df3 <- psm_res$af_matching_df

stopifnot(all(c("Group", "day28_los", "day28_outcome") %in% names(df3)))

to_factor01 <- function(x, labels = c("No", "Yes")) {
  if (is.factor(x)) return(x)
  if (is.character(x)) return(factor(x))
  if (is.logical(x)) return(factor(ifelse(x, labels[2], labels[1]), levels = labels))
  if (is.numeric(x)) {
    ux <- sort(unique(x[!is.na(x)]))
    if (length(ux) <= 2 && all(ux %in% c(0, 1))) {
      return(factor(ifelse(x == 1, labels[2], labels[1]), levels = labels))
    }
  }
  x
}

df3 <- df3 %>%
  mutate(
    Group = factor(Group),
    day28_outcome = as.integer(day28_outcome),
    day28_los = as.numeric(day28_los),
    Age = as.numeric(Age),
    Age_cat = factor(ifelse(Age < 65, "<65year", "≥65year"), levels = c("<65year", "≥65year")),
    Gender = factor(Gender)
  )

var_subgroups <- c(
  "Age_cat", "Gender",
  "Mechanical_ventilation_use",
  "Renal_Replacement_Therapy_use",
  "Vasopressor_use",
  "Sedative_use",
  "Heart_Failure",
  "Atrial_fibrillation",
  "Renal", "Liver", "COPD", "liver_cirrhosis", "Stroke", "CAD"
)

missing_subg <- setdiff(var_subgroups, names(df3))
if (length(missing_subg) > 0) stop("Missing subgroup variables: ", paste(missing_subg, collapse = ", "))

df3 <- df3 %>% mutate(across(all_of(var_subgroups), to_factor01))

var_cov <- c(
  "Age", "Gender", "weight",
  "sofa",
  "charlson",
  "vs_24h_map_first",
  "lab_24h_lactate_first",
  "lab_24h_creatinine_first",
  "lab_24h_k_first",
  "lab_24h_hemoglobin_first",
  "lab_24h_inr_first"
)

var_cov <- intersect(var_cov, names(df3))
if (length(var_cov) == 0) stop("No covariates found in df3. Check variable names.")

cluster_var <- if ("subclass" %in% names(df3)) "subclass" else NULL

res <- TableSubgroupMultiCox(
  formula        = Surv(day28_los, day28_outcome) ~ Group,
  var_subgroups  = var_subgroups,
  var_cov        = var_cov,
  data           = df3,
  cluster        = cluster_var,
  decimal.hr     = 2,
  decimal.pvalue = 3
)

plot_df <- as_tibble(res) %>%
  mutate(
    Variable = if_else(Variable == "NA", "", Variable),
    `Point Estimate` = suppressWarnings(as.numeric(`Point Estimate`)),
    Lower = suppressWarnings(as.numeric(Lower)),
    Upper = suppressWarnings(as.numeric(Upper)),
    `HR (95% CI)` = if_else(
      is.na(`Point Estimate`), "",
      sprintf("%.2f (%.2f–%.2f)", `Point Estimate`, Lower, Upper)
    ),
    ci = paste(rep(" ", 60), collapse = ""),
    Count   = if_else(is.na(Count), " ", as.character(Count)),
    Percent = if_else(is.na(Percent), " ", as.character(Percent)),
    `P value` = if_else(is.na(`P value`), " ", as.character(`P value`)),
    `P for interaction` = if_else(is.na(`P for interaction`), " ", as.character(`P for interaction`))
  ) %>%
  filter(!(Variable == "" & is.na(`Point Estimate`) & is.na(Lower) & is.na(Upper)))

display_df <- plot_df %>%
  select(Variable, Count, Percent, ci, `HR (95% CI)`, `P value`, `P for interaction`)

xrange <- range(c(plot_df$Lower, plot_df$Upper), na.rm = TRUE)
xlim_use <- c(max(0.01, xrange[1] * 0.8), xrange[2] * 1.2)

p <- forest(
  data      = display_df,
  est       = plot_df$`Point Estimate`,
  lower     = plot_df$Lower,
  upper     = plot_df$Upper,
  ci_column = which(names(display_df) == "ci"),
  ref_line  = 1,
  x_trans   = "log",
  xlim      = xlim_use,
  ticks_at  = c(0.25, 0.5, 1, 2, 4),
  theme     = forest_theme(
    ci_lwd = 1.2,
    ci_Theight = unit(0.15, "cm")
  )
)

print(res)
plot(p)
