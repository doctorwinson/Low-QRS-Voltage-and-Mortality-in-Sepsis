rm(list = ls())

library(strong)
library(tidyverse)

load("./load/matched_cohort.Rdata", .GlobalEnv)

df_psm <- psm_res$af_matching_df

table(df_psm$Group)
table(df_psm$Group, df_psm$day28_outcome)

# Kaplan–Meier curve after propensity score matching (PSM)
ok_unadj_km(
  time = "day28_los",
  outcome = "day28_outcome",
  col_name_gp = "Group",
  risk.table = TRUE,
  df = df_psm
)

df_iptw <- psm_res$bf_matching_df

table(df_iptw$Group, df_iptw$day28_outcome)

# IPTW-adjusted Kaplan–Meier curve
ok_adj_km(
  time = "day28_los",
  outcome = "day28_outcome",
  col_name_gp = "Group",
  df = df_iptw,
  risk.table = TRUE,
  conf.int = TRUE,
  psm_wt = "iptw"
)

