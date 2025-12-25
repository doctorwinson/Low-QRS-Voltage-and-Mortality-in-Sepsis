# Clear workspace
rm(list = ls())
options(stringsAsFactors = FALSE, scipen = 999)

library(strong)
library(tidyverse)

# Starting indices for tables and plots
tbl_num <- 4
plot_num <- 8

source("./script/cus_funs.R")

# Cohort mortality analysis parameters
pars <- list(
  mort_fml = c("day28"),
  mort_title = c("28-day"),
  cohort_attrib = "",
  psm_wt = "iptw",
  tbl_num = tbl_num,
  plot_num = plot_num
)

for (i in seq_along(pars$mort_title)) {
  starttime <- Sys.time()
  cat(sprintf("\n=== Total %d tasks; running task %d ===\n", length(pars$mort_title), i))
  cat(sprintf("\n=== Analyzing %s mortality models for the cohort ===\n", pars$mort_title[i]))
  cohort_mod(
    pars$mort_fml[i],
    pars$mort_title[i],
    pars$cohort_attrib,
    pars$psm_wt,
    pars$tbl_num,
    pars$plot_num
  )
  endtime <- Sys.time()
  print(endtime - starttime)
  cat("\n\n\n")
  gc()
}
