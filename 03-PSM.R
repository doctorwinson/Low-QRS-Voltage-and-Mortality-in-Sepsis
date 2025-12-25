rm(list = ls())
options(stringsAsFactors = FALSE)

library(strong)
library(tidyverse)

# Load custom functions
source("./script/cus_funs.R")

# Data preparation
load("./load/label_df.Rdata")

df <- as.data.frame(df)

psm_cohorts(
  # Shared cohort attribute
  cohort_attrib = "",
  # Cohort name (cohort 1)
  cohort_num = "",
  # Starting table index
  tbl_num = 1,
  # Starting plot index
  plot_num = 0
)
