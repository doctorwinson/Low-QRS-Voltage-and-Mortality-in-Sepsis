rm(list = ls())

library(strong)
library(tidyverse)
library(flextable)
library(officer)

source("./script/cus_funs.R")

# Load table objects
load_rdata_obj("./load/tbl/")

# Load table captions
tbl_cap <- ls() %>%
  str_subset("^tbl_cap_\\d+$") %>%
  str_sort(numeric = TRUE)

tbl_name <- rdata_name

# Supplementary table list
sup_tbl_ls <- data.frame(
  tbl = tbl_name,
  cap = sapply(tbl_cap, function(x) get(x))
) %>%
  flextable() %>%
  width(j = c(1, 2), width = c(1, 12)) %>%
  set_header_labels(tbl = "Table", cap = "Caption")

sup_tbl_ls

save_as_html(
  sup_tbl_ls,
  path = "./docx/sup_tbl_ls.html",
  title = "Supplementary Table List"
)

# Table 1 before vs after matching
main_tbl1 <- ok_bf_af_tbl1(tbl_3, tbl_4)
main_tbl1
main_tbl_cap1 <- "Table 1. Baseline characteristics before and after propensity score matching of two cohorts."

# Merge main model result tables
tbl_res <- ls() %>%
  str_subset("^tbl_res_\\d+$") %>%
  str_sort(numeric = TRUE)

mod_res_tbl(tbl_res) %>%
  assign("main_tbl2", ., .GlobalEnv)

main_tbl2
main_tbl_cap2 <- "Table 2. Primary outcome with different models for cohort."

main_tbl_name <- paste0("main_tbl", 1:2)
main_tbl_cap  <- paste0("main_tbl_cap", 1:2)

# Write main tables to a single DOCX
write_tbl_to_docx(
  write_path = "./docx/main_tbl.docx",
  tbl_name = main_tbl_name,
  cap_name = main_tbl_cap
)

# Write supplementary tables to a single DOCX
write_tbl_to_docx(
  write_path = "./docx/sup_tbl.docx",
  tbl_name = tbl_name,
  cap_name = tbl_cap
)
