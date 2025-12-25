rm(list = ls())

library(strong)
library(tidyverse)

tbl_num <- 0
cohort_attrib <- "original"

load("./load/label_df.Rdata")

add_tbl_obj(tbl_num)

tbl <- ok_icm_tbl1(
  vars, col_name_gp, df,
  var_header_name = header_name,
  var_header_loc = header_loc
)

tbl_cap <- paste0(
  "Supplementary Table ", update_obj(tbl_num),
  ". Basic demographic characteristics of the ", cohort_attrib, " cohort"
)

assign(tbl_name, tbl, envir = .GlobalEnv)
assign(tbl_cap_name, tbl_cap, envir = .GlobalEnv)

save_tbl_obj(tbl_name, tbl_cap_name)
