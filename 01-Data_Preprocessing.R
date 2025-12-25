rm(list = ls())

library(tidyverse)
library(lubridate)
library(strong)

df <- read.csv("./csv/20251215.csv")

df <- df %>%
  mutate(
    admittime = parse_date_time(admittime, orders = c("ymd HMS", "ymd HMSOS"), tz = "UTC"),
    intime    = parse_date_time(intime,    orders = c("ymd HMS", "ymd HMSOS"), tz = "UTC"),
    ecg_time  = parse_date_time(ecg_time,  orders = c("ymd HMS", "ymd HMSOS"), tz = "UTC")
  ) %>%
  filter(ecg_time >= intime - hours(24), ecg_time <= intime + hours(24)) %>%
  mutate(group = if_else(!is.na(low_qrs) & low_qrs == 1, 1L, 0L)) %>%
  arrange(stay_id, desc(group), ecg_time) %>%
  distinct(stay_id, .keep_all = TRUE) %>%
  filter(
    crtr_sepsis3 == 1,
    day28_los > 0,
    age > 18,
    icd_pregnancy == 0,
    icd_malignancy == 0,
    icu_hadm_order == 1
  ) %>%
  mutate(Group = group)

write.csv(df, "./data/ecg.csv", row.names = FALSE)
save(df, file = "./data/df.Rdata")

rm(list = ls())

library(tidyverse)
library(strong)

load("./data/df.Rdata")

ok_unadj_km(time = "day28_los", outcome = "day28_outcome", col_name_gp = "Group", df = df)

col_name_gp <- "Group"
gp_name <- c("NQRSV", "LQRSV")
header_name <- c(
  "Interventions (boolean for 1st 24 h)",
  "Comorbidities (boolean)",
  "Vital signs (1st 24 h)",
  "Laboratory tests (1st 24 h)"
)
header_loc <- c(7, 12, 21, 26)

df <- df %>%
  select(
    stay_id, intime, outtime, dischtime,
    day28_los, day28_outcome, icu_los, icu_outcome, hos_los, hos_outcome,
    age, gender, weight,
    sapsii, sofa, charlson,
    itvtn_24h_vent_tag, itvtn_24h_rrt_tag, drug_24h_vaso_tag, drug_24h_sedative_tag,
    icd_hf, icd_afib, icd_renal, icd_liver, icd_copd, icd_lc, icd_stroke, icd_cad,
    vs_24h_map_first, vs_24h_heart_rate_first, vs_24h_temp_first, vs_24h_resp_rate_first,
    lab_24h_wbc_first, lab_24h_rbc_first, lab_24h_hemoglobin_first, lab_24h_hct_first, lab_24h_platelet_first,
    lab_24h_na_first, lab_24h_k_first, lab_24h_hco3_first, lab_24h_cl_first,
    lab_24h_bun_first, lab_24h_lactate_first, lab_24h_creatinine_first,
    lab_24h_inr_first, lab_24h_pt_first,
    lab_24h_ph_first,
    Group
  ) %>%
  mutate(day28_los = day28_los / 24) %>%
  mutate(
    Group = factor(Group, levels = c(0, 1), labels = gp_name),
    gender = factor(gender, levels = c("F", "M"), labels = c("Female", "Male")),
    itvtn_24h_vent_tag = factor(itvtn_24h_vent_tag, levels = c(1, 0), labels = c("YES", "NO")),
    itvtn_24h_rrt_tag  = factor(itvtn_24h_rrt_tag,  levels = c(1, 0), labels = c("YES", "NO")),
    drug_24h_vaso_tag  = factor(drug_24h_vaso_tag,  levels = c(1, 0), labels = c("YES", "NO")),
    drug_24h_sedative_tag = factor(drug_24h_sedative_tag, levels = c(1, 0), labels = c("YES", "NO")),
    icd_hf = factor(icd_hf, levels = c(1, 0), labels = c("YES", "NO")),
    icd_afib = factor(icd_afib, levels = c(1, 0), labels = c("YES", "NO")),
    icd_renal = factor(icd_renal, levels = c(1, 0), labels = c("YES", "NO")),
    icd_liver = factor(icd_liver, levels = c(1, 0), labels = c("YES", "NO")),
    icd_copd = factor(icd_copd, levels = c(1, 0), labels = c("YES", "NO")),
    icd_cad = factor(icd_cad, levels = c(1, 0), labels = c("YES", "NO")),
    icd_stroke = factor(icd_stroke, levels = c(1, 0), labels = c("YES", "NO")),
    icd_lc = factor(icd_lc, levels = c(1, 0), labels = c("YES", "NO"))
  ) %>%
  rename(
    Age = age,
    Gender = gender,
    Weight = weight,
    SAPSii = sapsii,
    SOFA_score = sofa,
    charlson = charlson,
    Mechanical_ventilation_use = itvtn_24h_vent_tag,
    Renal_Replacement_Therapy_use = itvtn_24h_rrt_tag,
    Vasopressor_use = drug_24h_vaso_tag,
    Sedative_use = drug_24h_sedative_tag,
    Heart_Failure = icd_hf,
    Atrial_fibrillation = icd_afib,
    Renal = icd_renal,
    Liver = icd_liver,
    COPD = icd_copd,
    CAD = icd_cad,
    Stroke = icd_stroke,
    liver_cirrhosis = icd_lc,
    MAP = vs_24h_map_first,
    Heart_rate = vs_24h_heart_rate_first,
    Temperature = vs_24h_temp_first,
    RR = vs_24h_resp_rate_first,
    WBC = lab_24h_wbc_first,
    RBC = lab_24h_rbc_first,
    Hemoglobin = lab_24h_hemoglobin_first,
    Hct = lab_24h_hct_first,
    Platelet = lab_24h_platelet_first,
    Sodium = lab_24h_na_first,
    Potassium = lab_24h_k_first,
    Bicarbonate = lab_24h_hco3_first,
    Chloride = lab_24h_cl_first,
    BUN = lab_24h_bun_first,
    Lactate = lab_24h_lactate_first,
    Creatinine = lab_24h_creatinine_first,
    INR = lab_24h_inr_first,
    PT = lab_24h_pt_first,
    pH = lab_24h_ph_first
  )

covars <- names(df)[11:47]
vars <- covars
cat_covar <- covars[c(2, 7:18)]
num_covar <- covars[c(1, 3, 4:6, 19:37)]

ck_num_aligned(df, covars, num_covar)
ck_factor_aligned(df, covars, cat_covar)


# Multiple imputation (MI) for missing covariate values using `strong::ok_mi`, 
# which is implemented based on the `mice` R package.
df <- ok_mi(df, covars, seed = 42)

save(
  df, vars, covars, cat_covar, num_covar, col_name_gp, gp_name, header_name, header_loc,
  file = "./load/label_df.Rdata"
)

write.csv(df, "./data/data.csv", row.names = FALSE)
