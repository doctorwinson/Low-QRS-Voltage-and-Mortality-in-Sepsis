/*
Purpose:
  Build the study base population directly from official MIMIC-IV v2.2 and
  MIMIC-IV derived tables. This replaces the earlier private dependency on
  a non-public local cohort table.

Output:
  lowvoltage.official_analysis_population

Run before:
  02_extract_sepsis_ecg_cohort.sql
  03_extract_post_icu_24h_cohort.sql
  04_cohort_flow_diagnostics.sql
  05_ecg_timing_and_alternative_exposures.sql

Implementation notes:
  - Sepsis-3 is read from mimiciv_derived.sepsis3.
  - Mortality follow-up starts at ICU intime and is truncated at 28 days.
  - First-24-hour covariates are the first non-missing value for each variable
    within 24 hours after ICU admission.
  - ICD flags are rebuilt explicitly from mimiciv_hosp.diagnoses_icd.
*/

CREATE SCHEMA IF NOT EXISTS lowvoltage;

DROP TABLE IF EXISTS lowvoltage.official_analysis_population;
DROP TABLE IF EXISTS tmp_lv_icu_base;
DROP TABLE IF EXISTS tmp_lv_dx_flags;
DROP TABLE IF EXISTS tmp_lv_first_vitals;
DROP TABLE IF EXISTS tmp_lv_first_cbc;
DROP TABLE IF EXISTS tmp_lv_first_chem;
DROP TABLE IF EXISTS tmp_lv_first_coag;
DROP TABLE IF EXISTS tmp_lv_first_bg;
DROP TABLE IF EXISTS tmp_lv_interventions;
DROP TABLE IF EXISTS tmp_lv_sedatives;

CREATE TEMP TABLE tmp_lv_icu_base AS
SELECT
    ie.subject_id,
    ie.hadm_id,
    ie.stay_id,
    ie.first_careunit,
    ad.admittime,
    ad.admission_type,
    ad.admission_location,
    ie.intime,
    ie.outtime,
    ad.deathtime AS in_hos_dod,
    ad.dischtime,
    pt.dod AS out_hos_dod,
    /* patients.dod is date-level. Assign DOD-only deaths to the end of the
       recorded day so a same-day date is not treated as preceding ICU entry. */
    COALESCE(ad.deathtime, pt.dod::timestamp + INTERVAL '23 hours 59 minutes 59 seconds') AS final_dod,
    CASE
        WHEN COALESCE(ad.deathtime, pt.dod::timestamp + INTERVAL '23 hours 59 minutes 59 seconds') IS NOT NULL
        THEN EXTRACT(EPOCH FROM (COALESCE(ad.deathtime, pt.dod::timestamp + INTERVAL '23 hours 59 minutes 59 seconds') - ie.intime)) / 3600.0
        ELSE NULL
    END AS surv_h,
    CASE
        WHEN ad.deathtime IS NOT NULL
         AND ad.deathtime > ie.intime
         AND ad.deathtime <= ie.outtime THEN 1
        ELSE 0
    END AS icu_outcome,
    EXTRACT(EPOCH FROM (ie.outtime - ie.intime)) / 3600.0 AS icu_los,
    ad.hospital_expire_flag::integer AS hos_outcome,
    EXTRACT(EPOCH FROM (ad.dischtime - ad.admittime)) / 3600.0 AS hos_los,
    /* Exact in-hospital death times are compared as timestamps. When only
       patients.dod is available, mortality is classified by calendar date:
       same-day DOD remains eligible and deaths on the 28th calendar day are
       included. This avoids imposing an unsupported midnight or end-of-day
       time on the primary binary outcome. */
    CASE
        WHEN ad.deathtime IS NOT NULL AND ad.deathtime <= ie.intime THEN NULL::integer
        WHEN ad.deathtime IS NOT NULL AND ad.deathtime <= ie.intime + INTERVAL '28 days' THEN 1
        WHEN ad.deathtime IS NOT NULL THEN 0
        WHEN pt.dod IS NOT NULL AND pt.dod < ie.intime::date THEN NULL::integer
        WHEN pt.dod IS NOT NULL AND pt.dod <= (ie.intime + INTERVAL '28 days')::date THEN 1
        ELSE 0
    END AS day28_outcome,
    CASE
        WHEN ad.deathtime IS NOT NULL AND ad.deathtime <= ie.intime THEN NULL::numeric
        WHEN ad.deathtime IS NOT NULL AND ad.deathtime <= ie.intime + INTERVAL '28 days'
        THEN EXTRACT(EPOCH FROM (ad.deathtime - ie.intime)) / 3600.0
        WHEN ad.deathtime IS NOT NULL THEN 28 * 24
        WHEN pt.dod IS NOT NULL AND pt.dod < ie.intime::date THEN NULL::numeric
        WHEN pt.dod IS NOT NULL AND pt.dod <= (ie.intime + INTERVAL '28 days')::date
        THEN LEAST(
            EXTRACT(EPOCH FROM (pt.dod::timestamp + INTERVAL '23 hours 59 minutes 59 seconds' - ie.intime)) / 3600.0,
            28 * 24
        )
        ELSE 28 * 24
    END AS day28_los,
    COALESCE(id.admission_age, ag.age) AS age,
    ROW_NUMBER() OVER (PARTITION BY ie.subject_id ORDER BY ie.intime, ie.stay_id) AS icu_subject_order,
    ROW_NUMBER() OVER (PARTITION BY ie.hadm_id ORDER BY ie.intime, ie.stay_id) AS icu_hadm_order,
    ad.race,
    pt.gender,
    pt.anchor_year,
    pt.anchor_year_group,
    CASE
        WHEN ad.deathtime IS NOT NULL AND ad.deathtime <= ie.intime THEN 1
        WHEN ad.deathtime IS NULL AND pt.dod IS NOT NULL AND pt.dod < ie.intime::date THEN 1
        ELSE 0
    END AS surv_h_error_tag
FROM mimiciv_icu.icustays AS ie
INNER JOIN mimiciv_hosp.admissions AS ad
    ON ie.subject_id = ad.subject_id
   AND ie.hadm_id = ad.hadm_id
INNER JOIN mimiciv_hosp.patients AS pt
    ON ie.subject_id = pt.subject_id
LEFT JOIN mimiciv_derived.icustay_detail AS id
    ON ie.stay_id = id.stay_id
LEFT JOIN mimiciv_derived.age AS ag
    ON ie.subject_id = ag.subject_id
   AND ie.hadm_id = ag.hadm_id;

CREATE UNIQUE INDEX tmp_lv_icu_base_stay_idx ON tmp_lv_icu_base (stay_id);
CREATE INDEX tmp_lv_icu_base_subject_hadm_idx ON tmp_lv_icu_base (subject_id, hadm_id);
CREATE INDEX tmp_lv_icu_base_subject_time_idx ON tmp_lv_icu_base (subject_id, intime);
ANALYZE tmp_lv_icu_base;

CREATE TEMP TABLE tmp_lv_dx_flags AS
SELECT
    d.subject_id,
    d.hadm_id,
    1 AS icd_avail_tag,
    MAX(CASE
        WHEN (d.icd_version = 9 AND d.icd_code ~ '^(63[0-9]|64[0-9]|65[0-9]|66[0-9]|67[0-9]|V22|V23|V24|V27)')
          OR (d.icd_version = 10 AND (d.icd_code LIKE 'O%' OR d.icd_code ~ '^(Z33|Z34|Z36|Z37|Z39)'))
        THEN 1 ELSE 0
    END) AS icd_pregnancy,
    MAX(CASE
        /* Malignant neoplasms only. Do not classify benign or uncertain-behavior
           neoplasms (ICD-9 235-239; ICD-10 D3A/D37-D48) as malignancy. */
        WHEN (d.icd_version = 9 AND (
                  d.icd_code ~ '^(14[0-9]|15[0-9]|16[0-9]|17[0-9]|18[0-9]|19[0-9]|20[0-8])'
                  OR d.icd_code ~ '^209[0-2]'
                  OR d.icd_code ~ '^209(3[0-6]|7[0-9])'
              ))
          OR (d.icd_version = 10 AND d.icd_code LIKE 'C%')
        THEN 1 ELSE 0
    END) AS icd_malignancy,
    MAX(CASE
        WHEN (d.icd_version = 9 AND d.icd_code IN ('42731', '42732'))
          OR (d.icd_version = 10 AND d.icd_code LIKE 'I48%')
        THEN 1 ELSE 0
    END) AS icd_afib,
    MAX(CASE
        WHEN (d.icd_version = 9 AND (d.icd_code ~ '^(410|411|412|413|414)' OR d.icd_code IN ('V4581', 'V4582')))
          OR (d.icd_version = 10 AND (d.icd_code ~ '^(I2[0-5])' OR d.icd_code IN ('Z951', 'Z955')))
        THEN 1 ELSE 0
    END) AS icd_cad,
    MAX(CASE
        WHEN (d.icd_version = 9 AND d.icd_code ~ '^(5712|5715|5716)')
          OR (d.icd_version = 10 AND d.icd_code ~ '^(K703|K717|K743|K744|K745|K746)')
        THEN 1 ELSE 0
    END) AS icd_lc
FROM mimiciv_hosp.diagnoses_icd AS d
INNER JOIN (SELECT DISTINCT subject_id, hadm_id FROM tmp_lv_icu_base) AS b
    ON d.subject_id = b.subject_id
   AND d.hadm_id = b.hadm_id
GROUP BY d.subject_id, d.hadm_id;

CREATE UNIQUE INDEX tmp_lv_dx_flags_idx ON tmp_lv_dx_flags (subject_id, hadm_id);
ANALYZE tmp_lv_dx_flags;

CREATE TEMP TABLE tmp_lv_first_vitals AS
SELECT
    b.stay_id,
    COUNT(v.mbp) AS vs_24h_map_tag,
    (ARRAY_AGG(v.mbp ORDER BY v.charttime) FILTER (WHERE v.mbp IS NOT NULL))[1] AS vs_24h_map_first,
    MIN(v.mbp) AS vs_24h_map_min,
    MAX(v.mbp) AS vs_24h_map_max,
    COUNT(v.heart_rate) AS vs_24h_heart_rate_tag,
    (ARRAY_AGG(v.heart_rate ORDER BY v.charttime) FILTER (WHERE v.heart_rate IS NOT NULL))[1] AS vs_24h_heart_rate_first,
    MIN(v.heart_rate) AS vs_24h_heart_rate_min,
    MAX(v.heart_rate) AS vs_24h_heart_rate_max,
    COUNT(v.temperature) AS vs_24h_temp_tag,
    (ARRAY_AGG(v.temperature ORDER BY v.charttime) FILTER (WHERE v.temperature IS NOT NULL))[1] AS vs_24h_temp_first,
    MIN(v.temperature) AS vs_24h_temp_min,
    MAX(v.temperature) AS vs_24h_temp_max,
    COUNT(v.resp_rate) AS vs_24h_resp_rate_tag,
    (ARRAY_AGG(v.resp_rate ORDER BY v.charttime) FILTER (WHERE v.resp_rate IS NOT NULL))[1] AS vs_24h_resp_rate_first,
    MIN(v.resp_rate) AS vs_24h_resp_rate_min,
    MAX(v.resp_rate) AS vs_24h_resp_rate_max,
    COUNT(v.spo2) AS vs_24h_spo2_tag,
    (ARRAY_AGG(v.spo2 ORDER BY v.charttime) FILTER (WHERE v.spo2 IS NOT NULL))[1] AS vs_24h_spo2_first,
    MIN(v.spo2) AS vs_24h_spo2_min,
    MAX(v.spo2) AS vs_24h_spo2_max
FROM tmp_lv_icu_base AS b
LEFT JOIN mimiciv_derived.vitalsign AS v
    ON b.stay_id = v.stay_id
   AND v.charttime >= b.intime
   AND v.charttime <= b.intime + INTERVAL '24 hours'
GROUP BY b.stay_id;

CREATE UNIQUE INDEX tmp_lv_first_vitals_idx ON tmp_lv_first_vitals (stay_id);
ANALYZE tmp_lv_first_vitals;

CREATE TEMP TABLE tmp_lv_first_cbc AS
SELECT
    b.stay_id,
    COUNT(cbc.wbc) AS lab_24h_wbc_tag,
    (ARRAY_AGG(cbc.wbc ORDER BY cbc.charttime) FILTER (WHERE cbc.wbc IS NOT NULL))[1] AS lab_24h_wbc_first,
    MIN(cbc.wbc) AS lab_24h_wbc_min,
    MAX(cbc.wbc) AS lab_24h_wbc_max,
    COUNT(cbc.rbc) AS lab_24h_rbc_tag,
    (ARRAY_AGG(cbc.rbc ORDER BY cbc.charttime) FILTER (WHERE cbc.rbc IS NOT NULL))[1] AS lab_24h_rbc_first,
    MIN(cbc.rbc) AS lab_24h_rbc_min,
    MAX(cbc.rbc) AS lab_24h_rbc_max,
    COUNT(cbc.hemoglobin) AS lab_24h_hemoglobin_tag,
    (ARRAY_AGG(cbc.hemoglobin ORDER BY cbc.charttime) FILTER (WHERE cbc.hemoglobin IS NOT NULL))[1] AS lab_24h_hemoglobin_first,
    MIN(cbc.hemoglobin) AS lab_24h_hemoglobin_min,
    MAX(cbc.hemoglobin) AS lab_24h_hemoglobin_max,
    COUNT(cbc.hematocrit) AS lab_24h_hct_tag,
    (ARRAY_AGG(cbc.hematocrit ORDER BY cbc.charttime) FILTER (WHERE cbc.hematocrit IS NOT NULL))[1] AS lab_24h_hct_first,
    MIN(cbc.hematocrit) AS lab_24h_hct_min,
    MAX(cbc.hematocrit) AS lab_24h_hct_max,
    COUNT(cbc.platelet) AS lab_24h_platelet_tag,
    (ARRAY_AGG(cbc.platelet ORDER BY cbc.charttime) FILTER (WHERE cbc.platelet IS NOT NULL))[1] AS lab_24h_platelet_first,
    MIN(cbc.platelet) AS lab_24h_platelet_min,
    MAX(cbc.platelet) AS lab_24h_platelet_max
FROM tmp_lv_icu_base AS b
LEFT JOIN mimiciv_derived.complete_blood_count AS cbc
    ON b.subject_id = cbc.subject_id
   AND b.hadm_id = cbc.hadm_id
   AND cbc.charttime >= b.intime
   AND cbc.charttime <= b.intime + INTERVAL '24 hours'
GROUP BY b.stay_id;

CREATE UNIQUE INDEX tmp_lv_first_cbc_idx ON tmp_lv_first_cbc (stay_id);
ANALYZE tmp_lv_first_cbc;

CREATE TEMP TABLE tmp_lv_first_chem AS
SELECT
    b.stay_id,
    COUNT(ch.sodium) AS lab_24h_na_tag,
    (ARRAY_AGG(ch.sodium ORDER BY ch.charttime) FILTER (WHERE ch.sodium IS NOT NULL))[1] AS lab_24h_na_first,
    MIN(ch.sodium) AS lab_24h_na_min,
    MAX(ch.sodium) AS lab_24h_na_max,
    COUNT(ch.potassium) AS lab_24h_k_tag,
    (ARRAY_AGG(ch.potassium ORDER BY ch.charttime) FILTER (WHERE ch.potassium IS NOT NULL))[1] AS lab_24h_k_first,
    MIN(ch.potassium) AS lab_24h_k_min,
    MAX(ch.potassium) AS lab_24h_k_max,
    COUNT(ch.bicarbonate) AS lab_24h_hco3_chem_tag,
    (ARRAY_AGG(ch.bicarbonate ORDER BY ch.charttime) FILTER (WHERE ch.bicarbonate IS NOT NULL))[1] AS lab_24h_hco3_chem_first,
    MIN(ch.bicarbonate) AS lab_24h_hco3_chem_min,
    MAX(ch.bicarbonate) AS lab_24h_hco3_chem_max,
    COUNT(ch.chloride) AS lab_24h_cl_tag,
    (ARRAY_AGG(ch.chloride ORDER BY ch.charttime) FILTER (WHERE ch.chloride IS NOT NULL))[1] AS lab_24h_cl_first,
    MIN(ch.chloride) AS lab_24h_cl_min,
    MAX(ch.chloride) AS lab_24h_cl_max,
    COUNT(ch.bun) AS lab_24h_bun_tag,
    (ARRAY_AGG(ch.bun ORDER BY ch.charttime) FILTER (WHERE ch.bun IS NOT NULL))[1] AS lab_24h_bun_first,
    MIN(ch.bun) AS lab_24h_bun_min,
    MAX(ch.bun) AS lab_24h_bun_max,
    COUNT(ch.creatinine) AS lab_24h_creatinine_tag,
    (ARRAY_AGG(ch.creatinine ORDER BY ch.charttime) FILTER (WHERE ch.creatinine IS NOT NULL))[1] AS lab_24h_creatinine_first,
    MIN(ch.creatinine) AS lab_24h_creatinine_min,
    MAX(ch.creatinine) AS lab_24h_creatinine_max
FROM tmp_lv_icu_base AS b
LEFT JOIN mimiciv_derived.chemistry AS ch
    ON b.subject_id = ch.subject_id
   AND b.hadm_id = ch.hadm_id
   AND ch.charttime >= b.intime
   AND ch.charttime <= b.intime + INTERVAL '24 hours'
GROUP BY b.stay_id;

CREATE UNIQUE INDEX tmp_lv_first_chem_idx ON tmp_lv_first_chem (stay_id);
ANALYZE tmp_lv_first_chem;

CREATE TEMP TABLE tmp_lv_first_coag AS
SELECT
    b.stay_id,
    COUNT(co.inr) AS lab_24h_inr_tag,
    (ARRAY_AGG(co.inr ORDER BY co.charttime) FILTER (WHERE co.inr IS NOT NULL))[1] AS lab_24h_inr_first,
    MIN(co.inr) AS lab_24h_inr_min,
    MAX(co.inr) AS lab_24h_inr_max,
    COUNT(co.pt) AS lab_24h_pt_tag,
    (ARRAY_AGG(co.pt ORDER BY co.charttime) FILTER (WHERE co.pt IS NOT NULL))[1] AS lab_24h_pt_first,
    MIN(co.pt) AS lab_24h_pt_min,
    MAX(co.pt) AS lab_24h_pt_max
FROM tmp_lv_icu_base AS b
LEFT JOIN mimiciv_derived.coagulation AS co
    ON b.subject_id = co.subject_id
   AND b.hadm_id = co.hadm_id
   AND co.charttime >= b.intime
   AND co.charttime <= b.intime + INTERVAL '24 hours'
GROUP BY b.stay_id;

CREATE UNIQUE INDEX tmp_lv_first_coag_idx ON tmp_lv_first_coag (stay_id);
ANALYZE tmp_lv_first_coag;

CREATE TEMP TABLE tmp_lv_first_bg AS
SELECT
    b.stay_id,
    COUNT(bg.lactate) AS lab_24h_lactate_tag,
    (ARRAY_AGG(bg.lactate ORDER BY bg.charttime) FILTER (WHERE bg.lactate IS NOT NULL))[1] AS lab_24h_lactate_first,
    MIN(bg.lactate) AS lab_24h_lactate_min,
    MAX(bg.lactate) AS lab_24h_lactate_max,
    COUNT(bg.ph) AS lab_24h_ph_tag,
    (ARRAY_AGG(bg.ph ORDER BY bg.charttime) FILTER (WHERE bg.ph IS NOT NULL))[1] AS lab_24h_ph_first,
    MIN(bg.ph) AS lab_24h_ph_min,
    MAX(bg.ph) AS lab_24h_ph_max,
    COUNT(bg.bicarbonate) AS lab_24h_hco3_bg_tag,
    (ARRAY_AGG(bg.bicarbonate ORDER BY bg.charttime) FILTER (WHERE bg.bicarbonate IS NOT NULL))[1] AS lab_24h_hco3_bg_first,
    MIN(bg.bicarbonate) AS lab_24h_hco3_bg_min,
    MAX(bg.bicarbonate) AS lab_24h_hco3_bg_max
FROM tmp_lv_icu_base AS b
LEFT JOIN mimiciv_derived.bg AS bg
    ON b.subject_id = bg.subject_id
   AND b.hadm_id = bg.hadm_id
   AND bg.charttime >= b.intime
   AND bg.charttime <= b.intime + INTERVAL '24 hours'
GROUP BY b.stay_id;

CREATE UNIQUE INDEX tmp_lv_first_bg_idx ON tmp_lv_first_bg (stay_id);
ANALYZE tmp_lv_first_bg;

CREATE TEMP TABLE tmp_lv_interventions AS
SELECT
    b.stay_id,
    MAX(CASE
        WHEN vent.ventilation_status IN ('InvasiveVent', 'NonInvasiveVent', 'Tracheostomy')
        THEN 1 ELSE 0
    END) AS itvtn_24h_vent_tag,
    MAX(CASE WHEN va.stay_id IS NOT NULL THEN 1 ELSE 0 END) AS drug_24h_vaso_tag
FROM tmp_lv_icu_base AS b
LEFT JOIN mimiciv_derived.ventilation AS vent
    ON b.stay_id = vent.stay_id
   AND vent.starttime < b.intime + INTERVAL '24 hours'
   AND COALESCE(vent.endtime, vent.starttime) >= b.intime
LEFT JOIN mimiciv_derived.vasoactive_agent AS va
    ON b.stay_id = va.stay_id
   AND va.starttime < b.intime + INTERVAL '24 hours'
   AND COALESCE(va.endtime, va.starttime) >= b.intime
GROUP BY b.stay_id;

CREATE UNIQUE INDEX tmp_lv_interventions_idx ON tmp_lv_interventions (stay_id);
ANALYZE tmp_lv_interventions;

CREATE TEMP TABLE tmp_lv_sedatives AS
SELECT
    b.stay_id,
    MAX(CASE WHEN lower(di.label) LIKE '%midazolam%' THEN 1 ELSE 0 END) AS drug_24h_mdz_tag,
    MAX(CASE WHEN lower(di.label) LIKE '%propofol%' THEN 1 ELSE 0 END) AS drug_24h_propofol_tag,
    MAX(CASE WHEN lower(di.label) LIKE '%dexmedetomidine%' THEN 1 ELSE 0 END) AS drug_24h_dex_tag,
    MAX(CASE WHEN lower(di.label) LIKE '%fentanyl%' THEN 1 ELSE 0 END) AS drug_24h_fen_tag,
    MAX(CASE WHEN lower(di.label) ~ '(midazolam|propofol|dexmedetomidine|fentanyl)' THEN 1 ELSE 0 END) AS drug_24h_sedative_tag
FROM tmp_lv_icu_base AS b
LEFT JOIN mimiciv_icu.inputevents AS inp
    ON b.stay_id = inp.stay_id
   AND inp.starttime < b.intime + INTERVAL '24 hours'
   AND COALESCE(inp.endtime, inp.starttime) >= b.intime
   AND COALESCE(inp.amount, inp.rate, inp.originalamount, inp.originalrate, 0) > 0
   AND COALESCE(inp.statusdescription, '') NOT ILIKE '%rewritten%'
LEFT JOIN mimiciv_icu.d_items AS di
    ON inp.itemid = di.itemid
   AND lower(di.label) ~ '(midazolam|propofol|dexmedetomidine|fentanyl)'
GROUP BY b.stay_id;

CREATE UNIQUE INDEX tmp_lv_sedatives_idx ON tmp_lv_sedatives (stay_id);
ANALYZE tmp_lv_sedatives;

CREATE TABLE lowvoltage.official_analysis_population AS
SELECT
    b.stay_id,
    b.hadm_id,
    b.subject_id,
    b.first_careunit,
    b.admittime,
    b.admission_type,
    b.admission_location,
    b.intime,
    b.outtime,
    b.in_hos_dod,
    b.dischtime,
    b.out_hos_dod,
    b.final_dod,
    b.surv_h,
    b.icu_outcome,
    b.icu_los,
    b.hos_outcome,
    b.hos_los,
    b.day28_outcome,
    b.day28_los,
    b.age,
    b.icu_subject_order,
    b.icu_hadm_order,
    b.race,
    COALESCE(w.weight, w.weight_admit) AS weight_raw,
    CASE
        WHEN COALESCE(w.weight, w.weight_admit) BETWEEN 20 AND 300
        THEN COALESCE(w.weight, w.weight_admit)
        ELSE NULL
    END AS weight,
    CASE
        WHEN COALESCE(w.weight, w.weight_admit) IS NOT NULL
         AND COALESCE(w.weight, w.weight_admit) NOT BETWEEN 20 AND 300
        THEN 1 ELSE 0
    END AS weight_out_of_range,
    NULL::numeric AS height,
    NULL::numeric AS bmi_omr,
    NULL::numeric AS bmi_cal,
    NULL::numeric AS bmi,
    b.gender,
    b.anchor_year,
    b.anchor_year_group,
    NULL::integer AS apsiii,
    sap.sapsii,
    NULL::integer AS oasis,
    NULL::integer AS lods,
    sofa.sofa,
    NULL::double precision AS gcs,
    ch.charlson_comorbidity_index AS charlson,
    NULL::numeric AS meld_original,
    NULL::double precision AS meld_2016,
    NULL::double precision AS meld3,
    COALESCE(dx.icd_avail_tag, 0) AS icd_avail_tag,
    b.surv_h_error_tag,
    COALESCE(se.sepsis3::integer, 0) AS crtr_sepsis3,
    se.suspected_infection_time AS crtr_sepsis3_suspected_infection_time,
    se.sofa_time AS crtr_sepsis3_sofa_time,
    LEAST(se.suspected_infection_time, se.sofa_time) AS crtr_sepsis3_time,
    EXTRACT(EPOCH FROM (se.suspected_infection_time - b.intime)) / 3600.0
        AS crtr_sepsis3_suspected_infection_hours_from_icu,
    COALESCE(iv.itvtn_24h_vent_tag, 0) AS itvtn_24h_vent_tag,
    COALESCE(rrt.dialysis_present, 0) AS itvtn_24h_rrt_tag,
    0 AS itvtn_24h_cabg_tag,
    0 AS itvtn_24h_pci_tag,
    0 AS itvtn_24h_iabp_tag,
    0 AS itvtn_24h_picco_tag,
    0 AS itvtn_24h_nicom_tag,
    COALESCE(sed.drug_24h_mdz_tag, 0) AS drug_24h_mdz_tag,
    COALESCE(sed.drug_24h_propofol_tag, 0) AS drug_24h_propofol_tag,
    COALESCE(sed.drug_24h_dex_tag, 0) AS drug_24h_dex_tag,
    COALESCE(sed.drug_24h_fen_tag, 0) AS drug_24h_fen_tag,
    COALESCE(sed.drug_24h_sedative_tag, 0) AS drug_24h_sedative_tag,
    0 AS drug_24h_albumin_tag,
    0 AS drug_24h_avp_tag,
    0 AS drug_24h_da_tag,
    0 AS drug_24h_dba_tag,
    0 AS drug_24h_epi_tag,
    0 AS drug_24h_mil_tag,
    0 AS drug_24h_ne_tag,
    0 AS drug_24h_pe_tag,
    0 AS drug_24h_aii_tag,
    COALESCE(iv.drug_24h_vaso_tag, 0) AS drug_24h_vaso_tag,
    0 AS drug_24h_insulin_tag,
    COALESCE(ch.congestive_heart_failure, 0) AS icd_hf,
    0 AS icd_hypertension,
    COALESCE(dx.icd_afib, 0) AS icd_afib,
    GREATEST(COALESCE(ch.diabetes_without_cc, 0), COALESCE(ch.diabetes_with_cc, 0)) AS icd_diabetes,
    COALESCE(ch.renal_disease, 0) AS icd_renal,
    GREATEST(COALESCE(ch.mild_liver_disease, 0), COALESCE(ch.severe_liver_disease, 0)) AS icd_liver,
    COALESCE(ch.chronic_pulmonary_disease, 0) AS icd_copd,
    COALESCE(dx.icd_cad, 0) AS icd_cad,
    COALESCE(ch.cerebrovascular_disease, 0) AS icd_stroke,
    COALESCE(dx.icd_malignancy, 0) AS icd_malignancy,
    0 AS icd_cs,
    0 AS icd_ss,
    COALESCE(dx.icd_pregnancy, 0) AS icd_pregnancy,
    0 AS icd_ccd,
    0 AS icd_pph,
    0 AS icd_ppcd,
    0 AS icd_ac,
    0 AS icd_hit,
    COALESCE(dx.icd_lc, 0) AS icd_lc,
    COALESCE(v.vs_24h_map_tag, 0) AS vs_24h_map_tag,
    v.vs_24h_map_first,
    v.vs_24h_map_min,
    v.vs_24h_map_max,
    COALESCE(v.vs_24h_heart_rate_tag, 0) AS vs_24h_heart_rate_tag,
    v.vs_24h_heart_rate_first,
    v.vs_24h_heart_rate_min,
    v.vs_24h_heart_rate_max,
    COALESCE(v.vs_24h_temp_tag, 0) AS vs_24h_temp_tag,
    v.vs_24h_temp_first,
    v.vs_24h_temp_min,
    v.vs_24h_temp_max,
    COALESCE(v.vs_24h_resp_rate_tag, 0) AS vs_24h_resp_rate_tag,
    v.vs_24h_resp_rate_first,
    v.vs_24h_resp_rate_min,
    v.vs_24h_resp_rate_max,
    COALESCE(v.vs_24h_spo2_tag, 0) AS vs_24h_spo2_tag,
    v.vs_24h_spo2_first,
    v.vs_24h_spo2_min,
    v.vs_24h_spo2_max,
    COALESCE(bg.lab_24h_ph_tag, 0) AS lab_24h_ph_tag,
    bg.lab_24h_ph_first,
    bg.lab_24h_ph_min,
    bg.lab_24h_ph_max,
    COALESCE(chem.lab_24h_hco3_chem_tag, 0) + COALESCE(bg.lab_24h_hco3_bg_tag, 0) AS lab_24h_hco3_tag,
    COALESCE(chem.lab_24h_hco3_chem_first, bg.lab_24h_hco3_bg_first) AS lab_24h_hco3_first,
    CASE
        WHEN chem.lab_24h_hco3_chem_min IS NULL THEN bg.lab_24h_hco3_bg_min
        WHEN bg.lab_24h_hco3_bg_min IS NULL THEN chem.lab_24h_hco3_chem_min
        ELSE LEAST(chem.lab_24h_hco3_chem_min, bg.lab_24h_hco3_bg_min)
    END AS lab_24h_hco3_min,
    CASE
        WHEN chem.lab_24h_hco3_chem_max IS NULL THEN bg.lab_24h_hco3_bg_max
        WHEN bg.lab_24h_hco3_bg_max IS NULL THEN chem.lab_24h_hco3_chem_max
        ELSE GREATEST(chem.lab_24h_hco3_chem_max, bg.lab_24h_hco3_bg_max)
    END AS lab_24h_hco3_max,
    COALESCE(bg.lab_24h_lactate_tag, 0) AS lab_24h_lactate_tag,
    bg.lab_24h_lactate_first,
    bg.lab_24h_lactate_min,
    bg.lab_24h_lactate_max,
    COALESCE(chem.lab_24h_creatinine_tag, 0) AS lab_24h_creatinine_tag,
    chem.lab_24h_creatinine_first,
    chem.lab_24h_creatinine_min,
    chem.lab_24h_creatinine_max,
    COALESCE(chem.lab_24h_bun_tag, 0) AS lab_24h_bun_tag,
    chem.lab_24h_bun_first,
    chem.lab_24h_bun_min,
    chem.lab_24h_bun_max,
    COALESCE(cbc.lab_24h_wbc_tag, 0) AS lab_24h_wbc_tag,
    cbc.lab_24h_wbc_first,
    cbc.lab_24h_wbc_min,
    cbc.lab_24h_wbc_max,
    COALESCE(cbc.lab_24h_rbc_tag, 0) AS lab_24h_rbc_tag,
    cbc.lab_24h_rbc_first,
    cbc.lab_24h_rbc_min,
    cbc.lab_24h_rbc_max,
    COALESCE(cbc.lab_24h_hct_tag, 0) AS lab_24h_hct_tag,
    cbc.lab_24h_hct_first,
    cbc.lab_24h_hct_min,
    cbc.lab_24h_hct_max,
    COALESCE(cbc.lab_24h_hemoglobin_tag, 0) AS lab_24h_hemoglobin_tag,
    cbc.lab_24h_hemoglobin_first,
    cbc.lab_24h_hemoglobin_min,
    cbc.lab_24h_hemoglobin_max,
    COALESCE(cbc.lab_24h_platelet_tag, 0) AS lab_24h_platelet_tag,
    cbc.lab_24h_platelet_first,
    cbc.lab_24h_platelet_min,
    cbc.lab_24h_platelet_max,
    COALESCE(coag.lab_24h_inr_tag, 0) AS lab_24h_inr_tag,
    coag.lab_24h_inr_first,
    coag.lab_24h_inr_min,
    coag.lab_24h_inr_max,
    COALESCE(coag.lab_24h_pt_tag, 0) AS lab_24h_pt_tag,
    coag.lab_24h_pt_first,
    coag.lab_24h_pt_min,
    coag.lab_24h_pt_max,
    COALESCE(chem.lab_24h_cl_tag, 0) AS lab_24h_cl_tag,
    chem.lab_24h_cl_first,
    chem.lab_24h_cl_min,
    chem.lab_24h_cl_max,
    COALESCE(chem.lab_24h_na_tag, 0) AS lab_24h_na_tag,
    chem.lab_24h_na_first,
    chem.lab_24h_na_min,
    chem.lab_24h_na_max,
    COALESCE(chem.lab_24h_k_tag, 0) AS lab_24h_k_tag,
    chem.lab_24h_k_first,
    chem.lab_24h_k_min,
    chem.lab_24h_k_max
FROM tmp_lv_icu_base AS b
LEFT JOIN mimiciv_derived.sepsis3 AS se
    ON b.stay_id = se.stay_id
LEFT JOIN mimiciv_derived.first_day_weight AS w
    ON b.stay_id = w.stay_id
LEFT JOIN mimiciv_derived.sapsii AS sap
    ON b.stay_id = sap.stay_id
LEFT JOIN mimiciv_derived.first_day_sofa AS sofa
    ON b.stay_id = sofa.stay_id
LEFT JOIN mimiciv_derived.charlson AS ch
    ON b.subject_id = ch.subject_id
   AND b.hadm_id = ch.hadm_id
LEFT JOIN tmp_lv_dx_flags AS dx
    ON b.subject_id = dx.subject_id
   AND b.hadm_id = dx.hadm_id
LEFT JOIN tmp_lv_first_vitals AS v
    ON b.stay_id = v.stay_id
LEFT JOIN tmp_lv_first_cbc AS cbc
    ON b.stay_id = cbc.stay_id
LEFT JOIN tmp_lv_first_chem AS chem
    ON b.stay_id = chem.stay_id
LEFT JOIN tmp_lv_first_coag AS coag
    ON b.stay_id = coag.stay_id
LEFT JOIN tmp_lv_first_bg AS bg
    ON b.stay_id = bg.stay_id
LEFT JOIN mimiciv_derived.first_day_rrt AS rrt
    ON b.stay_id = rrt.stay_id
LEFT JOIN tmp_lv_interventions AS iv
    ON b.stay_id = iv.stay_id
LEFT JOIN tmp_lv_sedatives AS sed
    ON b.stay_id = sed.stay_id;

CREATE UNIQUE INDEX official_analysis_population_stay_idx
    ON lowvoltage.official_analysis_population (stay_id);
CREATE INDEX official_analysis_population_subject_time_idx
    ON lowvoltage.official_analysis_population (subject_id, intime);
CREATE INDEX official_analysis_population_hadm_idx
    ON lowvoltage.official_analysis_population (hadm_id);
ANALYZE lowvoltage.official_analysis_population;

SELECT
    COUNT(*) AS n_icu_stays,
    COUNT(DISTINCT subject_id) AS n_subjects,
    COUNT(*) FILTER (WHERE crtr_sepsis3 = 1) AS n_sepsis3_stays,
    COUNT(*) FILTER (
        WHERE crtr_sepsis3 = 1
          AND day28_los IS NOT NULL
          AND day28_los > 0
          AND day28_outcome IS NOT NULL
          AND age >= 18
          AND COALESCE(icd_pregnancy, 0) = 0
          AND COALESCE(icd_malignancy, 0) = 0
          AND icu_hadm_order = 1
    ) AS n_non_ecg_eligible_stays
FROM lowvoltage.official_analysis_population;
