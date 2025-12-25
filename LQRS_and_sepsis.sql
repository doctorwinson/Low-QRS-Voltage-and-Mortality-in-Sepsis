/*
================================================================================
SCRIPT NAME: LQRS_and_sepsis.sql
DESCRIPTION:
    SQL script for extracting the analytical cohort from the MIMIC-IV database.
    The script identifies patients with Low QRS Voltage, retrieves ECG records.

STRUCTURE:
    1. low_qrs_ecg:       Identify patients with 'Low Voltage' or 'Low QRS Voltage'.
    2. all_ecg_records:   Identify presence of any valid ECG measurement.
    3. population_cohort: Merge clinical features with the base population.

NOTES:
    - Database: MIMIC-IV
    - ECG Definition: Text search in machine reports.
    
    
ACKNOWLEDGEMENT & DATA PROVENANCE NOTE
    We express our gratitude to the Ascetic Practitioners in Critical Care (APCC)
    team and the easy Data Science for Medicine (easyDSM) team for their generous
    sharing of expertise and code related to critical care big data, as well as
    the cross-platform Big Data Master of Critical Care (BDMCC) software
    (available at https://github.com/ningyile/BDMCC_APP).

    IMPORTANT:
    The table `bdmcc.bdmcc_population` referenced in this SQL script is a
    BDMCC-derived base cohort table, generated/provided by the BDMCC software as
    part of its standardized cohort construction workflow on MIMIC-IV.

    We also especially appreciate the MIMIC official team’s efforts to
    open-source the database and supporting code.
================================================================================
*/

-- 1. Identify patients with 'Low Voltage' or 'Low QRS Voltage' ECG reports
WITH low_qrs_ecg AS (
    SELECT
        subject_id,
        ecg_time AS low_qrs_time, -- Timestamp of the specific low voltage event
        1 AS has_low_qrs          -- Binary flag: 1 = Presence of Low QRS Voltage
    FROM
        mimiciv_ecg.machine_measurements
    WHERE
        -- Case-insensitive regex match for clinical keywords
        report_1 ~* 'Low Voltage' OR report_1 ~* 'Low QRS Voltage'
),

-- 2. Identify all patients who have any ECG measurement on file
-- Serves as a control/baseline check for ECG availability
all_ecg_records AS (
    SELECT
        subject_id,
        ecg_time
    FROM
        mimiciv_ecg.machine_measurements
    WHERE
        report_1 IS NOT NULL
),

-- 3. Create the final population by joining all datasets
population_cohort AS (
    SELECT
        base_population.*,
        
        -- ECG specific columns
        low_qrs_ecg.low_qrs_time,
        low_qrs_ecg.has_low_qrs,
        
        -- General ECG columns
        all_ecg_records.ecg_time AS any_ecg_time
        
    FROM
        (
            -- Load the base sepsis population
            SELECT *
            FROM bdmcc.bdmcc_population
        ) AS base_population
        
    -- Link to Low QRS Voltage events
    -- Note: LEFT JOIN retains all patients from the base cohort
    LEFT JOIN low_qrs_ecg
        USING (subject_id)
        
    -- Link to All ECG records
    -- Note: Joining on subject_id may result in multiple rows per patient
    -- corresponding to each ECG record found.
    LEFT JOIN all_ecg_records
        USING (subject_id)
)

-- Final Output: Select all columns for analysis
SELECT *
FROM population_cohort;
