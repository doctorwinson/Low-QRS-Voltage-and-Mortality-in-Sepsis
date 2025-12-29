WITH ecg_data AS (
    SELECT 
        subject_id,
        ecg_time,
        CASE WHEN ARRAY[report_1, report_2, report_3, report_4, report_5, report_6, report_7] 
                  ~* ANY(ARRAY['Low Voltage', 'Low QRS Voltage'])
             THEN 1 
             ELSE 0 
        END AS low_qrs,
        CASE WHEN report_1 IS NOT NULL 
             THEN 1 
             ELSE 0 
        END AS has_valid_ecg
    FROM mimiciv_ecg.machine_measurements
    WHERE report_1 IS NOT NULL 
),
low_voltage_ecg AS (
    SELECT 
        subject_id,
        ecg_time AS lowqrstime,
        1 AS low_qrs
    FROM ecg_data
    WHERE low_qrs = 1
),
valid_ecg AS (
    SELECT 
        subject_id,
        ecg_time
    FROM ecg_data
    WHERE has_valid_ecg = 1
),
aki_patients AS (
    SELECT DISTINCT 
        stay_id,
        1 AS aki
    FROM mimiciv_derived.kdigo_stages
    WHERE aki_stage_smoothed > 0
),
population AS (
    SELECT 
        pop.*,
        lv.lowqrstime,
        lv.low_qrs,
        ve.ecg_time,
        aki.aki
    FROM bdmcc.bdmcc_population pop
    LEFT JOIN low_voltage_ecg lv 
        USING (subject_id)
    LEFT JOIN valid_ecg ve 
        USING (subject_id)
    LEFT JOIN aki_patients aki 
        USING (stay_id)
)
SELECT * 
FROM population;

