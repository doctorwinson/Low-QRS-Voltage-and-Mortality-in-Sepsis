-- Create a reproducible, local-only waveform validation sample.
--
-- This table contains restricted MIMIC-IV-ECG record paths and study IDs.
-- It must not be included in a manuscript submission or a public repository.

DROP TABLE IF EXISTS lowvoltage.waveform_validation_frame_20260804;

CREATE TABLE lowvoltage.waveform_validation_frame_20260804 AS
WITH eligible_ecgs AS (
    SELECT DISTINCT ON (cohort.ecg_study_id)
        cohort.ecg_study_id,
        cohort.low_qrs AS machine_report_class,
        records.path AS waveform_path,
        cohort.ecg_time,
        cohort.low_qrs_report_fields,
        cohort.full_report_text
    FROM lowvoltage.sepsis_ecg_lqrs_post_icu_24h_malignancy_inclusive AS cohort
    INNER JOIN lowvoltage.landmark_extension_inclusive_20260804 AS extension
        ON extension.stay_id = cohort.stay_id
       AND extension.landmark_selection_group = 'included_landmark'
    INNER JOIN mimiciv_ecg.record_list AS records
        ON records.study_id = cohort.ecg_study_id
    WHERE cohort.ecg_study_id IS NOT NULL
      AND cohort.low_qrs IN (0, 1)
    ORDER BY cohort.ecg_study_id, cohort.stay_id
),
ranked AS (
    SELECT
        eligible_ecgs.*,
        ROW_NUMBER() OVER (
            PARTITION BY machine_report_class
            ORDER BY md5(ecg_study_id::text || '|lowvoltage-validation-20260804')
        ) AS stratum_rank
    FROM eligible_ecgs
),
sampled AS (
    SELECT *
    FROM ranked
    WHERE stratum_rank <= 200
),
blinded AS (
    SELECT
        sampled.*,
        ROW_NUMBER() OVER (
            ORDER BY md5(ecg_study_id::text || '|lowvoltage-blinding-20260804')
        ) AS blind_rank
    FROM sampled
)
SELECT
    'V' || LPAD(blind_rank::text, 3, '0') AS sample_code,
    ecg_study_id,
    waveform_path,
    ecg_time,
    machine_report_class,
    low_qrs_report_fields,
    full_report_text
FROM blinded
ORDER BY blind_rank;

ALTER TABLE lowvoltage.waveform_validation_frame_20260804
    ADD PRIMARY KEY (sample_code);

DO $$
DECLARE
    positive_n integer;
    negative_n integer;
BEGIN
    SELECT COUNT(*) FILTER (WHERE machine_report_class = 1),
           COUNT(*) FILTER (WHERE machine_report_class = 0)
    INTO positive_n, negative_n
    FROM lowvoltage.waveform_validation_frame_20260804;

    IF positive_n <> 200 OR negative_n <> 200 THEN
        RAISE EXCEPTION
            'Validation sample incomplete: positive %, negative %',
            positive_n, negative_n;
    END IF;
END $$;

COMMENT ON TABLE lowvoltage.waveform_validation_frame_20260804 IS
    'Restricted local-only 1:1 stratified validation sample; never upload publicly.';
