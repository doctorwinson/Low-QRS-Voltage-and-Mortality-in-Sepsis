/*
Purpose:
  Build the complete non-oncologic landmark-selection population and derive a
  deterministic first post-ICU-admission ECG exposure for the 24-hour window.

Source versions:
  MIMIC-IV v2.2 and MIMIC-IV-ECG v1.0.

Key rules:
  - Search report_0 through report_17 for non-negated low-voltage phrases.
  - Match ECGs from ICU admission through 24 hours, no later than hospital
    discharge.
  - Assign each ECG to at most one eligible ICU stay.
  - Break first-ECG timestamp ties by study_id.
  - Preserve all 25,196 otherwise eligible stays so selection into the
    landmark cohort can be described directly.
*/

BEGIN;

DROP TABLE IF EXISTS lowvoltage.landmark_extension_20260720;
DROP TABLE IF EXISTS tmp_landmark_ecg_lines;
DROP TABLE IF EXISTS tmp_landmark_ecg_reports;
DROP TABLE IF EXISTS tmp_landmark_base;
DROP TABLE IF EXISTS tmp_landmark_ecg_candidates;
DROP TABLE IF EXISTS tmp_landmark_ecg_summary;

CREATE TEMP TABLE tmp_landmark_ecg_lines AS
SELECT
    mm.subject_id,
    mm.study_id,
    mm.ecg_time,
    r.report_seq,
    r.report_field,
    r.report_text,
    CASE
        WHEN (
            r.report_text ~* 'low[[:space:]-]+qrs[[:space:]-]+voltages?'
            OR r.report_text ~* 'low[[:space:]-]+voltages?'
        )
        AND r.report_text !~* '(no|not|without|absence[[:space:]-]+of|absent)[[:space:]-]+(low[[:space:]-]+)?(qrs[[:space:]-]+)?voltages?'
        THEN 1 ELSE 0
    END AS is_low_qrs_line
FROM mimiciv_ecg.machine_measurements AS mm
CROSS JOIN LATERAL (
    SELECT
        v.report_seq,
        v.report_field,
        btrim(v.report_text_raw) AS report_text
    FROM (
        VALUES
            (0,  'report_0',  mm.report_0),
            (1,  'report_1',  mm.report_1),
            (2,  'report_2',  mm.report_2),
            (3,  'report_3',  mm.report_3),
            (4,  'report_4',  mm.report_4),
            (5,  'report_5',  mm.report_5),
            (6,  'report_6',  mm.report_6),
            (7,  'report_7',  mm.report_7),
            (8,  'report_8',  mm.report_8),
            (9,  'report_9',  mm.report_9),
            (10, 'report_10', mm.report_10),
            (11, 'report_11', mm.report_11),
            (12, 'report_12', mm.report_12),
            (13, 'report_13', mm.report_13),
            (14, 'report_14', mm.report_14),
            (15, 'report_15', mm.report_15),
            (16, 'report_16', mm.report_16),
            (17, 'report_17', mm.report_17)
    ) AS v(report_seq, report_field, report_text_raw)
    WHERE NULLIF(btrim(v.report_text_raw), '') IS NOT NULL
) AS r;

CREATE INDEX tmp_landmark_ecg_lines_study_idx
    ON tmp_landmark_ecg_lines (study_id);
ANALYZE tmp_landmark_ecg_lines;

CREATE TEMP TABLE tmp_landmark_ecg_reports AS
SELECT
    subject_id,
    study_id,
    ecg_time,
    MAX(is_low_qrs_line)::integer AS low_qrs,
    COUNT(*) AS n_report_lines,
    string_agg(report_field, ', ' ORDER BY report_seq)
        FILTER (WHERE is_low_qrs_line = 1) AS low_qrs_report_fields,
    string_agg(report_text, ' | ' ORDER BY report_seq)
        FILTER (WHERE is_low_qrs_line = 1) AS low_qrs_report_text
FROM tmp_landmark_ecg_lines
GROUP BY subject_id, study_id, ecg_time;

CREATE INDEX tmp_landmark_ecg_reports_subject_time_idx
    ON tmp_landmark_ecg_reports (subject_id, ecg_time);
ANALYZE tmp_landmark_ecg_reports;

CREATE TEMP TABLE tmp_landmark_base AS
SELECT *
FROM lowvoltage.official_analysis_population
WHERE crtr_sepsis3 = 1
  AND day28_los IS NOT NULL
  AND day28_los > 0
  AND day28_outcome IS NOT NULL
  AND age >= 18
  AND COALESCE(icd_pregnancy, 0) = 0
  AND COALESCE(icd_malignancy, 0) = 0
  AND icu_hadm_order = 1;

CREATE UNIQUE INDEX tmp_landmark_base_stay_idx
    ON tmp_landmark_base (stay_id);
CREATE INDEX tmp_landmark_base_subject_time_idx
    ON tmp_landmark_base (subject_id, intime);
ANALYZE tmp_landmark_base;

CREATE TEMP TABLE tmp_landmark_ecg_candidates AS
SELECT *
FROM (
    SELECT
        b.stay_id,
        b.subject_id,
        e.study_id,
        e.ecg_time,
        e.low_qrs,
        e.n_report_lines,
        e.low_qrs_report_fields,
        e.low_qrs_report_text,
        EXTRACT(EPOCH FROM (e.ecg_time - b.intime)) / 3600.0 AS ecg_hours_from_icu_intime,
        CASE WHEN e.ecg_time > b.outtime THEN 1 ELSE 0 END AS ecg_after_icu_outtime,
        ROW_NUMBER() OVER (
            PARTITION BY e.study_id
            ORDER BY
                ABS(EXTRACT(EPOCH FROM (e.ecg_time - b.intime))),
                b.intime,
                b.stay_id
        ) AS ecg_to_stay_rank
    FROM tmp_landmark_base AS b
    INNER JOIN tmp_landmark_ecg_reports AS e
        ON e.subject_id = b.subject_id
       AND e.ecg_time >= b.intime
       AND e.ecg_time <= b.intime + INTERVAL '24 hours'
       AND e.ecg_time <= b.dischtime
) AS matched
WHERE ecg_to_stay_rank = 1;

CREATE INDEX tmp_landmark_ecg_candidates_stay_time_idx
    ON tmp_landmark_ecg_candidates (stay_id, ecg_time, study_id);
ANALYZE tmp_landmark_ecg_candidates;

CREATE TEMP TABLE tmp_landmark_ecg_summary AS
SELECT
    stay_id,
    COUNT(*) AS n_ecg_post_24h,
    MAX(low_qrs)::integer AS any_low_qrs_post_24h,
    SUM(low_qrs)::integer AS n_low_qrs_ecg_post_24h,
    (array_agg(study_id ORDER BY ecg_time, study_id))[1] AS first_post_ecg_study_id,
    (array_agg(ecg_time ORDER BY ecg_time, study_id))[1] AS first_post_ecg_time,
    (array_agg(low_qrs ORDER BY ecg_time, study_id))[1]::integer AS first_post_ecg_low_qrs,
    (array_agg(ecg_hours_from_icu_intime ORDER BY ecg_time, study_id))[1] AS first_post_ecg_hours_from_icu,
    (array_agg(ecg_after_icu_outtime ORDER BY ecg_time, study_id))[1]::integer AS first_post_ecg_after_icu_outtime,
    COUNT(*) FILTER (WHERE ecg_after_icu_outtime = 1) AS n_ecg_after_icu_outtime
FROM tmp_landmark_ecg_candidates
GROUP BY stay_id;

CREATE UNIQUE INDEX tmp_landmark_ecg_summary_stay_idx
    ON tmp_landmark_ecg_summary (stay_id);
ANALYZE tmp_landmark_ecg_summary;

CREATE TABLE lowvoltage.landmark_extension_20260720 AS
SELECT
    b.*,
    s.n_ecg_post_24h,
    s.any_low_qrs_post_24h,
    s.n_low_qrs_ecg_post_24h,
    s.first_post_ecg_study_id,
    s.first_post_ecg_time,
    s.first_post_ecg_low_qrs,
    s.first_post_ecg_hours_from_icu,
    s.first_post_ecg_after_icu_outtime,
    s.n_ecg_after_icu_outtime,
    CASE
        WHEN b.day28_outcome = 1 AND b.day28_los <= 24
            THEN 'death_by_24h'
        WHEN b.day28_los > 24 AND s.stay_id IS NOT NULL
            THEN 'included_landmark'
        WHEN b.day28_los > 24 AND s.stay_id IS NULL
            THEN 'alive_24h_no_post_ecg'
        ELSE 'other_exclusion'
    END AS landmark_selection_group
FROM tmp_landmark_base AS b
LEFT JOIN tmp_landmark_ecg_summary AS s
    ON s.stay_id = b.stay_id;

CREATE UNIQUE INDEX landmark_extension_20260720_stay_idx
    ON lowvoltage.landmark_extension_20260720 (stay_id);
CREATE INDEX landmark_extension_20260720_group_idx
    ON lowvoltage.landmark_extension_20260720 (landmark_selection_group);
ANALYZE lowvoltage.landmark_extension_20260720;

COMMIT;

SELECT
    landmark_selection_group,
    COUNT(*) AS stays,
    COUNT(DISTINCT subject_id) AS patients,
    COUNT(*) FILTER (WHERE any_low_qrs_post_24h = 1) AS any_positive_stays,
    COUNT(*) FILTER (WHERE first_post_ecg_low_qrs = 1) AS first_ecg_positive_stays,
    SUM(day28_outcome) AS deaths
FROM lowvoltage.landmark_extension_20260720
GROUP BY landmark_selection_group
ORDER BY landmark_selection_group;

SELECT
    COUNT(*) FILTER (WHERE landmark_selection_group = 'included_landmark') AS landmark_stays,
    COUNT(DISTINCT subject_id) FILTER (WHERE landmark_selection_group = 'included_landmark') AS landmark_patients,
    COUNT(*) FILTER (
        WHERE landmark_selection_group = 'included_landmark'
          AND n_ecg_post_24h = 1
    ) AS one_ecg_stays,
    COUNT(*) FILTER (
        WHERE landmark_selection_group = 'included_landmark'
          AND first_post_ecg_after_icu_outtime = 1
    ) AS first_ecg_after_icu_outtime
FROM lowvoltage.landmark_extension_20260720;
