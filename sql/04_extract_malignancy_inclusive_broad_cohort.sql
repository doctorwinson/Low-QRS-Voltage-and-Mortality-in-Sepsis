/*
Purpose:
  Build a malignancy-inclusive broad-window sensitivity cohort from the
  official MIMIC-derived analysis population and MIMIC-IV-ECG reports.

Key implementation details:
  1. Search all machine_measurements report_0 through report_17, not report_1
     alone.
  2. Keep ECG records at study_id level before applying the ICU time window.
  3. Join ECG records to ICU stays by subject_id and +/-24 hours around intime.
  4. Select one ECG per stay after time-window matching.
  5. Output field names expected by downstream R: ecg_time and low_qrs.

Exposure definition:
  low_qrs = 1 if any non-negated ECG machine report line within the selected ECG
  contains low QRS voltage / low voltage text.

Selection rule for one ECG per stay:
  - Prefer a low-QRS ECG within the ICU +/-24h window.
  - If multiple low-QRS ECGs exist, choose the earliest ECG time.
  - If no low-QRS ECG exists, choose the earliest ECG time.

This makes the stay-level group equivalent to "any machine-reported low QRS
voltage within the ICU +/-24h window".
*/

DROP TABLE IF EXISTS tmp_lqrs_ecg_report_lines;
DROP TABLE IF EXISTS tmp_lqrs_ecg_reports;
DROP TABLE IF EXISTS tmp_lqrs_eligible_base;
DROP TABLE IF EXISTS tmp_lqrs_ecg_candidates;
DROP TABLE IF EXISTS tmp_lqrs_ecg_window_summary;
DROP TABLE IF EXISTS tmp_lqrs_selected_ecg;
DROP TABLE IF EXISTS lowvoltage.sepsis_ecg_lqrs_broad_malignancy_inclusive;

CREATE TEMP TABLE tmp_lqrs_ecg_report_lines AS
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

CREATE INDEX tmp_lqrs_ecg_report_lines_study_idx
    ON tmp_lqrs_ecg_report_lines (study_id);
CREATE INDEX tmp_lqrs_ecg_report_lines_subject_time_idx
    ON tmp_lqrs_ecg_report_lines (subject_id, ecg_time);

ANALYZE tmp_lqrs_ecg_report_lines;

CREATE TEMP TABLE tmp_lqrs_ecg_reports AS
    SELECT
        subject_id,
        study_id,
        ecg_time,
        COUNT(*) AS n_report_lines,
        MAX(is_low_qrs_line) AS low_qrs,
        string_agg(report_text, ' | ' ORDER BY report_seq) AS full_report_text,
        string_agg(report_field, ', ' ORDER BY report_seq)
            FILTER (WHERE is_low_qrs_line = 1) AS low_qrs_report_fields,
        string_agg(report_text, ' | ' ORDER BY report_seq)
            FILTER (WHERE is_low_qrs_line = 1) AS low_qrs_report_text
    FROM tmp_lqrs_ecg_report_lines
    GROUP BY subject_id, study_id, ecg_time;

CREATE INDEX tmp_lqrs_ecg_reports_subject_time_idx
    ON tmp_lqrs_ecg_reports (subject_id, ecg_time);
CREATE INDEX tmp_lqrs_ecg_reports_study_idx
    ON tmp_lqrs_ecg_reports (study_id);

ANALYZE tmp_lqrs_ecg_reports;

CREATE TEMP TABLE tmp_lqrs_eligible_base AS
    SELECT *
    FROM lowvoltage.official_analysis_population
    WHERE crtr_sepsis3 = 1
      AND day28_los IS NOT NULL
      AND day28_los > 0
      AND day28_outcome IS NOT NULL
      AND age >= 18
      AND COALESCE(icd_pregnancy, 0) = 0
      AND icu_hadm_order = 1;

CREATE INDEX tmp_lqrs_eligible_base_subject_time_idx
    ON tmp_lqrs_eligible_base (subject_id, intime);
CREATE INDEX tmp_lqrs_eligible_base_stay_idx
    ON tmp_lqrs_eligible_base (stay_id);

ANALYZE tmp_lqrs_eligible_base;

CREATE TEMP TABLE tmp_lqrs_ecg_candidates AS
SELECT *
FROM (
    SELECT
        b.stay_id,
        b.subject_id,
        e.study_id AS ecg_study_id,
        e.ecg_time,
        e.low_qrs,
        e.n_report_lines,
        e.full_report_text,
        e.low_qrs_report_fields,
        e.low_qrs_report_text,
        EXTRACT(EPOCH FROM (e.ecg_time - b.intime)) / 3600.0 AS ecg_hours_from_icu_intime,
        ROW_NUMBER() OVER (
            PARTITION BY e.study_id
            ORDER BY ABS(EXTRACT(EPOCH FROM (e.ecg_time - b.intime))), b.intime, b.stay_id
        ) AS ecg_to_stay_rank
    FROM tmp_lqrs_eligible_base AS b
    INNER JOIN tmp_lqrs_ecg_reports AS e
        ON b.subject_id = e.subject_id
       AND e.ecg_time >= b.intime - INTERVAL '24 hours'
       AND e.ecg_time <= b.intime + INTERVAL '24 hours'
       AND e.ecg_time <= b.dischtime
) AS matched
WHERE ecg_to_stay_rank = 1;

CREATE INDEX tmp_lqrs_ecg_candidates_stay_idx
    ON tmp_lqrs_ecg_candidates (stay_id);

ANALYZE tmp_lqrs_ecg_candidates;

CREATE TEMP TABLE tmp_lqrs_ecg_window_summary AS
    SELECT
        stay_id,
        COUNT(*) AS n_ecg_in_window,
        SUM(low_qrs)::integer AS n_low_qrs_ecg_in_window,
        MIN(ecg_time) AS first_ecg_time_in_window,
        MIN(ecg_time) FILTER (WHERE low_qrs = 1) AS first_low_qrs_time_in_window
    FROM tmp_lqrs_ecg_candidates
    GROUP BY stay_id;

CREATE INDEX tmp_lqrs_ecg_window_summary_stay_idx
    ON tmp_lqrs_ecg_window_summary (stay_id);

ANALYZE tmp_lqrs_ecg_window_summary;

CREATE TEMP TABLE tmp_lqrs_selected_ecg AS
SELECT *
FROM (
    SELECT
        c.*,
        s.n_ecg_in_window,
        s.n_low_qrs_ecg_in_window,
        s.first_ecg_time_in_window,
        s.first_low_qrs_time_in_window,
        ROW_NUMBER() OVER (
            PARTITION BY c.stay_id
            ORDER BY c.low_qrs DESC, c.ecg_time ASC, c.ecg_study_id ASC
        ) AS ecg_rank_within_stay
    FROM tmp_lqrs_ecg_candidates AS c
    INNER JOIN tmp_lqrs_ecg_window_summary AS s
        ON c.stay_id = s.stay_id
    ) AS ranked
WHERE ecg_rank_within_stay = 1;

CREATE INDEX tmp_lqrs_selected_ecg_stay_idx
    ON tmp_lqrs_selected_ecg (stay_id);

ANALYZE tmp_lqrs_selected_ecg;

CREATE TABLE lowvoltage.sepsis_ecg_lqrs_broad_malignancy_inclusive AS
SELECT
    b.*,
    se.ecg_study_id,
    se.ecg_time,
    se.low_qrs,
    se.ecg_hours_from_icu_intime,
    se.n_report_lines,
    se.n_ecg_in_window,
    se.n_low_qrs_ecg_in_window,
    se.first_ecg_time_in_window,
    se.first_low_qrs_time_in_window,
    se.low_qrs_report_fields,
    se.low_qrs_report_text,
    se.full_report_text,
    'prefer_low_qrs_then_earliest_ecg_in_window' AS ecg_selection_rule
FROM tmp_lqrs_eligible_base AS b
INNER JOIN tmp_lqrs_selected_ecg AS se
    ON b.stay_id = se.stay_id;

CREATE UNIQUE INDEX sepsis_ecg_lqrs_broad_malignancy_inclusive_stay_idx
    ON lowvoltage.sepsis_ecg_lqrs_broad_malignancy_inclusive (stay_id);
ANALYZE lowvoltage.sepsis_ecg_lqrs_broad_malignancy_inclusive;

SELECT
    COUNT(*) AS stays,
    COUNT(DISTINCT subject_id) AS patients,
    SUM(low_qrs) AS exposed_stays,
    SUM(day28_outcome) AS deaths
FROM lowvoltage.sepsis_ecg_lqrs_broad_malignancy_inclusive;
