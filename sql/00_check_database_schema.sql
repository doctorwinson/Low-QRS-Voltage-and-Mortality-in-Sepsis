/*
Purpose:
  Check the available report_* columns in mimiciv_ecg.machine_measurements.

Why this exists:
  MIMIC-IV-ECG machine interpretation text is stored across report_0,
  report_1, ..., report_17. The previous scripts only searched report_1,
  which can miss low-voltage diagnoses recorded in other report fields.

Run this first and confirm that report_0 through report_17 exist before
running the extraction scripts.
*/

SELECT
    table_schema,
    table_name,
    column_name,
    data_type,
    ordinal_position
FROM information_schema.columns
WHERE table_schema = 'mimiciv_ecg'
  AND table_name = 'machine_measurements'
  AND column_name ~ '^report_[0-9]+$'
ORDER BY
    substring(column_name from '[0-9]+')::integer;
