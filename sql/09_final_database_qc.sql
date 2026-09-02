DO $$
DECLARE
    source_n integer;
    sepsis_n integer;
    eligible_n integer;
    landmark_n integer;
    exposed_n integer;
    deaths_n integer;
BEGIN
    SELECT COUNT(*) INTO source_n FROM mimiciv_icu.icustays;
    SELECT COUNT(*) INTO sepsis_n
    FROM lowvoltage.official_analysis_population
    WHERE crtr_sepsis3 = 1;
    SELECT COUNT(*) INTO eligible_n
    FROM lowvoltage.landmark_extension_inclusive_20260804;
    SELECT COUNT(*),
           COUNT(*) FILTER (WHERE any_low_qrs_post_24h = 1),
           SUM(day28_outcome)
    INTO landmark_n, exposed_n, deaths_n
    FROM lowvoltage.landmark_extension_inclusive_20260804
    WHERE landmark_selection_group = 'included_landmark';

    IF source_n <> 73181 OR sepsis_n <> 32970 OR eligible_n <> 30007
       OR landmark_n <> 14129 OR exposed_n <> 3213 OR deaths_n <> 2346 THEN
        RAISE EXCEPTION
            'QC failure: source %, sepsis %, eligible %, landmark %, exposed %, deaths %',
            source_n, sepsis_n, eligible_n, landmark_n, exposed_n, deaths_n;
    END IF;
END $$;

SELECT
    landmark_selection_group,
    COUNT(*) AS stays,
    COUNT(DISTINCT subject_id) AS patients,
    COUNT(*) FILTER (WHERE any_low_qrs_post_24h = 1) AS exposed_stays,
    SUM(day28_outcome) AS deaths
FROM lowvoltage.landmark_extension_inclusive_20260804
GROUP BY landmark_selection_group
ORDER BY landmark_selection_group;
