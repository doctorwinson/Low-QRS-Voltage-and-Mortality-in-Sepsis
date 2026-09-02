# Low QRS Voltage and Mortality in Sepsis

Reproducible analysis code for the cohort study **"Machine-reported low QRS voltage and 28-day mortality among 24-hour survivors with sepsis: a MIMIC-IV-ECG cohort study."**

## Scope

This package reproduces the study-specific cohort extraction, machine-report text phenotype, 24-hour landmark analysis, multiple imputation, statistical models, tables, figures, and quality-control checks for the manuscript on machine-reported low QRS voltage and 28-day mortality in sepsis.

Patient-level MIMIC data are not included. The `data/input` and `local_cache` directories are local working directories and must remain excluded from any public or submission archive.

## Data versions

- MIMIC-IV version 2.2
- MIMIC-IV-ECG version 1.0
- PostgreSQL schemas: `mimiciv_hosp`, `mimiciv_icu`, `mimiciv_derived`, and `mimiciv_ecg`

Required official derived tables include `age`, `bg`, `charlson`, `chemistry`, `coagulation`, `complete_blood_count`, `first_day_rrt`, `first_day_sofa`, `first_day_weight`, `icustay_detail`, `sapsii`, `sepsis3`, `vasoactive_agent`, `ventilation`, and `vitalsign`.

The exact historical MIMIC Code commit used to create the locally installed derived tables could not be recovered. See `provenance/Sepsis3_upstream_provenance.md`. The package therefore supports complete reproduction of the study-specific downstream workflow from compatible official derived tables, but it does not claim a bit-for-bit reconstruction of the historical upstream derived-table build.

## Software

The verified run used PostgreSQL 16, R 4.5.2, and Python 3. Package versions from the verified R run are recorded in `environment/sessionInfo.txt`. The main R dependencies are `broom`, `car`, `dplyr`, `ggplot2`, `lattice`, `MASS`, `mice`, `mitml`, `patchwork`, `readr`, `sandwich`, `splines`, `survival`, `survminer`, `tibble`, and `tidyr`. Figure 1 additionally requires Matplotlib; editable PPTX export requires `python-pptx`. The optional LZW TIFF conversion uses the R package `magick`.

## Reproduction from an empty working directory

1. Copy this code package to an empty directory.
2. Install MIMIC-IV v2.2 and MIMIC-IV-ECG v1.0 in PostgreSQL and build the required official derived concepts.
3. Set the PostgreSQL password through the standard `PGPASSWORD` environment variable or a local password file. Do not place credentials in this package.
4. Run the PowerShell entry point from the package root:

```powershell
.\run\run_pipeline.ps1 -Database mimic4_v22 -ForceImputation
```

If `psql`, `Rscript`, or `python` is not on `PATH`, pass the executable paths through `-Psql`, `-Rscript`, or `-Python`. The script resolves every project file relative to its own location; no local absolute path is required.

For an R-only rerun after the five local CSV inputs have already been exported, use:

```powershell
.\run\run_pipeline.ps1 -SkipSql -ForceImputation
```

## Execution order

The master script runs the SQL files in numeric order, exports the local analysis inputs, then runs `R/run_all.R`. The R workflow performs 40 imputations with 10 iterations and predictive mean matching with five donors, fits the primary and sensitivity models, calculates marginal standardized risks, tests the lactate functional form with an MI-pooled D2 comparison, creates diagnostics and tables, and generates Figures 2 and 3. The bundled fixed-coordinate flowchart generator creates Figure 1 and runs arithmetic and text-overflow QC. `R/07_postrun_qc.R` performs final numerical checks and creates the white-background RGB, LZW-compressed Figure 1 TIFF.

## Main inputs and outputs

Local inputs are written to `data/input`. They contain restricted patient-level data and must not be uploaded. Aggregate results are written to `outputs/results`, diagnostics to `outputs/diagnostics`, and figures to `outputs/figures`. The final post-run report is `outputs/postrun_QC_report.md`.

The `expected_outputs` directory contains the non-identifying aggregate CSV summaries and QC report from the verified clean-directory run. These files provide reference values for checking a new installation; they do not contain patient-level rows or identifier values.

The waveform validation workflow created a deterministic 200-positive/200-negative sampling frame. Two blinded readers reviewed all 400 ECGs, and a third physician adjudicated 23 discordant classifications. The completed protocol, data dictionary, aggregate confusion matrix, and aggregate performance estimates are provided in `validation`. The restricted linkage key, waveforms, neutral sample codes, reader forms, and record-level classifications are intentionally excluded.

## Analysis definitions

The primary estimand is the conditional odds ratio for death after the 24-hour landmark through day 28 among patients alive at 24 hours who had at least one ECG obtained in routine care during the first 24 hours after ICU admission. Exposure is any non-negated machine-report mention of low voltage across `report_0` through `report_17`; it is not a waveform-confirmed or cardiologist-confirmed diagnosis. The main model uses fixed natural-spline bases for age, weight, and SOFA score and subject-level cluster-robust standard errors. Lactate is linear in the main model; fixed-spline and log1p sensitivity analyses are included. Rubin pooling uses finite-sample degrees of freedom. Cox models are secondary time-to-event summaries.

## Submission safety

Before creating a submission ZIP, exclude `data/input`, `local_cache`, `local_only_validation`, all RDS files, and every file containing a MIMIC `subject_id`, `hadm_id`, `stay_id`, ECG `study_id`, or restricted waveform path. Only scripts, documentation, aggregate output summaries, and non-identifying figures may be shared.

## Data access

The source data are available from PhysioNet to credentialed users who complete the required training and data-use agreement. This repository does not redistribute MIMIC-IV, MIMIC-IV-ECG, waveforms, record-level validation classifications, linkage keys, or patient-level analysis extracts.

## Repository version

This public package corresponds to the analysis finalized on September 2, 2026. The `expected_outputs` directory contains aggregate reference results from the verified run so that users can compare a local reproduction against the manuscript analysis.
