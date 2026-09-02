param(
    [string]$Database = "mimic4_v22",
    [string]$HostName = "localhost",
    [int]$Port = 5432,
    [string]$UserName = "postgres",
    [string]$Psql = "psql",
    [string]$Rscript = "Rscript",
    [string]$Python = "python",
    [switch]$SkipSql,
    [switch]$ForceImputation
)

$ErrorActionPreference = "Stop"
$scriptDirectory = Split-Path -Parent $MyInvocation.MyCommand.Path
$projectRoot = (Resolve-Path (Join-Path $scriptDirectory "..")).Path
$inputDirectory = Join-Path $projectRoot "data\input"
New-Item -ItemType Directory -Force -Path $inputDirectory | Out-Null

function Invoke-PsqlFile([string]$FileName) {
    & $Psql -X -v ON_ERROR_STOP=1 -h $HostName -p $Port -U $UserName -d $Database -f $FileName
    if ($LASTEXITCODE -ne 0) { throw "PostgreSQL script failed: $FileName" }
}

function Export-Query([string]$Query, [string]$OutputName) {
    $outputPath = (Join-Path $inputDirectory $OutputName).Replace("\", "/")
    $copyCommand = "\copy ($Query) TO '$outputPath' WITH (FORMAT CSV, HEADER TRUE)"
    & $Psql -X -v ON_ERROR_STOP=1 -h $HostName -p $Port -U $UserName -d $Database -c $copyCommand
    if ($LASTEXITCODE -ne 0) { throw "CSV export failed: $OutputName" }
}

if (-not $SkipSql) {
    $sqlFiles = @(
        "00_check_database_schema.sql",
        "01_build_official_analysis_population.sql",
        "02_extract_nononcologic_broad_cohort.sql",
        "03_extract_nononcologic_post_icu_24h_cohort.sql",
        "04_extract_malignancy_inclusive_broad_cohort.sql",
        "05_extract_malignancy_inclusive_post_icu_24h_cohort.sql",
        "06_build_nononcologic_landmark_extension.sql",
        "07_build_malignancy_inclusive_landmark_extension.sql",
        "08_create_waveform_validation_frame.sql",
        "09_final_database_qc.sql"
    )
    foreach ($sqlFile in $sqlFiles) {
        Invoke-PsqlFile (Join-Path $projectRoot "sql\$sqlFile")
    }

    Export-Query "SELECT * FROM lowvoltage.sepsis_ecg_lqrs_post_icu_24h_official ORDER BY subject_id, stay_id" "nononc_post_icu_24h.csv"
    Export-Query "SELECT * FROM lowvoltage.landmark_extension_20260720 ORDER BY subject_id, stay_id" "landmark_extension.csv"
    Export-Query "SELECT * FROM lowvoltage.sepsis_ecg_lqrs_post_icu_24h_malignancy_inclusive ORDER BY subject_id, stay_id" "malignancy_inclusive_post_icu_24h.csv"
    Export-Query "SELECT * FROM lowvoltage.sepsis_ecg_lqrs_broad_malignancy_inclusive ORDER BY subject_id, stay_id" "malignancy_inclusive_broad.csv"
    Export-Query "SELECT * FROM lowvoltage.landmark_extension_inclusive_20260804 ORDER BY subject_id, stay_id" "landmark_extension_inclusive.csv"
}

$env:LOWVOLTAGE_PROJECT_ROOT = $projectRoot
$env:LOWVOLTAGE_INPUT_DIR = $inputDirectory
$env:LOWVOLTAGE_FORCE_IMPUTE = if ($ForceImputation) { "1" } else { "0" }

& $Rscript (Join-Path $projectRoot "R\run_all.R")
if ($LASTEXITCODE -ne 0) { throw "The R analysis pipeline failed." }

& $Python (Join-Path $projectRoot "figure1\generate_scientific_cohort_flowchart.py") `
    --config (Join-Path $projectRoot "figure1\Figure_1_cohort_selection_config.json") `
    --output-dir (Join-Path $projectRoot "outputs\figures\figure1_bundle")
if ($LASTEXITCODE -ne 0) { throw "Figure 1 generation or QC failed." }

& $Rscript (Join-Path $projectRoot "R\07_postrun_qc.R")
if ($LASTEXITCODE -ne 0) { throw "Post-run QC failed." }

Write-Host "Pipeline completed successfully: $projectRoot"
