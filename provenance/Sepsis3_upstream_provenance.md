# Sepsis-3 upstream provenance

The analysis reads the `mimiciv_derived.sepsis3` table and related official derived concepts installed in the local MIMIC-IV v2.2 PostgreSQL database. The available database did not retain the Git commit hash or view-definition checksum of the MIMIC Code checkout that originally generated these derived tables. Consequently, the exact historical upstream commit cannot be recovered from the database alone.

The code package fully specifies the study-specific cohort construction, ECG linkage, text phenotype, landmark eligibility, and downstream statistical analysis. Reproducing the derived concepts from raw MIMIC-IV tables first requires a compatible MIMIC Code build for MIMIC-IV v2.2. The package must not be described as proving bit-for-bit reproduction of the historical derived-table build.

The local pre-analysis checks require the expected schemas and derived tables, and `sql/09_final_database_qc.sql` verifies the study-specific row counts after reconstruction.
