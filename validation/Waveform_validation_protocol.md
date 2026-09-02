# Waveform validation protocol

## Purpose

This substudy estimated the classification performance of the machine-report text phenotype against an adjudicated waveform-based low-QRS-voltage reference standard.

## Sampling

The study SQL selected 200 machine-report-positive and 200 machine-report-negative classifying ECGs from the primary 24-hour landmark cohort. Selection was deterministic within each stratum using an MD5 ordering seed. A second deterministic shuffle assigned neutral `V001`-`V400` codes that did not reveal machine-report status. The local linkage key contains MIMIC-IV-ECG study IDs and restricted waveform paths and is not included in this package.

## Reference standard

Readers assessed 12-lead waveforms rendered from the raw WFDB records at standard calibration. Low QRS voltage was present when peak-to-peak QRS amplitude was less than 0.5 mV in every limb lead or less than 1.0 mV in every precordial lead. Limb-lead and precordial-lead criteria were recorded separately. ECGs with inadequate signal quality, missing leads, nonstandard calibration, or artifacts that prevented reliable amplitude measurement could be classified as uninterpretable rather than forced into a positive or negative category.

## Blinding and adjudication

Two readers independently assessed all 400 ECGs while blinded to machine-report classification, automated screening, clinical covariates, and outcome. They recorded waveform-defined classification, lead territory, interpretability, calibration, rhythm or pacing, and a brief reason for uncertainty. A third physician independently reviewed the 23 discordant ECGs while blinded to both initial classifications and the same index-test and clinical information; the adjudicator's classification became the final reference for those recordings.

Every reader who accessed a raw waveform or rendered ECG did so under the PhysioNet data-use agreement. The coordinator provided files labelled only by the neutral sample code. The report-class linkage key, clinical outcomes, and full machine-report text remained unavailable until the reference classifications were locked.

## Statistical analysis

The machine-report phenotype was the index test and the adjudicated waveform classification was the reference standard. Positive predictive value and negative predictive value were estimated directly within the sampled machine-report-positive and machine-report-negative strata. Because the validation sample contained equal numbers of index-test-positive and index-test-negative ECGs rather than the source-cohort ratio, sensitivity, specificity, overall agreement, and index-test-versus-reference Cohen's kappa were standardized to the 3,213 machine-report-positive and 10,916 machine-report-negative classifications in the primary cohort.

A stratified nonparametric bootstrap with 10,000 repetitions, resampling separately within the two machine-report strata, provided 95% confidence intervals for the standardized metrics. Exact binomial intervals were used for the predictive values. Pre-adjudication inter-reader agreement and three-category Cohen's kappa were reported separately from agreement between the machine report and the final waveform reference. No ECG was classified as uninterpretable, so the prespecified uninterpretable-record bounds were identical to the complete-reference analysis.

## Aggregate results

The two readers agreed on 377 of 400 ECGs (94.3%) before adjudication; three-category Cohen's kappa was 0.878 (95% bootstrap confidence interval, 0.829-0.925). After adjudication, 143 ECGs were waveform positive and 257 were waveform negative. The sample confusion matrix contained 140 true positives, 60 false positives, 3 false negatives, and 197 true negatives. Source-standardized sensitivity was 93.2%, specificity 91.8%, overall agreement 92.0%, and machine-report-versus-waveform kappa 0.751. Positive and negative predictive values were 70.0% and 98.5%, respectively. Aggregate numerical results are supplied in the two CSV files in this directory. Record-level classifications, neutral sample codes, waveforms, identifiers, and the linkage key are intentionally excluded.
