# Simultaneous Jackknife Empirical Likelihood

R code accompanying the manuscript **“Simultaneous Jackknife Empirical Likelihood Inference for Comparing Multiple Classifiers in Ordered Three-Class Problems.”**

## Objective

When several fitted classifiers are evaluated on the same held-out subjects, their performance estimates are dependent. Ranking the observed values does not establish whether classifier differences exceed test-sample uncertainty, and repeated unadjusted pairwise tests can inflate the familywise error rate.

This repository implements simultaneous inference for classifier-specific volumes under three-class ROC surfaces (VUSs). It treats the vector of VUS estimators as a vector-valued three-sample U-statistic and provides:

- An omnibus jackknife empirical likelihood (JEL) test of equal VUSs.
- Pairwise JEL tests with Holm and Bonferroni adjustment.
- Dependence-aware max-JEL inference using Gaussian multiplier calibration.
- Simultaneous profile-JEL confidence intervals for pairwise VUS differences.

Inference is conditional on the fitted scoring rules. Variability from drawing new training samples, retraining, or repeating model selection is outside this inferential target.

## Repository structure

| File | Purpose |
|---|---|
| `R/00_self_check.R` | Script parsing, core arithmetic, biomarker-data integrity, and revision-helper checks |
| `R/01_jel_core.R` | VUS estimation, pooled pseudo-values, JEL tests, class-centered Wald inference, multiplier calibration, and certified profile interval inversion |
| `R/02_models.R` | Base-R classifiers, preprocessing, and training-only hyperparameter selection |
| `R/03_real_data_application.R` | Cohort-separated urinary-biomarker analysis |
| `R/04_simulation_study.R` | Ten complete-null and seven alternative scenarios |
| `R/05_interval_simulation.R` | Five scenarios evaluating simultaneous interval coverage and length |
| `R/06_make_figures.R` | Main figures and wide-format plotting tables derived from the aggregate results |
| `R/07_revision_simulations.R` | Twelve additional imbalance, dependence, and classifier-count scenarios, paired power comparisons, and controlled timing |
| `R/08_revision_figures.R` | Additional figures in PDF and PNG, and a TeX results table |
| `R/09_validate_revision.R` | Independent numerical checks against enumeration, exact likelihood examples, and alternative reference calculations |
| `R/99_run_all.R` | Terminal and RStudio workflow driver |

## Software requirements

- R version 3.6.0 or later. The completed manuscript run used R 4.4.2.
- Standard packages supplied with R, including `parallel`. No contributed R packages are required.

Run the commands below from the repository root. On macOS and Linux, `MC_CORES` controls parallel simulation workers. Windows uses serial execution. Classifier fitting, validation, figure generation, and controlled timing are serial.

## Data

The application uses urinary-biomarker data associated with:

> Debernardi S, O'Brien H, Algahmdi AS, et al. (2020). A combination of urinary biomarker panel and PancRISK score for earlier detection of pancreatic cancer: A case-control study. *PLOS Medicine*, 17(12), e1003489. https://doi.org/10.1371/journal.pmed.1003489

The analysis CSV is included at:

```text
data/urinary_biomarkers_pancreatic_cancer.csv
```

Source data are available through the [official S1 Table](https://doi.org/10.1371/journal.pmed.1003489.s009) and the [Kaggle dataset](https://www.kaggle.com/datasets/johnjdavisiv/urinary-biomarkers-for-pancreatic-cancer). Any replacement CSV must preserve the expected columns and coding:

```text
sample_id, patient_cohort, age, diagnosis,
creatinine, LYVE1, REG1B, TFF1
```

`patient_cohort` is coded 1 or 2. `diagnosis` is coded 1 for healthy controls, 2 for benign hepatobiliary disease, and 3 for pancreatic ductal adenocarcinoma (PDAC). The five predictors must be complete and finite; the four biomarkers must be strictly positive for log transformation.

| Cohort | Controls | Benign disease | PDAC | Total | Use |
|---|---:|---:|---:|---:|---|
| 1 | 81 | 89 | 162 | 332 | Preprocessing, five-fold tuning, and model fitting |
| 2 | 102 | 119 | 37 | 258 | Common held-out evaluation sample |
| Total | 183 | 208 | 199 | 590 | |

Preprocessing is estimated within each training fold during tuning, then refitted on all of cohort 1. Cohort-2 outcomes are not used for fitting or tuning. Unique sample identifiers are checked, but do not by themselves establish person-level independence across records.

Every driver mode runs the data integrity check, so the CSV is required even for `check` and `revision`.

## Running the analysis

Start with the integrity and independent numerical checks:

```bash
Rscript R/99_run_all.R check
```

Run a small end-to-end pipeline check:

```bash
Rscript R/99_run_all.R quick
```

Quick mode uses 10 replications per main scenario, 5 per interval scenario, 5 per additional scenario, 2 repetitions per timing configuration, and 199 multiplier draws for each analysis. Training-only tuning remains enabled by default. Quick results are pipeline checks and must not be used as manuscript estimates.

Run the application or selected experiment groups:

```bash
Rscript R/99_run_all.R real
Rscript R/99_run_all.R original
Rscript R/99_run_all.R revision
```

| Mode | Work performed after integrity and numerical validation |
|---|---|
| `check` | No application or full simulation study |
| `quick` | Application and all experiment groups with reduced budgets |
| `real` | Application only |
| `original` | Application, 17 main scenarios, interval study, and main figures |
| `revision` | Twelve additional scenarios, timing study, and additional figures |
| `full` | Application and all experiment groups and figures |

Run the complete workflow with the manuscript output suffix on macOS or Linux:

```bash
OUTPUT_SUFFIX=_full_R MC_CORES=1 Rscript R/99_run_all.R full
```

To parallelize the simulation replications, increase the worker count, for example:

```bash
OUTPUT_SUFFIX=_full_R MC_CORES=8 Rscript R/99_run_all.R full
```

Controlled timing remains serial. Timing results depend on the hardware and software environment.

From R or RStudio on any platform:

```r
# Set the working directory to the repository root first.
Sys.setenv(JEL_PROJECT_ROOT = normalizePath("."))
source("R/99_run_all.R")

run_jel_workflow("check")
run_jel_workflow("quick")

Sys.setenv(OUTPUT_SUFFIX = "_full_R", MC_CORES = "1")
run_jel_workflow("full")
```

Sourcing the driver defines `run_jel_workflow()` without starting a long run. Existing environment variables override mode defaults; check that overrides from earlier experiments are appropriate before launching a full run. Leave `INT_REPS` unset to preserve the scenario-specific interval budgets.

## Full-run budgets

With the default replication and multiplier settings:

| Component | Settings | Replications | Multiplier draws per replication |
|---|---:|---:|---:|
| Main simulations | 17 | 5,000 per setting; 85,000 total | 2,000 |
| Endpoint-inversion study | 5 | N3: 1,000; N4, N6, N9, A3: 500 each; 3,000 total | 2,000 |
| Additional simulations | 12 | 5,000 per setting; 60,000 total | 1,999 |
| Controlled timing | 6 | 500 per configuration; 3,000 total | 1,999 |

The application uses 20,000 multiplier draws. The inference studies contain 148,000 replications in total, excluding the separate timing repetitions and the validator's small deterministic checks. The full workflow is computationally intensive.

## Configuration

These are driver defaults outside quick mode, except where noted. Existing environment variables take precedence.

| Variable | Default | Description |
|---|---|---|
| `RETUNE` | `true` | Repeat five-fold training-cohort hyperparameter selection; `false` uses fixed hyperparameters and reports CV losses as `NA` |
| `REAL_SEED` | `2026` | Application seed |
| `REAL_B_MULT` | `20000` | Application multiplier draws |
| `R_REPS` | `5000` | Replications per main scenario |
| `INT_REPS` | Unset | Optional override for every interval scenario; unset preserves 1,000/500 budgets |
| `B_MULT` | `2000` | Main and interval-study multiplier draws |
| `SIM_SEED` | `2026` | Main simulation seed |
| `INT_SEED` | `2026` | Interval-study seed |
| `REV_REPS` | `5000` | Replications per additional scenario |
| `REV_B_MULT` | `1999` | Additional-simulation and controlled-timing multiplier draws |
| `REV_SEED` | `20261016` | Additional-simulation seed |
| `REV_RUNTIME_REPS` | `500` | Repetitions per controlled timing configuration |
| `REV_RESUME` | `true` | Resume compatible additional-simulation checkpoints |
| `REV_CHECKPOINT_EVERY` | `100` | Replications per additional-simulation checkpoint batch |
| `MC_CORES` | `1` | Simulation workers on non-Windows systems |
| `ALPHA` | `0.05` | Test significance level and nominal interval error probability |
| `OUTPUT_SUFFIX` | `_R` | Result and figure suffix; quick mode defaults to `_quick`, and check mode to `_check` |
| `JEL_PROJECT_ROOT` | Auto-detected | Optional absolute repository path |

The completed manuscript outputs use `_full_R`, which is an explicit override rather than the driver's default suffix. Use a separate suffix for runs with different settings.

Only the additional-simulation scenario replications support checkpoint/resume. Their saved configuration must match the requested settings, R version, and relevant source checksums. If it does not, choose a new suffix or deliberately restart with `REV_RESUME=false`. Earlier workflow stages and controlled timing are rerun when their mode is invoked.

For a complete `revision` or `full` workflow, retain all 12 additional scenarios: the additional figure script requires the complete scenario set.

## Outputs and reproducibility records

Scripts create `results/` and `figures/` as needed. Outputs include:

- Scenario definitions, seeds, replicate records, valid-replication counts, and numerical-failure indicators.
- Aggregate rejection probabilities, Monte Carlo standard errors, and failure bounds where implemented.
- Simultaneous interval coverage and mean length.
- Cross-validation folds and losses, fitted models, held-out probabilities and scores, and classifier performance summaries.
- Omnibus and pairwise JEL results, adjusted p-values, and simultaneous profile intervals.
- Additional-simulation summaries, paired power gains, and timing measurements.
- Four main PNG figures and three additional figures in both PDF and PNG, plus `revision_tables{suffix}.tex`.

`R/06_make_figures.R` constructs the wide null and alternative plotting tables directly from `simulation_aggregate_results{suffix}.csv`. No manual conversion is required. Plotting inputs must use the selected output suffix.

Each driver run writes configuration, session information, source MD5 checksums, stage logs, and exit statuses under a unique directory in `results/logs/`. A `COMPLETED.txt` marker is written only after all selected stages finish successfully. The numerical validator writes `results/validation_report_R.txt`; this filename is fixed and is replaced when validation is rerun.

## Monte Carlo uncertainty

The main null-rejection and power figures generated by `R/06_make_figures.R` use **pointwise 95% normal-approximation Monte Carlo intervals**: the estimated rate plus or minus 1.96 times its Monte Carlo standard error, clipped to the probability range [0, 1].

The additional probability figures generated by `R/08_revision_figures.R` use **pointwise 95% Wilson Monte Carlo intervals**. Intervals for max-JEL minus Holm-JEL power gains use the standard error of the paired replicate-wise indicator differences.

These intervals describe Monte Carlo uncertainty in estimated simulation probabilities or their differences. They are separate from the JEL confidence intervals for VUS differences. Neither overlapping nor nonoverlapping marginal error bars replace a paired comparison of methods.

## Statistical interpretation

- The diagnostic classes must have a scientifically meaningful order. The application uses healthy control < benign disease < PDAC as a diagnostic ordering, not a biological progression claim.
- All classifiers must be evaluated on the same held-out subjects, with independence across subjects and the class-specific samples assumed by the theory. Cross-validation folds must not be treated as independent test samples.
- VUS measures ordered-class discrimination. It does not establish probability calibration, clinical utility, fairness, or transportability.
- With K classifiers, the complete pairwise contrast matrix has rank K - 1, while the covariance can have lower rank if additional degeneracy is present. Gaussian multiplier calibration permits structural singularity of the full pairwise covariance and requires positive retained marginal variances.
- The theoretical guarantees are pointwise asymptotic results for fixed K, positive limiting class proportions, and the stated nondegeneracy conditions. They are not finite-sample guarantees or guarantees for increasing K.
- Severe class imbalance can produce material size inflation and undercoverage. A class size of 10 fails in the reported examples; no universal minimum class size is established.
- Passing `09_validate_revision.R` supports the numerical identities and data flow tested by that script. It does not prove statistical calibration for every population or validate clinical usefulness.
