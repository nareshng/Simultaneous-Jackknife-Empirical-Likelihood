# Simultaneous Jackknife Empirical Likelihood

R code accompanying the manuscript **“Simultaneous Jackknife Empirical Likelihood Inference for Comparing Multiple Classifiers in Ordered Three-Class Problems.”**

## Objective

When several fitted classifiers are evaluated on the same held-out subjects, their performance estimates are dependent. Ranking the observed values does not establish whether classifier differences exceed test-sample uncertainty, and performing all pairwise tests without adjustment inflates the familywise error rate.

This repository implements simultaneous inference for classifier-specific volumes under three-class ROC surfaces (VUSs). The method treats the vector of VUS estimators as a vector-valued three-sample U-statistic and provides:

- an omnibus JEL test of equal VUSs;
- pairwise JEL tests with Holm and Bonferroni adjustment;
- dependence-aware max-JEL inference using Gaussian multiplier calibration; and
- simultaneous profile-JEL confidence intervals for all pairwise VUS differences.

The inferential target is conditional on the fitted scoring rules. Training-sample and model-selection uncertainty are not included.

## Repository structure

| File | Purpose |
|---|---|
| `R/00_self_check.R` | Core pseudo-value and biomarker-data integrity checks |
| `R/01_jel_core.R` | VUS estimation, pseudo-values, JEL tests, max-JEL calibration, and confidence intervals |
| `R/02_models.R` | Base-R classifiers used in the real-data application |
| `R/03_real_data_application.R` | Cohort-separated urinary-biomarker analysis |
| `R/04_simulation_study.R` | Ten null and seven alternative simulation scenarios |
| `R/05_interval_simulation.R` | Simultaneous interval coverage and length experiments |
| `R/06_make_figures.R` | Regeneration of the four manuscript figures from result tables |
| `R/99_run_all.R` | Driver for real-data, diagnostic, and manuscript-scale workflows |

The eight R files are complementary; none is an obsolete duplicate.

## Software requirements

- R;
- the standard `parallel` package supplied with R; and
- no contributed R packages.

Run all commands from the repository root.

## Data

The real-data application uses the urinary-biomarker data from:

> Debernardi S, O'Brien H, Algahmdi AS, et al. (2020). A combination of urinary biomarker panel and PancRISK score for earlier detection of pancreatic cancer: A case-control study. *PLOS Medicine*, 17(12), e1003489. https://doi.org/10.1371/journal.pmed.1003489

Download the [official S1 Table](https://doi.org/10.1371/journal.pmed.1003489.s009) or the [Kaggle mirror](https://www.kaggle.com/datasets/johnjdavisiv/urinary-biomarkers-for-pancreatic-cancer), export it as a plain CSV, and save it as:

```text
data/urinary_biomarkers_pancreatic_cancer.csv
```

The required columns are:

```text
sample_id, patient_cohort, age, diagnosis,
creatinine, LYVE1, REG1B, TFF1
```

The integrity check expects 590 observations, cohort sizes 332 and 258, and cohort-2 class sizes 102 controls, 119 benign cases, and 37 PDAC cases. Cohort 1 is used for preprocessing, tuning, and model fitting; cohort 2 is the common held-out evaluation sample.

## Running the analysis

Clone the repository and install the data before running the checks:

```bash
git clone https://github.com/nareshng/Simultaneous-Jackknife-Empirical-Likelihood.git
cd Simultaneous-Jackknife-Empirical-Likelihood
Rscript R/00_self_check.R
```

Run only the real-data application:

```bash
Rscript R/99_run_all.R real
```

Run a reduced diagnostic analysis:

```bash
Rscript R/99_run_all.R quick
```

The diagnostic mode uses 200 replications, 400 multiplier draws, 50 interval replications, one core, and the manuscript-selected classifier hyperparameters. It tests the workflow but does not repeat hyperparameter selection at full scale.

Run the manuscript-scale analysis:

```bash
MC_CORES=8 Rscript R/99_run_all.R full
```

The full workflow uses 5,000 replications for each main simulation scenario, 2,000 multiplier draws per simulation replication, and 20,000 multiplier draws for the real-data analysis. It is computationally intensive. On Windows, the simulation runs serially because `parallel::mclapply()` is unavailable.

Individual scripts can also be run directly. For example:

```bash
R_REPS=200 B_MULT=400 MC_CORES=2 OUTPUT_SUFFIX=_check \
  Rscript R/04_simulation_study.R
```

## Configuration

| Variable | Default | Description |
|---|---:|---|
| `RETUNE` | `true` | Repeat five-fold training-cohort hyperparameter selection |
| `REAL_B_MULT` | `20000` | Multiplier draws for real-data max-JEL inference |
| `R_REPS` | `5000` | Replications per main simulation scenario |
| `INT_REPS` | scenario-specific | Override interval replications for all interval scenarios |
| `B_MULT` | `2000` | Multiplier draws per simulation replication |
| `MC_CORES` | `1` | Parallel workers on non-Windows systems |
| `ALPHA` | `0.05` | Nominal test level |
| `SIM_SEED` | `20260826` | Main simulation seed |
| `INT_SEED` | `20260827` | Interval simulation seed |
| `OUTPUT_SUFFIX` | `_R` | Suffix added to regenerated output files |

## Outputs

Scripts create `results/` and `figures/` as needed. Principal outputs include:

- scenario definitions and replicate-level simulation diagnostics;
- aggregate rejection probabilities, Monte Carlo standard errors, and numerical-failure counts;
- simultaneous interval coverage and mean length;
- held-out classifier scores and performance summaries;
- omnibus and pairwise JEL results; and
- PNG versions of the four manuscript figures.

Large replicate-level result files should normally be archived separately rather than committed to the main source tree. The exact aggregate tables, metadata, and final figures used in the manuscript should nevertheless be retained in a tagged release or archival deposit.

## Important implementation note

At present, `R/06_make_figures.R` expects `simulation_null_results*.csv` and `simulation_alternative_results*.csv`, whereas `R/04_simulation_study.R` writes only the long-format `simulation_aggregate_results*.csv`. Consequently, `quick` and `full` modes will fail at the figure-generation step unless the two wide-format tables are supplied or generated from the aggregate file. This mismatch should be corrected before claiming one-command reproducibility.

## Statistical interpretation

- All classifiers must be evaluated on the same independent held-out subjects.
- Cross-validation folds must not be treated as independent test samples.
- The method compares ordered three-class discrimination; it does not evaluate calibration, clinical utility, fairness, or transportability.
- The complete pairwise contrast vector has rank `K - 1`, so its covariance is structurally singular. The multiplier procedure avoids inverting this full covariance matrix.
- Severe class imbalance can invalidate asymptotic calibration. The manuscript's setting with only 10 observations in one class shows material size distortion.

## Citation

Please cite the accompanying simultaneous-inference manuscript when a public bibliographic record becomes available. The underlying multivariate three-sample JEL construction is described in:

> Garg N, Mathew L, Dewan I, Kattumannil SK. (2024; revised 2025). Jackknife empirical likelihood method for U statistics based on multivariate samples and its applications. arXiv:2408.14038v2. https://doi.org/10.48550/arXiv.2408.14038

## License

No software license is currently included. Add an explicit license before publication or archival; without one, reuse rights remain restricted by default copyright law.
