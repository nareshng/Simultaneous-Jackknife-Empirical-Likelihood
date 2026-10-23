# 03_real_data_application.R
# Urinary-biomarker application:
# cohort 1 is used only for preprocessing, tuning, and model fitting;
# cohort 2 is retained as a common independent evaluation sample.
#


.jel_find_root <- function() {
  explicit <- Sys.getenv("JEL_PROJECT_ROOT", unset = "")
  ofiles <- unlist(lapply(sys.frames(), function(frame) {
    value <- frame$ofile
    if (is.null(value)) character() else as.character(value)
  }), use.names = FALSE)
  cli <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  # Rscript represents spaces as ~+~ in the --file command-line argument.
  cli <- gsub("~+~", " ", cli, fixed = TRUE)
  starts <- c(explicit, dirname(rev(ofiles)), dirname(cli), getwd())
  for (start in starts[nzchar(starts)]) {
    candidate <- normalizePath(start, mustWork = FALSE)
    for (up in 0:2) {
      if (file.exists(file.path(candidate, "R", "01_jel_core.R"))) {
        return(normalizePath(candidate, mustWork = TRUE))
      }
      candidate <- dirname(candidate)
    }
  }
  stop("Cannot locate the project root. Set JEL_PROJECT_ROOT to its absolute path.")
}
package_root <- .jel_find_root()
source(file.path(package_root, "R", "01_jel_core.R"))
source(file.path(package_root, "R", "02_models.R"))

data_file <- file.path(package_root, "data", "urinary_biomarkers_pancreatic_cancer.csv")
results_dir <- file.path(package_root, "results")
figures_dir <- file.path(package_root, "figures")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

boolean_env <- function(name, default) {
  value <- tolower(Sys.getenv(name, default))
  if (!value %in% c("1", "true", "yes", "0", "false", "no")) {
    stop(name, " must be true or false.")
  }
  value %in% c("1", "true", "yes")
}
integer_env <- function(name, default, lower = 1L) {
  raw <- Sys.getenv(name, default)
  number <- suppressWarnings(as.numeric(raw))
  if (length(number) != 1L || !is.finite(number) ||
      number != floor(number) || number < lower || number > .Machine$integer.max) {
    stop(name, " must be an integer >= ", lower, ".")
  }
  as.integer(number)
}
retune <- boolean_env("RETUNE", "true")
Bmult <- integer_env("REAL_B_MULT", "20000")
suffix <- Sys.getenv("OUTPUT_SUFFIX", "_R")
seed <- integer_env("REAL_SEED", "2026", lower = 0L)
alpha <- suppressWarnings(as.numeric(Sys.getenv("ALPHA", "0.05")))
if (length(alpha) != 1L || !is.finite(alpha) || alpha <= 0 || alpha >= 1) {
  stop("ALPHA must lie strictly between zero and one.")
}
if (getRversion() < "3.6.0") stop("R >= 3.6.0 is required for the pinned sample RNG.")
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
if (!retune) {
  warning("RETUNE=false: fixed hyperparameters are used; no CV loss is recomputed. ",
          "Training_CV_log_loss will be NA, not a manuscript reference constant.")
}
writeLines(capture.output(sessionInfo()),
  file.path(results_dir, paste0("real_data_session_info", suffix, ".txt")))

dat <- utils::read.csv(data_file, stringsAsFactors = FALSE, check.names = FALSE)
required <- c(
  "sample_id", "patient_cohort", "age", "diagnosis",
  "creatinine", "LYVE1", "REG1B", "TFF1"
)
if (!all(required %in% names(dat))) {
  stop("Data file is missing: ", paste(setdiff(required, names(dat)), collapse = ", "))
}

predictors <- c("age", "creatinine", "LYVE1", "REG1B", "TFF1")
if (anyNA(dat[, predictors]) || any(!is.finite(as.matrix(dat[, predictors])))) {
  stop("The five prespecified predictors must be complete.")
}
if (any(dat[, c("creatinine", "LYVE1", "REG1B", "TFF1")] <= 0)) {
  stop("Biomarker variables to be log transformed must be strictly positive.")
}

if (anyNA(dat$sample_id) || any(!nzchar(trimws(dat$sample_id))) || anyDuplicated(dat$sample_id)) {
  stop("sample_id must be complete and unique across both cohorts.")
}
train <- dat$patient_cohort == 1
test <- dat$patient_cohort == 2
if (!all(train | test) || any(train & test)) stop("Unexpected cohort coding.")
if (length(intersect(dat$sample_id[train], dat$sample_id[test])) > 0L) {
  stop("Sample identifiers overlap between training and test cohorts.")
}

train_df <- dat[train, predictors, drop = FALSE]
test_df <- dat[test, predictors, drop = FALSE]
y_train <- as.integer(dat$diagnosis[train]) - 1L
y_test <- as.integer(dat$diagnosis[test]) - 1L
if (!identical(sort(unique(y_train)), 0:2) ||
    !identical(sort(unique(y_test)), 0:2)) {
  stop("Both cohorts must contain the three ordered diagnostic groups.")
}

methods <- c("ridge_multinomial", "regularized_qda", "rbf_krr", "weighted_knn")
display_names <- c(
  ridge_multinomial = "Ridge multinomial",
  regularized_qda = "Regularized QDA",
  rbf_krr = "RBF kernel ridge",
  weighted_knn = "Distance-weighted kNN"
)
fixed_hyper <- list(
  ridge_multinomial = list(lambda = 1),
  regularized_qda = list(lambda = 1),
  rbf_krr = list(gamma = 0.08, lambda = 0.1),
  weighted_knn = list(k = 25L)
)
if (retune) {
  folds <- make_stratified_folds(y_train, v = 5L, seed = seed)
  utils::write.csv(data.frame(sample_id = dat$sample_id[train],
    diagnosis = dat$diagnosis[train], fold = folds),
    file.path(results_dir, paste0("real_data_training_folds", suffix, ".csv")),
    row.names = FALSE)
}

fits <- list()
probabilities <- list()
cv_summaries <- list()
performance <- vector("list", length(methods))

for (m in methods) {
  message("Processing ", display_names[[m]])
  if (retune) {
    cv <- cross_validate_classifier(
      train_df, y_train, m,
      grid = default_parameter_grid(m),
      v = 5L, seed = seed, verbose = TRUE
    )
    hyper <- cv$best_hyperparameters
    cv_loss <- cv$best_cv_log_loss
    utils::write.csv(
      cv$table,
      file.path(results_dir, paste0("cv_", m, suffix, ".csv")),
      row.names = FALSE
    )
    cv_summaries[[m]] <- cv
  } else {
    hyper <- fixed_hyper[[m]]
    cv_loss <- NA_real_
  }

  fitted <- fit_classifier(train_df, y_train, m, hyper)
  prob <- predict_classifier(fitted, test_df)
  if (any(!is.finite(prob)) || max(abs(rowSums(prob) - 1)) > 1e-8) {
    stop("Invalid predicted probabilities for ", m)
  }
  score <- ordinal_score_from_probabilities(prob)

  fits[[m]] <- fitted
  probabilities[[m]] <- prob
  performance[[which(methods == m)]] <- data.frame(
    Classifier = display_names[[m]],
    Selected_hyperparameters = paste(
      paste(names(hyper), unlist(hyper), sep = "="), collapse = "; "
    ),
    Training_CV_log_loss = cv_loss,
    Hyperparameter_source = if (retune) "five-fold cohort-1 CV" else "fixed; CV not run",
    Heldout_accuracy = classification_accuracy(prob, y_test),
    Heldout_log_loss = multiclass_log_loss(prob, y_test),
    stringsAsFactors = FALSE
  )
}

score_matrix <- do.call(cbind, lapply(probabilities, ordinal_score_from_probabilities))
colnames(score_matrix) <- unname(display_names[methods])
scores_by_class <- lapply(0:2, function(cls) score_matrix[y_test == cls, , drop = FALSE])
jel <- vus_and_pseudovalues(scores_by_class, ties = TRUE)
global <- jel_global_test(jel$theta_hat, jel$pseudo_values, alpha = alpha)
global_wald <- wald_global_test(jel$theta_hat, jel$pseudo_values, alpha = alpha, ns = jel$ns)
pair <- jel_pairwise_tests(jel$theta_hat, jel$pseudo_values)
maxcal <- max_jel_calibration(
  jel$theta_hat, jel$pseudo_values, jel$ns,
  alpha = alpha, B = Bmult, seed = seed
)
ci <- jel_simultaneous_intervals(
  jel$theta_hat, jel$pseudo_values, jel$ns,
  alpha = alpha, B = Bmult, seed = seed, method = "maxJEL"
)

perf <- do.call(rbind, performance)
perf$Heldout_VUS <- jel$theta_hat
utils::write.csv(
  perf,
  file.path(results_dir, paste0("real_data_model_performance", suffix, ".csv")),
  row.names = FALSE
)

score_output <- data.frame(
  sample_id = dat$sample_id[test],
  diagnosis = dat$diagnosis[test],
  ordered_class = y_test,
  score_matrix,
  check.names = FALSE
)
utils::write.csv(
  score_output,
  file.path(results_dir, paste0("real_data_heldout_scores", suffix, ".csv")),
  row.names = FALSE
)

prob_output <- data.frame(sample_id = dat$sample_id[test], ordered_class = y_test)
for (m in methods) {
  for (cls in 0:2) prob_output[[paste0(m, "_p", cls)]] <- probabilities[[m]][, cls + 1L]
}
utils::write.csv(prob_output,
  file.path(results_dir, paste0("real_data_heldout_probabilities", suffix, ".csv")),
  row.names = FALSE)
saveRDS(list(fits = fits, cross_validation = cv_summaries,
             training_sample_ids = dat$sample_id[train], seed = seed,
             alpha = alpha, multiplier_draws = Bmult, retuned = retune),
  file.path(results_dir, paste0("real_data_fitted_models", suffix, ".rds")))

pair$Classifier_A <- unname(display_names[methods[pair$classifier_a]])
pair$Classifier_B <- unname(display_names[methods[pair$classifier_b]])
pair$max_JEL_adjusted_p <- maxcal$adjusted_p
pair$max_JEL_reject <- maxcal$reject
pair$max_JEL_CI_lower <- ci$lower
pair$max_JEL_CI_upper <- ci$upper
pair_out <- pair[, c(
  "Classifier_A", "Classifier_B", "difference", "statistic",
  "p_value", "holm_p", "bonferroni_p", "max_JEL_adjusted_p",
  "max_JEL_reject", "max_JEL_CI_lower", "max_JEL_CI_upper"
)]
names(pair_out)[3:7] <- c(
  "VUS_difference_A_minus_B", "JEL_statistic", "Raw_JEL_p",
  "Holm_adjusted_p", "Bonferroni_adjusted_p"
)
utils::write.csv(
  pair_out,
  file.path(results_dir, paste0("real_data_pairwise_inference", suffix, ".csv")),
  row.names = FALSE
)

class_counts <- as.list(table(factor(y_test, levels = 0:2)))
names(class_counts) <- c("control", "benign", "pdac")
global_lines <- c(
  sprintf("Training cohort size: %d", sum(train)),
  sprintf("Held-out cohort size: %d", sum(test)),
  sprintf("Held-out class sizes: control=%d, benign=%d, PDAC=%d",
          class_counts$control, class_counts$benign, class_counts$pdac),
  sprintf("Global JEL statistic: %.10f", global$statistic),
  sprintf("Global JEL df: %d", global$df),
  sprintf("Global JEL p-value: %.10g", global$p_value),
  sprintf("Global Wald statistic: %.10f", global_wald$statistic),
  sprintf("Global Wald df: %d", global_wald$df),
  sprintf("Global Wald p-value: %.10g", global_wald$p_value),
  sprintf("Training-only CV rerun: %s", retune),
  sprintf("Application seed: %d", seed),
  sprintf("Nominal alpha: %.10g", alpha),
  sprintf("Max-JEL multiplier draws: %d", Bmult),
  sprintf("Max-JEL critical value: %.10f", unique(maxcal$critical_value)),
  "Inference target: conditional on the fitted classifiers; held-out evaluation distribution.",
  "Unique sample IDs do not establish person-level independence."
)
writeLines(
  global_lines,
  file.path(results_dir, paste0("real_data_global_test", suffix, ".txt"))
)

writeLines(c(global_lines,
  sprintf("RUN_ID=%s", Sys.getenv("JEL_RUN_ID", "direct-script")),
  sprintf("DATA_MD5=%s", unname(tools::md5sum(data_file))),
  "RNGkind=Mersenne-Twister; Inversion; Rejection",
  "CV criterion=unweighted mean of five fold-average clipped multiclass log losses",
  "RETUNE=false produces NA CV losses; reference constants are never reported.",
  capture.output(sessionInfo())),
  file.path(results_dir, paste0("real_data_run_metadata", suffix, ".txt")))

# Data-audit table documents missingness and the strict cohort separation.
audit <- data.frame(
  variable = names(dat),
  missing = vapply(dat, function(x) sum(is.na(x) | trimws(as.character(x)) == ""), integer(1L)),
  stringsAsFactors = FALSE
)
utils::write.csv(
  audit,
  file.path(results_dir, paste0("real_data_missingness_audit", suffix, ".csv")),
  row.names = FALSE
)

# Figure: held-out VUS estimates.
png(
  file.path(figures_dir, paste0("fig3_real_vus", suffix, ".png")),
  width = 1800, height = 1050, res = 180
)
op <- par(mar = c(8, 5, 3, 1))
plot(
  seq_along(jel$theta_hat), jel$theta_hat,
  ylim = range(c(1 / 6, jel$theta_hat)) + c(-0.03, 0.03),
  xaxt = "n", xlab = "", ylab = "Held-out VUS",
  pch = 19, cex = 1.4,
  main = "Classifier VUS on independent cohort 2"
)
axis(1, at = seq_along(jel$theta_hat), labels = colnames(score_matrix), las = 2)
abline(h = 1 / 6, lty = 2)
mtext("Chance ordering = 1/6", side = 4, line = -1.5, at = 1 / 6, cex = 0.8)
par(op)
dev.off()

# Figure: simultaneous max-JEL intervals.
png(
  file.path(figures_dir, paste0("fig4_pairwise_intervals", suffix, ".png")),
  width = 1900, height = 1250, res = 180
)
op <- par(mar = c(5, 16, 3, 1))
labels <- paste(pair_out$Classifier_A, "-", pair_out$Classifier_B)
ord <- rev(seq_len(nrow(pair_out)))
xlim <- range(c(pair_out$max_JEL_CI_lower, pair_out$max_JEL_CI_upper, 0))
plot(
  pair_out$VUS_difference_A_minus_B[ord], seq_along(ord),
  xlim = xlim + c(-0.02, 0.02), ylim = c(0.5, length(ord) + 0.5),
  yaxt = "n", ylab = "", xlab = "VUS difference (A - B)",
  pch = 19, main = sprintf("%.1f%% simultaneous max-JEL intervals", 100 * (1 - alpha))
)
segments(
  pair_out$max_JEL_CI_lower[ord], seq_along(ord),
  pair_out$max_JEL_CI_upper[ord], seq_along(ord),
  lwd = 2
)
axis(2, at = seq_along(ord), labels = labels[ord], las = 2, cex.axis = 0.75)
abline(v = 0, lty = 2)
par(op)
dev.off()

message("Real-data analysis completed. Outputs written to ", results_dir)

