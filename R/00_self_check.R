# 00_self_check.R
# Fast checks of code, core arithmetic, supplied data and revision helpers.
# No Monte Carlo simulations, hyperparameter search, or manuscript estimates.

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

for (script in list.files(file.path(package_root, "R"), pattern = "[.]R$", full.names = TRUE)) {
  parse(file = script)
}
source(file.path(package_root, "R", "01_jel_core.R"))
source(file.path(package_root, "R", "02_models.R"))
jel_self_check()

# Direct, unequal-class, tie-aware enumeration checks the optimized VUS kernel.
small_scores <- list(
  cbind(c(0, 1, 1), c(3, 2, 1)),
  cbind(c(1, 1, 2, 3), c(2, 2, 1, 0)),
  cbind(c(1, 2, 2, 3, 4), c(3, 2, 1, 1, 0))
)
fit <- vus_and_pseudovalues(small_scores, ties = TRUE)
tie_kernel <- function(x, y, z) {
  if (x < y && y < z) 1 else if ((x == y && y < z) || (x < y && y == z)) {
    0.5
  } else if (x == y && y == z) 1 / 6 else 0
}
brute <- vapply(seq_len(ncol(small_scores[[1L]])), function(k) {
  total <- 0
  for (x in small_scores[[1L]][, k]) {
    for (y in small_scores[[2L]][, k]) {
      for (z in small_scores[[3L]][, k]) total <- total + tie_kernel(x, y, z)
    }
  }
  total / prod(vapply(small_scores, nrow, integer(1L)))
}, numeric(1L))
stopifnot(max(abs(fit$theta_hat - brute)) < 1e-12,
          max(abs(colMeans(fit$pseudo_values) - brute)) < 1e-12)

path <- file.path(package_root, "data", "urinary_biomarkers_pancreatic_cancer.csv")
if (!file.exists(path)) stop("Missing supplied biomarker CSV: ", path)
dat <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("sample_id", "patient_cohort", "diagnosis", "age", "creatinine", "LYVE1", "REG1B", "TFF1")
if (!all(required %in% names(dat))) stop("Missing CSV columns: ", paste(setdiff(required, names(dat)), collapse = ", "))
stopifnot(nrow(dat) == 590L, sum(dat$patient_cohort == 1) == 332L,
          sum(dat$patient_cohort == 2) == 258L,
          !anyNA(dat$sample_id), !anyDuplicated(dat$sample_id),
          all(nzchar(trimws(dat$sample_id))))
counts <- table(factor(dat$patient_cohort, levels = 1:2), factor(dat$diagnosis, levels = 1:3))
stopifnot(identical(as.integer(counts[1L, ]), c(81L, 89L, 162L)),
          identical(as.integer(counts[2L, ]), c(102L, 119L, 37L)))
predictors <- c("age", "creatinine", "LYVE1", "REG1B", "TFF1")
stopifnot(!anyNA(dat[, predictors]), all(is.finite(as.matrix(dat[, predictors]))),
          all(as.matrix(dat[, setdiff(predictors, "age")]) > 0))
train <- dat$patient_cohort == 1L
y <- as.integer(dat$diagnosis[train]) - 1L
RNGkind(kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
folds <- make_stratified_folds(y, v = 5L, seed = 2026L)
stopifnot(identical(folds, make_stratified_folds(y, v = 5L, seed = 2026L)),
          all(table(factor(y, levels = 0:2), factor(folds, levels = 1:5)) > 0))
prep <- fit_preprocessor(dat[train, predictors, drop = FALSE])
X <- apply_preprocessor(dat[train, predictors, drop = FALSE], prep)
stopifnot(max(abs(colMeans(X))) < 1e-10,
          max(abs(apply(X, 2L, stats::sd) - 1)) < 1e-10)

# 07 is sourced into its own environment under an explicit no-simulation guard.
helper_path <- file.path(package_root, "R", "07_revision_simulations.R")
if (!file.exists(helper_path)) stop("Missing revision helper script: ", helper_path)
old_guard <- Sys.getenv("REVISION_HELPERS_ONLY", unset = NA_character_)
Sys.setenv(REVISION_HELPERS_ONLY = "true")
helper_env <- new.env(parent = globalenv())
tryCatch(source(helper_path, local = helper_env), finally = {
  if (is.na(old_guard)) Sys.unsetenv("REVISION_HELPERS_ONLY") else {
    Sys.setenv(REVISION_HELPERS_ONLY = old_guard)
  }
})
if (exists("rev_helper_self_check", envir = helper_env, inherits = FALSE)) {
  helper_env$rev_helper_self_check()
} else stop("Revision helper self-check function is missing.")

results_dir <- file.path(package_root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)
suffix <- Sys.getenv("OUTPUT_SUFFIX", "_R")
report <- c("PASS: all R files parse.",
  "PASS: core self-check and direct unequal-class tie-aware VUS enumeration.",
  "PASS: supplied sample IDs, class/cohort counts, five-predictor completeness/positivity.",
  "PASS: deterministic stratified folds and training-only preprocessing arithmetic.",
  "PASS: revision helpers sourced without executing a simulation.",
  sprintf("DATA_MD5=%s", unname(tools::md5sum(path))),
  "Unique sample IDs do not prove distinct patients; no separate patient identifier is provided.",
  "This check does not reproduce manuscript tables or establish asymptotic calibration.",
  capture.output(sessionInfo()))
writeLines(report, file.path(results_dir, paste0("self_check", suffix, ".txt")))
cat(paste(report[1:5], collapse = "\n"), "\n")
