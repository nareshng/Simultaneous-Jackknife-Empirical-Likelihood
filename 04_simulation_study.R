# 04_simulation_study.R
# Manuscript-grade stress simulations for global and simultaneous JEL inference.
#
# Full manuscript run (computationally intensive):
#   R_REPS=5000 B_MULT=1999 MC_CORES=8 Rscript R/04_simulation_study.R
#
# Quick diagnostic run:
#   R_REPS=200 B_MULT=399 MC_CORES=2 Rscript R/04_simulation_study.R
#
# The script saves replicate-level diagnostics, aggregate estimates, Monte Carlo
# standard errors, numerical-failure counts, and the full scenario registry.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
if (is.na(script_file) || !nzchar(script_file)) {
  package_root <- normalizePath(getwd())
} else {
  package_root <- normalizePath(file.path(dirname(script_file), ".."))
}
source(file.path(package_root, "R", "01_jel_core.R"))

results_dir <- file.path(package_root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

Rreps <- as.integer(Sys.getenv("R_REPS", "5000"))
Bmult <- as.integer(Sys.getenv("B_MULT", "1999"))
cores <- as.integer(Sys.getenv("MC_CORES", "1"))
alpha <- as.numeric(Sys.getenv("ALPHA", "0.05"))
master_seed <- as.integer(Sys.getenv("SIM_SEED", "20260826"))
suffix <- Sys.getenv("OUTPUT_SUFFIX", "_R")

scenario <- function(id, description, K, ns, rho, distribution, d, alternative = FALSE) {
  list(id = id, description = description, K = K, ns = as.integer(ns),
       rho = rho, distribution = distribution, d = as.numeric(d),
       alternative = alternative)
}

scenarios <- list(
  scenario("N1", "Gaussian, small balanced", 4, c(20,20,20), .5, "gaussian", rep(.80,4)),
  scenario("N2", "Gaussian, low dependence", 4, c(40,40,40), .2, "gaussian", rep(.80,4)),
  scenario("N3", "Gaussian, high dependence", 4, c(40,40,40), .8, "gaussian", rep(.80,4)),
  scenario("N4", "Gaussian, cohort-matched imbalance", 4, c(102,119,37), .5, "gaussian", rep(.80,4)),
  scenario("N5", "Gaussian, severe imbalance", 4, c(10,40,160), .5, "gaussian", rep(.80,4)),
  scenario("N6", "Heavy-tailed t(3)", 4, c(40,40,40), .5, "t3", rep(.80,4)),
  scenario("N7", "Log-normal margins", 4, c(40,40,40), .5, "lognormal", rep(.80,4)),
  scenario("N8", "10% contaminated normal", 4, c(40,40,40), .5, "contaminated", rep(.80,4)),
  scenario("N9", "Rounded scores with ties", 4, c(40,40,40), .5, "rounded", rep(.80,4)),
  scenario("N10", "Six classifiers", 6, c(60,60,60), .5, "gaussian", rep(.80,6)),
  scenario("A1", "Dense Gaussian alternative", 4, c(40,40,40), .5, "gaussian", c(.65,.75,.85,.95), TRUE),
  scenario("A2", "Dense alternative, high dependence", 4, c(40,40,40), .8, "gaussian", c(.65,.75,.85,.95), TRUE),
  scenario("A3", "Sparse Gaussian alternative", 4, c(40,40,40), .5, "gaussian", c(.80,.80,.80,1.00), TRUE),
  scenario("A4", "Two-cluster Gaussian alternative", 4, c(40,40,40), .5, "gaussian", c(.72,.72,.92,.92), TRUE),
  scenario("A5", "Dense heavy-tailed alternative", 4, c(40,40,40), .5, "t3", c(.65,.75,.85,.95), TRUE),
  scenario("A6", "Sparse cohort-matched alternative", 4, c(102,119,37), .5, "gaussian", c(.80,.80,.80,1.00), TRUE),
  scenario("A7", "Sparse rounded-score alternative", 4, c(40,40,40), .5, "rounded", c(.80,.80,.80,1.00), TRUE)
)

registry <- do.call(rbind, lapply(scenarios, function(s) {
  data.frame(
    ScenarioID = s$id, Scenario = s$description, K = s$K,
    n1 = s$ns[1], n2 = s$ns[2], n3 = s$ns[3],
    rho = s$rho, Distribution = s$distribution,
    d = paste(format(s$d, nsmall = 2), collapse = ","),
    Alternative = s$alternative,
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(
  registry,
  file.path(results_dir, paste0("simulation_scenario_registry", suffix, ".csv")),
  row.names = FALSE
)

generate_correlated <- function(n, K, rho) {
  Sigma <- matrix(rho, K, K)
  diag(Sigma) <- 1
  matrix(stats::rnorm(n * K), n, K) %*% chol(Sigma)
}

generate_scores <- function(s) {
  out <- vector("list", 3L)
  for (cc in 0:2) {
    nclass <- s$ns[cc + 1L]
    base <- generate_correlated(nclass, s$K, s$rho)

    if (s$distribution == "gaussian" || s$distribution == "rounded") {
      E <- base
      score <- sweep(E, 2L, cc * s$d, "+")
      if (s$distribution == "rounded") score <- round(score / 0.5) * 0.5
    } else if (s$distribution == "t3") {
      E <- base / sqrt(stats::rchisq(nclass, df = 3) / 3)
      score <- sweep(E, 2L, cc * s$d, "+")
    } else if (s$distribution == "lognormal") {
      score <- exp(sweep(base, 2L, cc * s$d, "+"))
    } else if (s$distribution == "contaminated") {
      E <- base
      bad <- stats::runif(nclass) < 0.10
      if (any(bad)) {
        E[bad, ] <- 5 * generate_correlated(sum(bad), s$K, s$rho)
      }
      score <- sweep(E, 2L, cc * s$d, "+")
    } else {
      stop("Unknown distribution: ", s$distribution)
    }
    out[[cc + 1L]] <- score
  }
  out
}

evaluate_replicate <- function(scores, s, replicate_seed) {
  K <- s$K
  pm <- pair_contrast_matrix(K)
  true_null <- abs(outer(s$d, s$d, "-")[lower.tri(matrix(0, K, K))]) < 1e-12
  # combn ordering differs from lower.tri extraction; compute directly.
  true_null <- vapply(seq_len(nrow(pm$pairs)), function(j) {
    abs(s$d[pm$pairs[j, 1L]] - s$d[pm$pairs[j, 2L]]) < 1e-12
  }, logical(1L))
  false_null <- !true_null

  ans <- tryCatch({
    jel <- vus_and_pseudovalues(scores, ties = TRUE)
    global <- jel_global_test(jel$theta_hat, jel$pseudo_values, alpha)
    global_wald <- wald_global_test(jel$theta_hat, jel$pseudo_values, alpha)
    pair <- jel_pairwise_tests(jel$theta_hat, jel$pseudo_values)
    maxcal <- max_jel_calibration(
      jel$theta_hat, jel$pseudo_values, jel$ns,
      alpha = alpha, B = Bmult, seed = replicate_seed
    )

    P <- pm$P
    delta <- as.vector(P %*% jel$theta_hat)
    a <- pooled_centering_vector(jel$ns)
    G <- jel$pseudo_values %*% t(P) - a %o% delta
    sdev <- sqrt(colMeans(G^2))
    wald_t <- sqrt(nrow(G)) * delta / sdev

    set.seed(replicate_seed)
    Xi <- matrix(stats::rnorm(Bmult * nrow(G)), Bmult, nrow(G))
    Zstar <- (Xi %*% G) / sqrt(nrow(G))
    Zstar <- sweep(Zstar, 2L, sdev, "/")
    max_star <- apply(abs(Zstar), 1L, max)
    wald_critical <- as.numeric(stats::quantile(max_star, 1 - alpha, type = 1L))
    wald_reject <- abs(wald_t) > wald_critical

    raw_reject <- pair$p_value < alpha
    holm_reject <- pair$holm_p < alpha
    max_reject <- maxcal$reject

    c(
      numerical_failure = 0,
      global_jel_reject = as.numeric(global$reject),
      global_wald_reject = as.numeric(global_wald$reject),
      unadjusted_any = as.numeric(any(raw_reject)),
      holm_any = as.numeric(any(holm_reject)),
      maxjel_any = as.numeric(any(max_reject)),
      waldmax_any = as.numeric(any(wald_reject)),
      holm_any_true_discovery = if (any(false_null)) as.numeric(any(holm_reject[false_null])) else NA_real_,
      maxjel_any_true_discovery = if (any(false_null)) as.numeric(any(max_reject[false_null])) else NA_real_,
      waldmax_any_true_discovery = if (any(false_null)) as.numeric(any(wald_reject[false_null])) else NA_real_,
      holm_true_null_fwer = if (any(true_null)) as.numeric(any(holm_reject[true_null])) else NA_real_,
      maxjel_true_null_fwer = if (any(true_null)) as.numeric(any(max_reject[true_null])) else NA_real_,
      waldmax_true_null_fwer = if (any(true_null)) as.numeric(any(wald_reject[true_null])) else NA_real_,
      global_jel_statistic = global$statistic,
      global_wald_statistic = global_wald$statistic,
      maxjel_critical = unique(maxcal$critical_value),
      waldmax_critical = wald_critical
    )
  }, error = function(e) {
    c(
      numerical_failure = 1,
      global_jel_reject = NA, global_wald_reject = NA,
      unadjusted_any = NA, holm_any = NA, maxjel_any = NA, waldmax_any = NA,
      holm_any_true_discovery = NA, maxjel_any_true_discovery = NA,
      waldmax_any_true_discovery = NA,
      holm_true_null_fwer = NA, maxjel_true_null_fwer = NA,
      waldmax_true_null_fwer = NA,
      global_jel_statistic = NA, global_wald_statistic = NA,
      maxjel_critical = NA, waldmax_critical = NA
    )
  })
  ans
}

mean_mcse <- function(x) {
  valid <- is.finite(x)
  n <- sum(valid)
  if (n == 0L) return(c(estimate = NA_real_, mcse = NA_real_, n_valid = 0))
  p <- mean(x[valid])
  c(estimate = p, mcse = sqrt(p * (1 - p) / n), n_valid = n)
}

run_scenario <- function(s, scenario_index) {
  message("Starting ", s$id, ": ", s$description)
  seeds <- master_seed + scenario_index * 1000000L + seq_len(Rreps)

  one <- function(rr) {
    set.seed(seeds[rr])
    scores <- generate_scores(s)
    c(replicate = rr, seed = seeds[rr],
      evaluate_replicate(scores, s, replicate_seed = seeds[rr] + 500000L))
  }

  if (cores > 1L && .Platform$OS.type != "windows") {
    values <- parallel::mclapply(seq_len(Rreps), one, mc.cores = cores,
                                mc.preschedule = TRUE)
  } else {
    values <- lapply(seq_len(Rreps), one)
  }
  M <- as.data.frame(do.call(rbind, values))
  M$ScenarioID <- s$id
  M$Scenario <- s$description
  M <- M[, c("ScenarioID", "Scenario", setdiff(names(M), c("ScenarioID", "Scenario")))]

  utils::write.csv(
    M,
    file.path(results_dir, paste0("replicates_", s$id, suffix, ".csv")),
    row.names = FALSE
  )

  binary_metrics <- c(
    "global_jel_reject", "global_wald_reject", "unadjusted_any",
    "holm_any", "maxjel_any", "waldmax_any",
    "holm_any_true_discovery", "maxjel_any_true_discovery",
    "waldmax_any_true_discovery", "holm_true_null_fwer",
    "maxjel_true_null_fwer", "waldmax_true_null_fwer"
  )
  agg <- do.call(rbind, lapply(binary_metrics, function(v) {
    z <- mean_mcse(M[[v]])
    data.frame(
      ScenarioID = s$id, Scenario = s$description, Metric = v,
      Estimate = unname(z["estimate"]), MCSE = unname(z["mcse"]),
      Valid_replicates = as.integer(z["n_valid"]),
      Total_replicates = Rreps,
      Numerical_failures = sum(M$numerical_failure == 1, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  message("Completed ", s$id)
  agg
}

aggregate_list <- vector("list", length(scenarios))
for (ss in seq_along(scenarios)) {
  aggregate_list[[ss]] <- run_scenario(scenarios[[ss]], ss)
}
aggregate <- do.call(rbind, aggregate_list)
utils::write.csv(
  aggregate,
  file.path(results_dir, paste0("simulation_aggregate_results", suffix, ".csv")),
  row.names = FALSE
)

metadata <- c(
  sprintf("R_REPS=%d", Rreps),
  sprintf("B_MULT=%d", Bmult),
  sprintf("MC_CORES=%d", cores),
  sprintf("ALPHA=%.8f", alpha),
  sprintf("SIM_SEED=%d", master_seed),
  sprintf("R_VERSION=%s", R.version.string),
  sprintf("PLATFORM=%s", R.version$platform),
  sprintf("TIMESTAMP_UTC=%s", format(Sys.time(), tz = "UTC", usetz = TRUE))
)
writeLines(
  metadata,
  file.path(results_dir, paste0("simulation_run_metadata", suffix, ".txt"))
)
message("All simulation scenarios completed.")
