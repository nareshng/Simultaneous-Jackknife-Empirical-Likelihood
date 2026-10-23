# 09_validate_revision.R
# Independent, deterministic numerical validation using base R only.
# Run from the extracted package: Rscript R/09_validate_revision.R
# This script does not fit classifiers or launch the simulation study.

validate_revision <- function() {
  argv <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  script <- if (length(argv)) sub("^--file=", "", argv[1L]) else ""
  source_files <- vapply(sys.frames(), function(frame) {
    if (exists("ofile", envir = frame, inherits = FALSE)) as.character(get("ofile", envir = frame))[1L] else ""
  }, character(1L))
  source_files <- source_files[nzchar(source_files)]
  explicit_root <- Sys.getenv("JEL_PROJECT_ROOT", "")
  root <- if (nzchar(explicit_root)) normalizePath(explicit_root, mustWork = TRUE) else if (length(source_files)) {
    normalizePath(file.path(dirname(tail(source_files, 1L)), ".."), mustWork = TRUE)
  } else if (nzchar(script)) normalizePath(file.path(dirname(script), ".."), mustWork = TRUE) else normalizePath(getwd(), mustWork = TRUE)
  if (!file.exists(file.path(root, "R", "01_jel_core.R"))) stop("Run from the extracted R-only package root.")
  prior <- Sys.getenv("REVISION_HELPERS_ONLY", unset = NA_character_)
  Sys.setenv(REVISION_HELPERS_ONLY = "true")
  on.exit(if (is.na(prior)) Sys.unsetenv("REVISION_HELPERS_ONLY") else Sys.setenv(REVISION_HELPERS_ONLY = prior), add = TRUE)
  env <- new.env(parent = globalenv())
  source(file.path(root, "R", "01_jel_core.R"), local = env)
  source(file.path(root, "R", "07_revision_simulations.R"), local = env)
  log <- c("Independent numerical validation of the R-only revision", R.version.string,
           "No classifier training or full simulation was run by this validation script.")
  pass <- function(message) { line <- paste0("PASS: ", message); cat(line, "\n"); log <<- c(log, line) }
  close <- function(actual, expected, tolerance = 1e-8, label = "values") {
    if (length(actual) != length(expected) || any(!is.finite(actual)) || any(!is.finite(expected)) ||
        max(abs(as.numeric(actual) - as.numeric(expected))) > tolerance) {
      stop(label, " disagree; maximum absolute error = ", suppressWarnings(max(abs(actual - expected))))
    }
    invisible(TRUE)
  }
  # Enumerate the defining kernel independently of the fast pairwise-array code.
  brute <- function(scores, ties = TRUE, return_H = FALSE) {
    sizes <- vapply(scores, nrow, integer(1L)); K <- ncol(scores[[1L]])
    index <- expand.grid(seq_len(sizes[1L]), seq_len(sizes[2L]), seq_len(sizes[3L]))
    h <- matrix(0, nrow(index), K)
    for (r in seq_len(nrow(index))) {
      x <- scores[[1L]][index[r, 1L], ]; y <- scores[[2L]][index[r, 2L], ]; z <- scores[[3L]][index[r, 3L], ]
      h[r, ] <- as.numeric(x < y & y < z)
      if (ties) h[r, ] <- h[r, ] + .5 * as.numeric((x == y & y < z) | (x < y & y == z)) + as.numeric(x == y & y == z) / 6
    }
    H <- if (return_H) do.call(rbind, lapply(seq_len(3L), function(s) {
      do.call(rbind, lapply(seq_len(sizes[s]), function(i) colMeans(h[index[, s] == i, , drop = FALSE])))
    })) else NULL
    list(theta = colMeans(h), sum = colSums(h), H = H)
  }
  # A one-dimensional EL oracle uses an independent uniroot solve.
  scalar_oracle <- function(g) {
    g <- as.numeric(g)
    if (all(g == 0)) return(0)
    if (min(g) >= 0 || max(g) <= 0) return(Inf)
    lo <- (-1 / max(g)) * (1 - 1e-11); hi <- (-1 / min(g)) * (1 - 1e-11)
    lambda <- uniroot(function(t) sum(g / (1 + t * g)), c(lo, hi), tol = 1e-12)$root
    max(0, 2 * sum(log1p(lambda * g)))
  }
  set.seed(20261018L)
  sizes <- c(4L, 5L, 7L)
  scores <- lapply(sizes, function(m) matrix(sample(-2:3, m * 3L, replace = TRUE), m, 3L))
  reference <- brute(scores, return_H = TRUE)
  fit <- env$vus_and_pseudovalues(scores, ties = TRUE)
  close(fit$theta_hat, reference$theta, label = "triplet VUS")
  n <- sum(sizes); offset <- 0L
  for (s in seq_len(3L)) {
    for (i in seq_len(sizes[s])) {
      smaller <- scores; smaller[[s]] <- smaller[[s]][-i, , drop = FALSE]
      remaining_sum <- brute(smaller)$sum
      # Original fixed transformed-kernel constant; no estimator helper reused.
      deleted_pooled <- choose(n, 3L) / prod(sizes) * remaining_sum / choose(n - 1L, 3L)
      expected <- n * reference$theta - (n - 1L) * deleted_pooled
      close(fit$pseudo_values[offset + i, ], expected, label = "fixed-kernel pooled deletion")
    }
    offset <- offset + sizes[s]
  }
  close(env$vus_and_pseudovalues(scores, ties = FALSE)$theta_hat, brute(scores, ties = FALSE)$theta, label = "strict-order kernel")
  a <- env$pooled_centering_vector(sizes)
  close(mean(a), 1); close(colMeans(fit$pseudo_values), fit$theta_hat)
  pass("fast tie-aware and strict VUS agree with full triplet enumeration; every unequal-class pooled delete-one identity holds")

  all_tied <- lapply(c(3L, 4L, 5L), function(m) matrix(0, m, 2L))
  tied <- env$vus_and_pseudovalues(all_tied)
  tied_a <- env$pooled_centering_vector(tied$ns)
  close(tied$theta_hat, rep(1 / 6, 2L))
  close(tied$pseudo_values, tied_a %o% rep(1 / 6, 2L))
  cubic <- env$vus_and_pseudovalues(lapply(scores, function(x) x^3))
  close(cubic$theta_hat, fit$theta_hat); close(cubic$pseudo_values, fit$pseudo_values)
  perm <- lapply(scores, function(x) x[rev(seq_len(nrow(x))), , drop = FALSE])
  perm_fit <- env$vus_and_pseudovalues(perm)
  rows <- unlist(lapply(seq_len(3L), function(s) sum(sizes[seq_len(s - 1L)]) + rev(seq_len(sizes[s]))))
  close(perm_fit$theta_hat, fit$theta_hat); close(perm_fit$pseudo_values, fit$pseudo_values[rows, ])
  pass("all-tied chance ordering is 1/6; monotone score transforms and subject permutations preserve the defining inference inputs")

  g <- c(rep(-1, 4L), rep(2, 6L))
  # Exact category-constrained probabilities: minus rows 1/6, plus rows 1/18.
  scalar_exact <- -2 * (4 * log(10 / 6) + 6 * log(10 / 18))
  close(scalar_oracle(g), scalar_exact)
  close(env$el_mean_statistic(matrix(g, ncol = 1L)), scalar_exact, tolerance = 1e-7)
  close(env$rev_scalar_jel(g), scalar_exact, tolerance = 1e-7)
  for (scale in c(.01, 7, -3)) {
    close(env$rev_scalar_jel(scale * g), scalar_exact, tolerance = 1e-7)
    close(env$el_mean_statistic(matrix(scale * g, ncol = 1L)), scalar_exact, tolerance = 1e-7)
  }
  stopifnot(is.infinite(env$rev_scalar_jel(c(0, 1, 2))), env$rev_scalar_jel(rep(0, 8)) == 0)
  pass("scalar solvers recover exact constrained EL, are invariant to nonzero rescaling, and the revision solver distinguishes zero from infeasible scores")

  G <- rbind(matrix(c(1, 0), 3L, 2L, byrow = TRUE), c(-1, 0),
             matrix(c(0, 1), 2L, 2L, byrow = TRUE), matrix(c(0, -1), 2L, 2L, byrow = TRUE))
  vector_exact <- 2 * (3 * log(1.5) + log(.5))
  transformed <- G %*% matrix(c(1.2, .7, -.4, 1.5), 2L, 2L)
  for (x in list(G, transformed)) {
    baseline <- env$el_mean_statistic(x, return_details = TRUE)
    revised <- env$rev_global_jel(x)
    stopifnot(baseline$converged, revised$converged)
    close(baseline$statistic, vector_exact, 1e-7); close(revised$statistic, vector_exact, 1e-7)
  }
  bad <- matrix(NA_real_, 12L, 3L)
  failure <- env$jel_global_test(rep(0, 3L), bad)
  stopifnot(is.na(failure$reject), !failure$converged, !is.finite(failure$statistic))
  infeasible_vector <- rbind(c(1, 0), c(0, 1), c(1, 1), c(2, 1), c(1, 2))
  no_interior_solution <- env$rev_global_jel(infeasible_vector)
  stopifnot(!no_interior_solution$converged, !is.finite(no_interior_solution$statistic))
  pass("global EL agrees with an exact two-coordinate solution and an invertible change of coordinates; nonfinite solves and an excluded-origin vector hull are not reported as converged decisions")

  C <- env$global_contrast_matrix(3L)
  delta <- as.vector(C %*% reference$theta)
  b <- rep(n * (n - 1L) / (sizes * (n - 3L)), sizes)
  independent_influence <- b * sweep(reference$H, 2L, reference$theta, "-")
  influence_C <- independent_influence %*% t(C)
  covariance <- crossprod(influence_C) / n
  wald_expected <- n * sum(delta * solve(covariance, delta))
  wald <- env$wald_global_test(fit$theta_hat, fit$pseudo_values, ns = sizes)
  stopifnot(wald$full_rank, identical(wald$centering, "class-factor a_s"))
  close(wald$covariance, covariance); close(wald$statistic, wald_expected)
  legacy <- suppressWarnings(env$wald_global_test(fit$theta_hat, fit$pseudo_values))
  stopifnot(max(abs(legacy$covariance - covariance)) > 1e-6)
  kappa <- mean((a - 1)^2)
  close(legacy$covariance, covariance + kappa * delta %o% delta)
  close(legacy$statistic, wald$statistic / (1 + kappa * wald$statistic / n))
  pass("corrected global Wald matches enumerated influences; the legacy covariance differs by its exact positive rank-one term and obeys the Sherman-Morrison statistic identity")

  # Xi=sqrt(n)*I gives an exact empirical second-moment quadrature for G.
  pm <- env$pair_contrast_matrix(3L)
  pair_influence <- independent_influence %*% t(pm$P)
  pair_sd <- sqrt(colMeans(pair_influence^2))
  maxima <- apply(abs(sweep(pair_influence, 2L, pair_sd, "/")), 1L, max)
  calibration <- env$max_jel_calibration(fit$theta_hat, fit$pseudo_values, sizes,
                                        B = n, multiplier_matrix = sqrt(n) * diag(n))
  close(unique(calibration$critical_value), unname(quantile(maxima, .95, type = 1L)))
  # For K-dimensional revision sampling, +/-sqrt(K)e_j has exact mean zero
  # and identity raw second moment. This validates the conditional Gaussian
  # covariance factorization without relying on any particular eigenvector sign.
  K <- 3L; Z <- rbind(sqrt(K) * diag(K), -sqrt(K) * diag(K))
  revised_cal <- env$rev_max_calibration(fit$theta_hat, fit$pseudo_values, a,
                                        B = 2L * K, standard_normals = Z)
  close(revised_cal$covariance, crossprod(independent_influence) / n)
  close(crossprod(revised_cal$gaussian_draws) / nrow(Z), revised_cal$covariance)
  expected_var <- diag(crossprod(pair_influence) / n)
  close(revised_cal$variance, expected_var)
  pass("subject-multiplier maximum and reduced-dimensional Gaussian factorization recover independently enumerated influence covariance")

  # Independent profile inversion for a regular scalar estimating equation.
  v <- rep(c(-.30, -.10, .20, .50), c(12L, 18L, 17L, 13L))
  estimate <- mean(v); critical <- qnorm(.975)
  profile <- env$jel_profile_interval(v, rep(1, length(v)), estimate, critical)
  target <- critical^2
  oracle_f <- function(d) scalar_oracle(v - d) - target
  lo <- uniroot(oracle_f, c(min(v) + 1e-7, estimate), tol = 1e-10)$root
  hi <- uniroot(oracle_f, c(estimate, max(v) - 1e-7), tol = 1e-10)$root
  close(profile, c(lo, hi), tolerance = 2e-6)
  close(c(scalar_oracle(v - profile[1L]), scalar_oracle(v - profile[2L])), rep(target, 2L), tolerance = 1e-4)
  stopifnot(profile[1L] < estimate, profile[2L] > estimate)
  tiny <- env$jel_profile_interval(v * 1e-4, rep(1, length(v)), estimate * 1e-4, critical)
  close(tiny, profile * 1e-4, tolerance = 1e-10)
  stopifnot(tiny[2L] > tiny[1L])
  symmetric_tiny <- env$jel_profile_interval(rep(c(-1e-4, 1e-4), 20L), rep(1, 40L), 0, critical)
  stopifnot(symmetric_tiny[1L] < 0, symmetric_tiny[2L] > 0)
  close(scalar_oracle(rep(c(-1e-4, 1e-4), 20L) - symmetric_tiny[1L]), target, tolerance = 1e-5)
  pass("regular and tiny-variance profile endpoints agree with independent inversion, attain their cutoff, and preserve interval scaling without grid collapse")

  # Independent Gaussian truth identity via Owen's T integral, rather than
  # reusing the middle-class density integral in 07.
  normal_oracle <- function(d) {
    h <- d / sqrt(2)
    owens_t <- integrate(function(t) exp(-.5 * h^2 * (1 + t^2)) / (1 + t^2),
                         0, sqrt(3), rel.tol = 1e-12, abs.tol = 1e-12)$value / (2 * pi)
    pnorm(h) - 2 * owens_t
  }
  ds <- c(0, .65, .8, 1, 2, 5)
  close(env$rev_gaussian_truth(ds), vapply(ds, normal_oracle, numeric(1L)), tolerance = 1e-10)
  close(env$rev_gaussian_truth(0), 1 / 6)
  stopifnot(all(diff(env$rev_gaussian_truth(ds)) > 0))
  scenarios <- env$rev_scenarios()
  stopifnot(length(scenarios) == 12L, !anyDuplicated(vapply(scenarios, `[[`, "", "id")))
  for (s in scenarios) {
    set.seed(1908L); generated <- env$rev_generate_scores(s)
    stopifnot(identical(vapply(generated$scores, nrow, integer(1L)), as.integer(s$ns)),
              all(vapply(generated$scores, ncol, integer(1L)) == s$K),
              all(vapply(generated$scores, function(x) all(is.finite(x)), logical(1L))))
    truth <- env$rev_gaussian_truth(generated$d)
    P <- env$pair_contrast_matrix(s$K)$P
    zeros <- abs(as.vector(P %*% truth)) < 1e-12
    expected_zero_count <- if (s$alternative) choose(s$K - 1L, 2L) else choose(s$K, 2L)
    stopifnot(sum(zeros) == expected_zero_count)
  }
  pass("Gaussian VUS truth agrees with independent bivariate-normal/Owen-T identity; all 12 scenario dimensions and partial-null contrast counts are correct")

  for (id in c("S10", "HA", "LA", "K12N")) {
    s <- scenarios[[which(vapply(scenarios, `[[`, "", "id") == id)]]
    first <- env$rev_evaluate_replicate(s, 20261018L, 127L)
    second <- env$rev_evaluate_replicate(s, 20261018L, 127L)
    stochastic <- setdiff(names(first), c("t_vus", "t_scalar", "t_multiplier", "t_global", "t_total"))
    stopifnot(identical(first[stochastic], second[stochastic]), first$numerical_failure == 0L)
    indicators <- unlist(first[env$rev_metrics], use.names = FALSE)
    stopifnot(all(indicators[is.finite(indicators)] %in% c(0, 1)))
    if (!s$alternative) close(first$maxjel_fwer + first$maxjel_coverage, 1)
    # True-null FWER cannot exceed any-discovery across all nulls; and
    # intersection coverage is no larger than its Bonferroni-same-LR coverage
    # only when the random max cutoff is <= the Bonferroni cutoff.
    if (first$critical^2 <= qchisq(1 - .05 / choose(s$K, 2L), 1L)) {
      stopifnot(first$maxjel_coverage <= first$bonf_coverage)
    }
  }
  pass("fixed-seed end-to-end S10/HA/LA/K12N indicators reproduce exactly, remain binary, and satisfy complete-null coverage/FWER and cutoff-order identities")

  # Summary denominators must reflect metric-specific missing numerical results.
  s <- scenarios[[1L]]
  M <- data.frame(matrix(0, 3L, length(env$rev_metrics), dimnames = list(NULL, env$rev_metrics)))
  M$global_jel <- c(1, 0, NA_real_); M$numerical_failure <- c(0L, 0L, 1L)
  summary <- env$rev_summarize(M, s, B = 127L)
  item <- summary[summary$metric == "global_jel", ]
  close(item$estimate, .5); stopifnot(item$valid_replicates == 2L, item$failed_replicates == 1L)
  close(item$mcse, sqrt(.5 * .5 / 2)); close(item$worst_case_low, 1 / 3); close(item$worst_case_high, 2 / 3)
  pass("metric-specific failed replicates are excluded transparently, with correct MCSE denominators and worst-case failure bounds")

  log <- c(log, "", "All independent numerical checks passed.",
           "These checks validate finite-sample calculations and data flow, not universal statistical calibration.",
           "R reruns use R's RNG and need not match historical Python seeds or rates bit for bit.",
           "The corrected class-factor global-Wald covariance and explicit failure accounting may also change historical R aggregates.")
  directory <- file.path(root, "results"); dir.create(directory, recursive = TRUE, showWarnings = FALSE)
  writeLines(log, file.path(directory, "validation_report_R.txt"))
  cat("All independent numerical checks passed.\n")
  invisible(TRUE)
}

validate_revision()
