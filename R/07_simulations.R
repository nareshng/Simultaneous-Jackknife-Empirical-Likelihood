# 07_revision_simulations.R
# Base-R port of the revision's 12 supplementary experiments.
# These are NEW R runs. R and Python RNGs/seeding do not produce identical data.
# No contributed R package, Python installation or online service is required.

rev_script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
rev_script_file <- if (length(rev_script_arg)) sub("^--file=", "", rev_script_arg[1L]) else ""
rev_script_file <- gsub("~+~", " ", rev_script_file, fixed = TRUE)
rev_source_files <- vapply(sys.frames(), function(frame) {
  if (exists("ofile", envir = frame, inherits = FALSE)) as.character(get("ofile", envir = frame))[1L] else ""
}, character(1L))
rev_source_files <- rev_source_files[nzchar(rev_source_files)]
rev_explicit_root <- Sys.getenv("JEL_PROJECT_ROOT", "")
rev_root <- if (nzchar(rev_explicit_root)) {
  normalizePath(rev_explicit_root, mustWork = TRUE)
} else if (length(rev_source_files)) {
  normalizePath(file.path(dirname(tail(rev_source_files, 1L)), ".."), mustWork = TRUE)
} else if (nzchar(rev_script_file)) {
  normalizePath(file.path(dirname(rev_script_file), ".."), mustWork = TRUE)
} else if (file.exists(file.path(getwd(), "R", "01_jel_core.R"))) {
  normalizePath(getwd(), mustWork = TRUE)
} else if (file.exists(file.path(getwd(), "01_jel_core.R"))) {
  normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
} else stop("Cannot locate package root; run from the R-only package root.")
if (!file.exists(file.path(rev_root, "R", "01_jel_core.R"))) stop("JEL project root does not contain R/01_jel_core.R: ", rev_root)
if (!exists("vus_and_pseudovalues", mode = "function", inherits = TRUE)) {
  source(file.path(rev_root, "R", "01_jel_core.R"), local = TRUE)
}

rev_scalar_jel <- function(G, iterations = 45L) {
  # A candidate outside the scalar convex hull has EL ratio zero: LR=Inf.
  # This is a valid infeasible constraint, not a numerical nonconvergence.
  G <- as.matrix(G)
  if (any(!is.finite(G))) stop("Non-finite scalar score matrix.")
  mn <- apply(G, 2L, min); mx <- apply(G, 2L, max)
  ans <- rep(Inf, ncol(G)); names(ans) <- colnames(G)
  zero <- pmax(abs(mn), abs(mx)) < 1e-14
  ans[zero] <- 0
  valid <- mn < 0 & mx > 0
  if (any(valid)) {
    g <- G[, valid, drop = FALSE]
    lower <- (-1 / mx[valid]) * (1 - 1e-12)
    upper <- (-1 / mn[valid]) * (1 - 1e-12)
    for (iteration in seq_len(iterations)) {
      lambda <- (lower + upper) / 2
      den <- 1 + sweep(g, 2L, lambda, "*")
      if (any(den <= 0)) stop("Scalar bisection left its feasible domain.")
      positive <- colSums(g / den) > 0
      lower[positive] <- lambda[positive]
      upper[!positive] <- lambda[!positive]
    }
    lambda <- (lower + upper) / 2
    ans[valid] <- pmax(0, 2 * colSums(log1p(sweep(g, 2L, lambda, "*"))))
  }
  ans
}

rev_global_jel <- function(G, max_iter = 100L) {
  # Concave dual; Newton decrement stops before floating-point line-search stalls.
  G <- as.matrix(G)
  failed <- function(reason, iteration = 0L) {
    list(statistic = NA_real_, converged = FALSE, iterations = iteration, reason = reason)
  }
  if (any(!is.finite(G))) return(failed("non-finite global score matrix"))
  scale <- sqrt(colMeans(G^2))
  if (any(!is.finite(scale)) || any(scale < 1e-12)) return(failed("degenerate global coordinate"))
  x <- sweep(G, 2L, scale, "/"); lambda <- rep(0, ncol(x))
  for (iteration in seq_len(max_iter)) {
    den <- as.vector(1 + x %*% lambda)
    if (any(!is.finite(den)) || any(den <= 0)) return(failed("left global dual domain", iteration))
    W <- x / den; grad <- colSums(W)
    value <- sum(log(den))
    # A small dual gradient alone also occurs at infinity for an infeasible
    # hull. A genuine EL solution must normalize its primal weights.
    primal_ok <- abs(mean(1 / den) - 1) < 1e-7
    if (max(abs(grad)) < 1e-7 && primal_ok) {
      return(list(statistic = max(0, 2 * value), converged = TRUE,
                  iterations = iteration, reason = "gradient tolerance"))
    }
    H <- crossprod(W)
    step <- tryCatch(solve(H, grad), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) return(failed("singular global Newton system", iteration))
    decrement <- sum(grad * step)
    if (decrement >= 0 && decrement < 1e-10 && primal_ok) {
      return(list(statistic = max(0, 2 * value), converged = TRUE,
                  iterations = iteration, reason = "Newton decrement tolerance"))
    }
    if (!is.finite(decrement) || decrement < 0) return(failed("invalid global Newton direction", iteration))
    accepted <- FALSE; rate <- 1
    for (ls in 0:40) {
      trial <- lambda + rate * step; d2 <- as.vector(1 + x %*% trial)
      if (all(is.finite(d2)) && all(d2 > 0) &&
          sum(log(d2)) >= value + 1e-4 * rate * decrement) {
        lambda <- trial; accepted <- TRUE; break
      }
      rate <- rate / 2
    }
    if (!accepted) return(failed("global line search failed", iteration))
  }
  failed("global maximum iterations reached", max_iter)
}

rev_gaussian_truth <- local({
  cache <- new.env(parent = emptyenv())
  function(d) {
    vapply(d, function(dd) {
      key <- sprintf("%.17g", dd)
      if (!exists(key, envir = cache, inherits = FALSE)) {
        val <- integrate(function(y) dnorm(y - dd) * pnorm(y) * pnorm(2 * dd - y),
                         lower = -Inf, upper = Inf, abs.tol = 1e-12, rel.tol = 1e-10,
                         subdivisions = 1000L, stop.on.error = TRUE)$value
        assign(key, val, envir = cache)
      }
      get(key, envir = cache, inherits = FALSE)
    }, numeric(1L), USE.NAMES = FALSE)
  }
})

rev_scenarios <- function() {
  make <- function(id, K, ns, dependence, alternative) {
    list(id = id, K = K, ns = as.integer(ns), dependence = dependence, alternative = alternative)
  }
  out <- lapply(c(10L, 20L, 40L, 80L), function(m) {
    make(paste0("S", m), 4L, c(m, 40, 160), "equicorrelated", FALSE)
  })
  for (dep in c("class_specific", "nonlinear")) {
    for (alt in c(FALSE, TRUE)) {
      out[[length(out) + 1L]] <- make(paste0(if (dep == "class_specific") "H" else "L", if (alt) "A" else "N"),
                                    4L, c(40, 40, 40), dep, alt)
    }
  }
  for (K in c(8L, 12L)) {
    for (alt in c(FALSE, TRUE)) {
      out[[length(out) + 1L]] <- make(paste0("K", K, if (alt) "A" else "N"), K,
                                    c(60, 60, 60), "class_specific", alt)
    }
  }
  out
}

rev_generate_scores <- function(s) {
  K <- s$K; d <- rep(0.8, K)
  if (isTRUE(s$alternative)) d[K] <- 1
  out <- vector("list", 3L)
  for (cc in 0:2) {
    n <- s$ns[cc + 1L]; E <- matrix(rnorm(n * K), n, K)
    if (s$dependence == "equicorrelated") {
      noise <- (E + rnorm(n)) / sqrt(2)
    } else if (s$dependence == "class_specific") {
      z <- ((0:(K - 1L)) - cc - 1L) %% K
      load <- cbind(.2 + .55 * z / (K - 1L), .25 * sin(2 * pi * (z + cc) / K))
      noise <- matrix(rnorm(2 * n), n, 2L) %*% t(load) +
        sweep(E, 2L, sqrt(1 - rowSums(load^2)), "*")
    } else if (s$dependence == "nonlinear") {
      shifts <- ((0:(K - 1L)) * (cc + 1L) * .173) %% 1
      uniforms <- outer(pnorm(rnorm(n)), shifts, "+") %% 1
      shared <- qnorm(pmax(pmin(uniforms, 1 - 1e-15), 1e-15))
      noise <- sqrt(.65) * shared + sqrt(.35) * E
    } else stop("Unknown revision dependence: ", s$dependence)
    out[[cc + 1L]] <- sweep(noise, 2L, cc * d, "+")
  }
  list(scores = out, d = d)
}

rev_max_calibration <- function(U, V, a, B, alpha = .05, standard_normals = NULL) {
  n <- nrow(V); K <- ncol(V); pm <- pair_contrast_matrix(K)
  first <- pm$pairs[, 1L]; second <- pm$pairs[, 2L]
  G <- V - a %o% U
  covariance <- crossprod(G) / n
  covariance <- (covariance + t(covariance)) / 2
  variance <- diag(covariance)[first] + diag(covariance)[second] -
    2 * covariance[cbind(first, second)]
  if (any(!is.finite(variance)) || any(variance <= 1e-24)) stop("Degenerate pairwise variance.")
  eig <- eigen(covariance, symmetric = TRUE)
  if (min(eig$values) < -1e-10 * max(1, max(abs(eig$values)))) stop("Non-positive-semidefinite covariance.")
  if (is.null(standard_normals)) standard_normals <- matrix(rnorm(B * K), B, K)
  if (!identical(dim(standard_normals), c(as.integer(B), as.integer(K)))) stop("standard_normals must be B by K.")
  draws <- sweep(standard_normals, 2L, sqrt(pmax(eig$values, 0)), "*") %*% t(eig$vectors)
  contrasted <- draws[, first, drop = FALSE] - draws[, second, drop = FALSE]
  maximum <- apply(abs(sweep(contrasted, 2L, sqrt(variance), "/")), 1L, max)
  list(critical = unname(quantile(maximum, 1 - alpha, type = 1L)), covariance = covariance,
       variance = variance, maxima = maximum, delta = U[first] - U[second],
       centered = G, pairs = pm$pairs, gaussian_draws = draws)
}

rev_metrics <- c("global_jel", "global_wald", "unadjusted_fwer", "holm_fwer", "maxjel_fwer", "waldmax_fwer",
                 "holm_discovery", "maxjel_discovery", "waldmax_discovery", "maxjel_coverage", "bonf_coverage")

rev_evaluate_replicate <- function(s, seed, B, alpha = .05) {
  set.seed(seed)
  row <- as.list(setNames(rep(NA_real_, length(rev_metrics)), rev_metrics))
  row <- c(list(seed = seed), row, list(numerical_failure = 0L, global_numerical_failure = 0L,
       wald_numerical_failure = 0L, scalar_infeasible = 0L, truth_scalar_infeasible = 0L,
       global_statistic = NA_real_, global_wald_statistic = NA_real_, critical = NA_real_,
       t_vus = NA_real_, t_scalar = NA_real_, t_multiplier = NA_real_, t_global = NA_real_,
       t_total = NA_real_, error_message = ""))
  tryCatch({
    generated <- rev_generate_scores(s); K <- s$K; d <- generated$d; n <- sum(s$ns)
    pairs <- pair_contrast_matrix(K)$pairs
    first <- pairs[, 1L]; second <- pairs[, 2L]; q <- nrow(pairs)
    t0 <- proc.time()[["elapsed"]]
    fit <- vus_and_pseudovalues(generated$scores, ties = TRUE)
    U <- fit$theta_hat; V <- fit$pseudo_values; a <- pooled_centering_vector(fit$ns)
    t1 <- proc.time()[["elapsed"]]
    Pv <- V[, first, drop = FALSE] - V[, second, drop = FALSE]
    LR <- rev_scalar_jel(Pv); t2 <- proc.time()[["elapsed"]]
    truth_theta <- rev_gaussian_truth(d); truth <- truth_theta[first] - truth_theta[second]
    LRtrue <- if (isTRUE(s$alternative)) rev_scalar_jel(Pv - a %o% truth) else LR
    maxcal <- rev_max_calibration(U, V, a, B = B, alpha = alpha)
    t3 <- proc.time()[["elapsed"]]
    p <- pchisq(LR, df = 1L, lower.tail = FALSE)
    reject_holm <- p.adjust(p, method = "holm") < alpha
    reject_max <- LR > maxcal$critical^2
    reject_wald <- abs(sqrt(n) * maxcal$delta / sqrt(maxcal$variance)) > maxcal$critical
    null <- abs(d[first] - d[second]) < 1e-12; alternative <- !null
    global <- rev_global_jel(V[, seq_len(K - 1L), drop = FALSE] - V[, K])
    row$global_jel <- if (global$converged) as.integer(global$statistic > qchisq(1 - alpha, K - 1L)) else NA_integer_
    row$global_numerical_failure <- as.integer(!global$converged)
    row$global_statistic <- global$statistic
    if (!global$converged) row$error_message <- global$reason
    C <- global_contrast_matrix(K); gd <- as.vector(C %*% U)
    gcov <- C %*% maxcal$covariance %*% t(C)
    gwald <- tryCatch(n * sum(gd * solve(gcov, gd)), error = function(e) NA_real_)
    row$global_wald_statistic <- gwald
    row$global_wald <- if (is.finite(gwald)) as.integer(gwald > qchisq(1 - alpha, K - 1L)) else NA_integer_
    row$wald_numerical_failure <- as.integer(!is.finite(gwald))
    if (!is.finite(gwald)) row$error_message <- paste(row$error_message, "global Wald failed")
    row$unadjusted_fwer <- as.integer(any((p < alpha)[null]))
    row$holm_fwer <- as.integer(any(reject_holm[null]))
    row$maxjel_fwer <- as.integer(any(reject_max[null]))
    row$waldmax_fwer <- as.integer(any(reject_wald[null]))
    if (any(alternative)) {
      row$holm_discovery <- as.integer(any(reject_holm[alternative]))
      row$maxjel_discovery <- as.integer(any(reject_max[alternative]))
      row$waldmax_discovery <- as.integer(any(reject_wald[alternative]))
    }
    row$maxjel_coverage <- as.integer(all(LRtrue <= maxcal$critical^2))
    row$bonf_coverage <- as.integer(all(LRtrue <= qchisq(1 - alpha / q, 1L)))
    row$scalar_infeasible <- sum(is.infinite(LR))
    row$truth_scalar_infeasible <- sum(is.infinite(LRtrue))
    row$numerical_failure <- as.integer(!global$converged || !is.finite(gwald))
    row$critical <- maxcal$critical
    t4 <- proc.time()[["elapsed"]]
    row$t_vus <- t1 - t0; row$t_scalar <- t2 - t1; row$t_multiplier <- t3 - t2
    row$t_global <- t4 - t3; row$t_total <- t4 - t0
    row
  }, error = function(e) {
    row$numerical_failure <- 1L
    row$error_message <- conditionMessage(e)
    row
  })
}

rev_wilson <- function(p, n) {
  if (n <= 0L || !is.finite(p)) return(c(low = NA_real_, high = NA_real_))
  z <- qnorm(.975); den <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / den
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / den
  c(low = max(0, center - half), high = min(1, center + half))
}

rev_summarize <- function(M, s, B) {
  applicable <- rev_metrics
  if (!isTRUE(s$alternative)) applicable <- setdiff(applicable, c("holm_discovery", "maxjel_discovery", "waldmax_discovery"))
  do.call(rbind, lapply(applicable, function(metric) {
    v <- M[[metric]]; valid <- is.finite(v); nv <- sum(valid); total <- nrow(M)
    estimate <- if (nv) mean(v[valid]) else NA_real_; interval <- rev_wilson(estimate, nv)
    successes <- sum(v[valid]); missing <- total - nv
    data.frame(scenario = s$id, K = s$K, ns = paste(s$ns, collapse = "/"), dependence = s$dependence,
      alternative = s$alternative, metric = metric, estimate = estimate,
      mcse = if (nv) sqrt(estimate * (1 - estimate) / nv) else NA_real_,
      low = unname(interval["low"]), high = unname(interval["high"]), reps = nv, B = B,
      valid_replicates = nv, failed_replicates = missing, total_replicates = total,
      numerical_failures = sum(M$numerical_failure != 0L),
      worst_case_low = successes / total, worst_case_high = (successes + missing) / total,
      stringsAsFactors = FALSE)
  }))
}

rev_write_csv <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  write.csv(x, tmp, row.names = FALSE, na = "NA")
  if (!file.rename(tmp, path)) {
    if (!file.copy(tmp, path, overwrite = TRUE)) stop("Cannot write ", path)
    unlink(tmp)
  }
}

rev_helper_self_check <- function() {
  stopifnot(abs(rev_gaussian_truth(0) - 1 / 6) < 1e-10)
  stopifnot(abs(rev_gaussian_truth(.8) - .4565862906658332) < 1e-10)
  stopifnot(abs(rev_scalar_jel(cbind(c(-1, 1), c(-2, 2)))) < 1e-10)
  stopifnot(is.infinite(rev_scalar_jel(c(1, 2, 3))))
  g <- rbind(diag(3), -diag(3)); fit <- rev_global_jel(g)
  stopifnot(fit$converged, abs(fit$statistic) < 1e-10, length(rev_scenarios()) == 12L)
  message("Revision helper self-check passed.")
  invisible(TRUE)
}

rev_run <- function() {
  envint <- function(name, default, minimum = 1L) {
    value <- suppressWarnings(as.integer(Sys.getenv(name, as.character(default))))
    if (length(value) != 1L || is.na(value) || value < minimum) stop("Invalid ", name)
    value
  }
  reps <- envint("REV_REPS", 5000L); B <- envint("REV_B_MULT", 1999L, 2L)
  seed <- envint("REV_SEED", 20261016L); cores <- envint("MC_CORES", 1L)
  runtime_reps <- envint("REV_RUNTIME_REPS", 500L)
  every <- envint("REV_CHECKPOINT_EVERY", 100L)
  alpha <- as.numeric(Sys.getenv("ALPHA", ".05"))
  if (!is.finite(alpha) || alpha <= 0 || alpha >= 1) stop("ALPHA must lie in (0,1).")
  suffix <- Sys.getenv("OUTPUT_SUFFIX", "_R")
  resume <- tolower(Sys.getenv("REV_RESUME", "true")) %in% c("1", "true", "yes")
  out <- file.path(rev_root, "results"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
  scenarios <- rev_scenarios(); all_ids <- vapply(scenarios, `[[`, "", "id")
  selected <- trimws(strsplit(Sys.getenv("REV_SCENARIOS", paste(all_ids, collapse = ",")), ",", fixed = TRUE)[[1L]])
  if (any(!selected %in% all_ids)) stop("Unknown REV_SCENARIOS identifier.")
  # Set and record RNG kind so a resume never silently changes the R stream.
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  config <- list(reps = reps, B = B, seed = seed, alpha = alpha, runtime_reps = runtime_reps,
    rng = RNGkind(), R_version = R.version.string, selected = selected,
    source_md5 = unname(tools::md5sum(file.path(rev_root, "R", c("01_jel_core.R", "07_revision_simulations.R")))))
  config_path <- file.path(out, paste0("revision_config", suffix, ".rds"))
  if (resume && file.exists(config_path) && !identical(readRDS(config_path), config)) {
    stop("Revision resume configuration/source changed. Choose a new OUTPUT_SUFFIX, or set REV_RESUME=false to rerun this suffix.")
  }
  saveRDS(config, config_path)
  writeLines(c(capture.output(sessionInfo()), "", capture.output(str(config))),
             file.path(out, paste0("revision_sessionInfo", suffix, ".txt")))
  registry <- do.call(rbind, lapply(seq_along(scenarios), function(j) {
    s <- scenarios[[j]]
    data.frame(scenario = s$id, index = j - 1L, K = s$K, n1 = s$ns[1L], n2 = s$ns[2L], n3 = s$ns[3L],
      dependence = s$dependence, alternative = s$alternative,
      seed_rule = paste0(seed, " + ", j - 1L, "*1000000 + (replicate-1)"), stringsAsFactors = FALSE)
  }))
  rev_write_csv(registry, file.path(out, paste0("revision_scenario_registry", suffix, ".csv")))
  rev_helper_self_check()
  summaries <- list(); runtimes <- list(); paired <- list()
  for (index in seq_along(scenarios)) {
    s <- scenarios[[index]]; if (!s$id %in% selected) next
    path <- file.path(out, paste0("revision_replicates_", s$id, suffix, ".csv"))
    M <- NULL
    if (resume && file.exists(path)) {
      M <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
      if (!all(c("replicate", "seed", rev_metrics) %in% names(M)) || nrow(M) > reps ||
          !identical(as.integer(M$replicate), seq_len(nrow(M)))) stop("Invalid checkpoint: ", path)
      expected_seeds <- seed + (index - 1L) * 1000000L + seq_len(nrow(M)) - 1L
      if (any(M$seed != expected_seeds)) stop("Checkpoint seed mismatch: ", path)
    }
    completed <- if (is.null(M)) 0L else nrow(M)
    message("Revision ", s$id, ": ", completed, "/", reps, " existing replications; B=", B)
    if (completed < reps) {
      for (start in seq.int(completed + 1L, reps, by = every)) {
        ids <- seq.int(start, min(reps, start + every - 1L))
        one <- function(rr) {
          ss <- seed + (index - 1L) * 1000000L + rr - 1L
          value <- rev_evaluate_replicate(s, ss, B, alpha)
          as.data.frame(c(list(replicate = rr), value), stringsAsFactors = FALSE)
        }
        values <- if (cores > 1L && .Platform$OS.type != "windows") {
          parallel::mclapply(ids, one, mc.cores = cores, mc.preschedule = TRUE, mc.set.seed = FALSE)
        } else lapply(ids, one)
        if (any(vapply(values, inherits, logical(1L), what = "try-error"))) stop("Parallel worker failed; previous checkpoint preserved.")
        added <- do.call(rbind, values)
        M <- if (is.null(M)) added else rbind(M, added)
        rev_write_csv(M, path)
        message("  ", s$id, " checkpoint ", nrow(M), "/", reps)
      }
    }
    summary <- rev_summarize(M, s, B); summaries[[length(summaries) + 1L]] <- summary
    times <- M$t_total[is.finite(M$t_total)]
    runtimes[[length(runtimes) + 1L]] <- data.frame(scenario = s$id, K = s$K, q = s$K * (s$K - 1L) / 2,
      n = sum(s$ns), reps = nrow(M), valid_timings = length(times),
      median_seconds = if (length(times)) median(times) else NA_real_,
      p95_seconds = if (length(times)) unname(quantile(times, .95)) else NA_real_,
      numerical_failures = sum(M$numerical_failure != 0L),
      global_numerical_failures = sum(M$global_numerical_failure != 0L),
      scalar_infeasible = sum(M$scalar_infeasible), truth_scalar_infeasible = sum(M$truth_scalar_infeasible))
    if (s$alternative) {
      valid <- is.finite(M$maxjel_discovery) & is.finite(M$holm_discovery)
      dif <- M$maxjel_discovery[valid] - M$holm_discovery[valid]
      gain <- if (length(dif)) mean(dif) else NA_real_
      half <- if (length(dif) > 1L) qnorm(.975) * sd(dif) / sqrt(length(dif)) else NA_real_
      paired[[length(paired) + 1L]] <- data.frame(scenario = s$id, gain = gain, low = gain - half, high = gain + half,
        reps = length(dif), failed_replicates = nrow(M) - length(dif))
    }
    rev_write_csv(do.call(rbind, summaries), file.path(out, paste0("revision_summary", suffix, ".csv")))
    rev_write_csv(do.call(rbind, runtimes), file.path(out, paste0("revision_runtime", suffix, ".csv")))
    if (length(paired)) rev_write_csv(do.call(rbind, paired), file.path(out, paste0("revision_paired_power_gains", suffix, ".csv")))
    key <- summary[summary$metric %in% c("maxjel_fwer", "maxjel_coverage", "maxjel_discovery"),
                   c("metric", "estimate", "low", "high", "failed_replicates")]
    print(key, row.names = FALSE)
    if (any(summary$failed_replicates)) warning(s$id, ": numerical failures; valid-only estimates and worst-case failure bounds recorded, not counted as nonrejections.")
  }
  # Controlled timing holds n and B fixed across K; no classifier fitting or
  # profile endpoint inversion is included. New R timings replace no Python timing.
  timing <- list()
  for (K in c(4L, 8L, 12L)) for (alternative in c(FALSE, TRUE)) {
    s <- list(id = paste0("B", K, if (alternative) "A" else "N"), K = K,
              ns = c(60L, 60L, 60L), dependence = "class_specific", alternative = alternative)
    values <- lapply(seq_len(runtime_reps), function(r) rev_evaluate_replicate(s,
      202610160L + K * 100000L + as.integer(alternative) * 10000L + r - 1L, B, alpha))
    M <- do.call(rbind, lapply(values, as.data.frame, stringsAsFactors = FALSE))
    valid <- is.finite(M$t_total); times <- M$t_total[valid]
    med <- function(name) if (any(is.finite(M[[name]]))) median(M[[name]][is.finite(M[[name]])]) * 1000 else NA_real_
    timing[[length(timing) + 1L]] <- data.frame(K = K, q = K * (K - 1L) / 2, n = 180L, B = B,
      alternative = alternative, repetitions = runtime_reps, valid_timings = length(times),
      median_ms = med("t_total"), p95_ms = if (length(times)) unname(quantile(times, .95)) * 1000 else NA_real_,
      vus_median_ms = med("t_vus"), scalar_null_median_ms = med("t_scalar"),
      multiplier_and_truth_median_ms = med("t_multiplier"), global_median_ms = med("t_global"),
      numerical_failures = sum(M$numerical_failure != 0L))
    rev_write_csv(do.call(rbind, timing), file.path(out, paste0("revision_controlled_runtime", suffix, ".csv")))
  }
  message("Revision R experiments completed. Do not substitute earlier Python rates or timings for these R outputs.")
  invisible(TRUE)
}

if (sys.nframe() == 0L && !tolower(Sys.getenv("REVISION_HELPERS_ONLY", "false")) %in% c("true", "1", "yes")) {
  rev_run()
}
