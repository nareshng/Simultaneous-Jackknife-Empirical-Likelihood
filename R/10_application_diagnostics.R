# Application diagnostics for the simultaneous VUS/JEL manuscript.
# Base R only: no classifier fitting, tuning, or simulation.
#
# Run from your analysis project directory:
#   Rscript R/10_application_diagnostics.R
# Or specify input and an optional output directory:
#   Rscript R/10_application_diagnostics.R path/to/scores.csv path/to/output
#
# Default input: results/real_data_heldout_scores{OUTPUT_SUFFIX}.csv
# Standalone OUTPUT_SUFFIX default: _full_R (the completed manuscript run).
# The workflow driver passes explicit input/output paths for its own suffix.
# Default output: application_diagnostics{input suffix}/ beside the input.
# The VUS routine below is copied unchanged from the supplied
# revised_R_analysis/R/01_jel_core.R.

vus_and_pseudovalues <- function(scores_by_class, ties = TRUE) {
  if (length(scores_by_class) != 3L) {
    stop("scores_by_class must be a list of three matrices.")
  }
  X <- as.matrix(scores_by_class[[1L]])
  Y <- as.matrix(scores_by_class[[2L]])
  Z <- as.matrix(scores_by_class[[3L]])
  if (ncol(X) != ncol(Y) || ncol(X) != ncol(Z)) {
    stop("All score matrices must have the same number of classifiers.")
  }
  ns <- c(nrow(X), nrow(Y), nrow(Z))
  if (any(ns <= 1L)) stop("Each class must contain at least two observations.")
  n1 <- ns[1L]; n2 <- ns[2L]; n3 <- ns[3L]
  n <- sum(ns); K <- ncol(X)

  U <- numeric(K)
  V <- matrix(NA_real_, nrow = n, ncol = K)
  colnames(V) <- colnames(X)

  for (k in seq_len(K)) {
    x <- X[, k]; y <- Y[, k]; z <- Z[, k]

    lt_xy <- outer(x, y, "<")
    gt_zy <- outer(z, y, ">")
    if (ties) {
      eq_xy <- outer(x, y, "==")
      eq_zy <- outer(z, y, "==")
      w_xy <- lt_xy + 0.5 * eq_xy
      w_zy <- gt_zy + 0.5 * eq_zy
      A <- colSums(w_xy)
      B <- colSums(w_zy)
      EX <- colSums(eq_xy)
      EZ <- colSums(eq_zy)
      middle_contribution <- A * B - (EX * EZ) / 12
      S <- sum(middle_contribution)

      remove_x <- rowSums(
        sweep(w_xy, 2L, B, "*") -
          sweep(eq_xy, 2L, EZ / 12, "*")
      )
      remove_y <- middle_contribution
      remove_z <- rowSums(
        sweep(w_zy, 2L, A, "*") -
          sweep(eq_zy, 2L, EX / 12, "*")
      )
    } else {
      A <- colSums(lt_xy)
      B <- colSums(gt_zy)
      S <- sum(A * B)
      remove_x <- rowSums(sweep(lt_xy, 2L, B, "*"))
      remove_y <- A * B
      remove_z <- rowSums(sweep(gt_zy, 2L, A, "*"))
    }

    u <- S / (n1 * n2 * n3)
    U[k] <- u

    raw_minus_x <- (S - remove_x) / ((n1 - 1L) * n2 * n3)
    raw_minus_y <- (S - remove_y) / (n1 * (n2 - 1L) * n3)
    raw_minus_z <- (S - remove_z) / (n1 * n2 * (n3 - 1L))

    # Delete-one value for the transformed pooled three-sample U-statistic.
    pooled_factor <- n / (n - 3L)
    pooled_minus_x <- pooled_factor * (n1 - 1L) / n1 * raw_minus_x
    pooled_minus_y <- pooled_factor * (n2 - 1L) / n2 * raw_minus_y
    pooled_minus_z <- pooled_factor * (n3 - 1L) / n3 * raw_minus_z

    V[, k] <- c(
      n * u - (n - 1L) * pooled_minus_x,
      n * u - (n - 1L) * pooled_minus_y,
      n * u - (n - 1L) * pooled_minus_z
    )
  }

  identity_error <- max(abs(colMeans(V) - U))
  if (!is.finite(identity_error) || identity_error > 5e-8) {
    stop(sprintf("Pseudo-value identity failed: max error = %.3e", identity_error))
  }

  list(theta_hat = U, pseudo_values = V, ns = ns)
}

run_application_diagnostics <- function(input_file, output_dir = NULL) {
  if (!file.exists(input_file)) stop("Saved score file not found: ", input_file)
  # Match companion tables to this score file, never to an older run.
  score_name <- basename(input_file)
  score_pattern <- "^real_data_heldout_scores(.*)[.]csv$"
  input_suffix <- if (grepl(score_pattern, score_name)) {
    sub(score_pattern, "\\1", score_name)
  } else NULL
  companion_file <- function(stem) {
    if (is.null(input_suffix)) return(NA_character_)
    file.path(dirname(input_file), paste0(stem, input_suffix, ".csv"))
  }
  if (is.null(input_suffix)) {
    message("Nonstandard score filename: companion-output consistency checks skipped.")
  }
  dat <- utils::read.csv(input_file, check.names = FALSE,
                         stringsAsFactors = FALSE)
  methods <- c("Ridge multinomial", "Regularized QDA",
               "RBF kernel ridge", "Distance-weighted kNN")
  short_names <- c("Ridge", "QDA", "RBF", "kNN")
  required <- c("ordered_class", methods)
  if (!all(required %in% names(dat))) {
    stop("Missing columns: ", paste(setdiff(required, names(dat)), collapse = ", "))
  }
  if (anyNA(dat$ordered_class)) stop("Missing ordered_class values.")
  levels_y <- sort(unique(dat$ordered_class))
  if (!(identical(as.numeric(levels_y), c(0, 1, 2)) ||
        identical(as.numeric(levels_y), c(1, 2, 3)))) {
    stop("ordered_class must contain all of 0,1,2 or all of 1,2,3.")
  }
  if ("sample_id" %in% names(dat) && anyDuplicated(dat$sample_id)) {
    stop("Duplicate sample_id values: inspect these before proceeding.")
  }
  if (!all(vapply(dat[methods], is.numeric, logical(1L)))) {
    stop("Every classifier-score column must be numeric.")
  }
  scores <- as.matrix(dat[methods])
  if (any(!is.finite(scores))) stop("Missing or nonfinite classifier scores.")
  groups <- lapply(levels_y, function(y) {
    scores[dat$ordered_class == y, , drop = FALSE]
  })
  ns <- vapply(groups, nrow, integer(1L))
  if (any(ns < 2L)) stop("Each diagnostic class needs at least two subjects.")
  n <- sum(ns)
  K <- length(methods)
  fit <- vus_and_pseudovalues(groups, ties = TRUE)
  theta <- fit$theta_hat
  V <- fit$pseudo_values

  # V, a, and class_index all use the same concatenated class order.
  a_s <- n * (n - 2 * ns - 1) / (ns * (n - 3))
  a <- rep(a_s, times = ns)
  H <- V - a %o% theta
  if (max(abs(colMeans(H))) > 5e-8 || abs(mean(a) - 1) > 1e-12) {
    stop("Class-centering identity failed.")
  }
  pairs <- t(utils::combn(seq_len(K), 2L))
  q <- nrow(pairs)
  P <- matrix(0, q, K)
  P[cbind(seq_len(q), pairs[, 1L])] <- 1
  P[cbind(seq_len(q), pairs[, 2L])] <- -1
  G <- H %*% t(P)
  delta <- as.vector(P %*% theta)
  s2 <- colMeans(G^2)
  positive <- is.finite(s2) & s2 > 0
  influence_max <- rep(NA_real_, q)
  for (j in which(positive)) {
    influence_max[j] <- max(abs(G[, j])) / sqrt(sum(G[, j]^2))
  }
  class_index <- rep(seq_along(ns), times = ns)
  max_class <- vapply(seq_len(q), function(j) {
    if (!positive[j]) return(NA_character_)
    c("Control", "Benign", "PDAC")[class_index[which.max(abs(G[, j]))]]
  }, character(1L))
  result <- data.frame(
    Contrast = paste(short_names[pairs[, 1L]], short_names[pairs[, 2L]], sep = " - "),
    Classifier_A = methods[pairs[, 1L]],
    Classifier_B = methods[pairs[, 2L]],
    VUS_difference = delta,
    Influence_variance_s2 = s2,
    Estimated_variance_of_difference = s2 / n,
    SE_of_difference = sqrt(s2 / n),
    Standardized_influence_max = influence_max,
    Largest_squared_influence_share = influence_max^2,
    Class_of_largest_absolute_influence = max_class,
    Positive_contrast_variance = positive,
    stringsAsFactors = FALSE
  )

  # Fixed omnibus basis: first K-1 classifiers minus the final classifier.
  # The full six-pair covariance is structurally singular.
  C <- cbind(diag(K - 1L), rep(-1, K - 1L))
  HC <- H %*% t(C)
  omega <- crossprod(HC) / n
  omnibus_names <- paste(short_names[seq_len(K - 1L)], "kNN", sep = " - ")
  dimnames(omega) <- list(omnibus_names, omnibus_names)
  eigenvalues <- eigen((omega + t(omega)) / 2, symmetric = TRUE,
                       only.values = TRUE)$values
  lambda_min <- min(eigenvalues)
  lambda_max <- max(eigenvalues)
  condition_number <- if (lambda_min > 0) lambda_max / lambda_min else Inf
  rank_tolerance <- if (lambda_max > 0) 1e-10 * lambda_max else 0
  numerical_rank <- sum(eigenvalues > rank_tolerance)

  # Verify consistency with saved application outputs when they are present.
  checks <- character()
  pair_file <- companion_file("real_data_pairwise_inference")
  if (!is.na(pair_file) && file.exists(pair_file)) {
    old <- utils::read.csv(pair_file, check.names = FALSE, stringsAsFactors = FALSE)
    required_old <- c("Classifier_A", "Classifier_B", "VUS_difference_A_minus_B")
    if (!all(required_old %in% names(old))) stop("Unexpected saved pairwise columns.")
    key <- function(A, B) paste(A, B, sep = " || ")
    idx <- match(key(result$Classifier_A, result$Classifier_B),
                 key(old$Classifier_A, old$Classifier_B))
    if (anyNA(idx)) stop("Saved pairwise output has a different classifier family.")
    err <- max(abs(delta - old$VUS_difference_A_minus_B[idx]))
    if (!is.finite(err) || err > 1e-9) {
      stop("Recomputed VUS differences do not match saved results.")
    }
    checks <- c(checks, sprintf("Six saved VUS differences matched; maximum error %.3e.", err))
  }
  perf_file <- companion_file("real_data_model_performance")
  if (!is.na(perf_file) && file.exists(perf_file)) {
    old <- utils::read.csv(perf_file, check.names = FALSE, stringsAsFactors = FALSE)
    if (!all(c("Classifier", "Heldout_VUS") %in% names(old))) {
      stop("Unexpected saved performance columns.")
    }
    idx <- match(methods, old$Classifier)
    if (anyNA(idx)) stop("Saved model-performance output has different classifiers.")
    err <- max(abs(theta - old$Heldout_VUS[idx]))
    if (!is.finite(err) || err > 1e-9) {
      stop("Recomputed VUS estimates do not match saved results.")
    }
    checks <- c(checks, sprintf("Four saved VUS estimates matched; maximum error %.3e.", err))
  }

  if (is.null(output_dir)) {
    output_name <- paste0("application_diagnostics",
                          if (is.null(input_suffix)) "" else input_suffix)
    output_dir <- file.path(dirname(input_file), output_name)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  utils::write.csv(result, file.path(output_dir, "application_diagnostics.csv"),
                   row.names = FALSE)
  utils::write.csv(omega, file.path(output_dir, "omnibus_covariance.csv"),
                   row.names = TRUE)
  C_export <- data.frame(Contrast = omnibus_names, C, check.names = FALSE)
  names(C_export)[-1L] <- short_names
  utils::write.csv(C_export, file.path(output_dir, "omnibus_contrast_matrix.csv"),
                   row.names = FALSE)
  summary_lines <- c(
    "APPLICATION DIAGNOSTICS FROM SAVED HELD-OUT SCORES",
    paste("Input:", normalizePath(input_file, mustWork = TRUE)),
    paste("Input MD5:", unname(tools::md5sum(input_file))),
    sprintf("n=%d; control=%d; benign=%d; PDAC=%d", n, ns[1L], ns[2L], ns[3L]),
    paste("Class order:", paste(levels_y, collapse = " < ")),
    "Tie convention: strict ordering=1; one adjacent tie=1/2; all equal=1/6.",
    "Centering: H_l = V_l - a_l theta_hat; G_l = P H_l.",
    "s2_j = mean(G_lj^2); estimated variance of Delta_hat_j = s2_j/n.",
    "M_j = max(abs(G_lj))/sqrt(sum(G_lj^2)).",
    "M_j^2 is the largest single subject's share of the squared centered influences.",
    "Omnibus covariance: mean((C H_l)(C H_l)^T); C compares Ridge, QDA, RBF against kNN.",
    paste("Omnibus eigenvalues:", paste(format(eigenvalues, digits = 12), collapse = ", ")),
    sprintf("Smallest eigenvalue: %.12g", lambda_min),
    sprintf("Largest eigenvalue: %.12g", lambda_max),
    sprintf("Omnibus condition number (lambda_max/lambda_min): %.12g", condition_number),
    sprintf("Numerical rank at relative tolerance 1e-10: %d of %d", numerical_rank, K - 1L),
    sprintf("Contrasts with nonpositive/nonfinite variance: %d", sum(!positive)),
    checks,
    "Condition numbers depend on the stated contrast basis and are descriptive.",
    "No directions are dropped and no diagnostic threshold certifies calibration.",
    "These diagnostics do not establish coverage or stability to deleting observations.",
    "Classifier fitting, tuning, simulations, tests and interval inversion were not rerun.",
    "Unique sample IDs alone do not establish person-level independence."
  )
  writeLines(c(summary_lines, "", capture.output(print(result, row.names = FALSE, digits = 6)),
               "", capture.output(sessionInfo())),
             file.path(output_dir, "application_diagnostics_summary.txt"))

  # Ready-to-insert LaTeX table; full precision remains in the CSV.
  rows <- vapply(seq_len(q), function(j) {
    paste0(sprintf("%s--%s & %.6f & %.6f & %.6f & %.4f ",
                   short_names[pairs[j, 1L]], short_names[pairs[j, 2L]],
                   delta[j], s2[j], sqrt(s2[j] / n), influence_max[j]), "\\\\")
  }, character(1L))
  tex <- c(
    "\\begin{table}[!t]",
    "\\caption{Diagnostics for the six held-out VUS contrasts.}",
    "\\label{tab:application-diagnostics}",
    "\\centering",
    "\\begin{tabular}{lrrrr}",
    "\\hline",
    paste0("Contrast & $\\widehat\\Delta$ & $\\widehat s^2$ & SE & $M$ ", "\\\\"),
    "\\hline", rows, "\\hline", "\\end{tabular}",
    "\\par\\smallskip",
    "\\begin{minipage}{\\columnwidth}\\footnotesize",
    sprintf("$n=%d$. Here $\\widehat s_j^2=n^{-1}\\sum_lG_{lj}^2$,", n),
    "$\\mathrm{SE}_j=\\widehat s_j/\\sqrt n$, and",
    "$M_j=\\max_l|G_{lj}|/(\\sum_lG_{lj}^2)^{1/2}$.",
    sprintf("For the three omnibus contrasts against kNN, the covariance condition number is %.4f", condition_number),
    sprintf("and the smallest eigenvalue is %.6g.", lambda_min),
    "These are descriptive diagnostics and do not certify finite-sample calibration.",
    "\\end{minipage}", "\\end{table}"
  )
  writeLines(tex, file.path(output_dir, "application_diagnostics_table.tex"))
  cat(paste(summary_lines, collapse = "\n"), "\n\n")
  print(result[, c("Contrast", "VUS_difference", "Influence_variance_s2",
                  "SE_of_difference", "Standardized_influence_max")],
        row.names = FALSE, digits = 6)
  cat("\nOutputs written to:", normalizePath(output_dir), "\n")
  invisible(list(diagnostics = result, omnibus_covariance = omega,
                 eigenvalues = eigenvalues, condition_number = condition_number))
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) > 2L) {
  stop("Usage: Rscript R/10_application_diagnostics.R [scores.csv] [output_directory]")
}
default_suffix <- Sys.getenv("OUTPUT_SUFFIX", unset = "_full_R")
default_input <- file.path("results", paste0(
  "real_data_heldout_scores", default_suffix, ".csv"))
if (!file.exists(default_input) && file.exists(basename(default_input))) {
  default_input <- basename(default_input)
}
input_file <- if (length(args) >= 1L) args[1L] else default_input
output_dir <- if (length(args) >= 2L) args[2L] else NULL
run_application_diagnostics(input_file, output_dir)

