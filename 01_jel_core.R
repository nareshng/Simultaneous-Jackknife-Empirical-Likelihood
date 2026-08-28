# 01_jel_core.R
# Core estimators for simultaneous JEL inference with ordered three-class VUS.

# Tie convention:
#   kappa(x,y,z) = 1        if x < y < z
#                = 1/2      if x = y < z or x < y = z
#                = 1/6      if x = y = z
#                = 0        otherwise.
# This is the standard randomized ordering convention for an ordered triplet.

pair_contrast_matrix <- function(K) {
  stopifnot(K >= 2L)
  cmb <- utils::combn(K, 2L)
  P <- matrix(0, nrow = ncol(cmb), ncol = K)
  for (j in seq_len(ncol(cmb))) {
    P[j, cmb[1L, j]] <- 1
    P[j, cmb[2L, j]] <- -1
  }
  rownames(P) <- paste0(cmb[1L, ], "-", cmb[2L, ])
  list(pairs = t(cmb), P = P)
}

global_contrast_matrix <- function(K) {
  stopifnot(K >= 2L)
  C <- matrix(0, nrow = K - 1L, ncol = K)
  for (j in seq_len(K - 1L)) {
    C[j, j] <- 1
    C[j, K] <- -1
  }
  C
}

pooled_centering_vector <- function(ns) {
  ns <- as.integer(ns)
  if (length(ns) != 3L || any(ns <= 1L)) {
    stop("ns must contain three class sizes, all larger than one.")
  }
  n <- sum(ns)
  if (n <= 3L) stop("Total sample size must exceed three.")
  unlist(lapply(ns, function(m) {
    rep(n * (n - 2 * m - 1) / (m * (n - 3)), m)
  }), use.names = FALSE)
}

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

.el_score <- function(lambda, G) {
  den <- as.vector(1 + G %*% lambda)
  if (any(!is.finite(den)) || any(den <= 0)) {
    return(rep(NA_real_, ncol(G)))
  }
  colSums(G / den)
}

el_mean_statistic <- function(G, tol = 1e-10, max_iter = 200L,
                              ridge = 1e-12, return_details = FALSE) {
  G <- as.matrix(G)
  storage.mode(G) <- "double"
  if (nrow(G) <= ncol(G)) {
    stop("The number of pseudo-values must exceed the contrast dimension.")
  }
  if (any(!is.finite(G))) {
    out <- list(statistic = NA_real_, lambda = rep(NA_real_, ncol(G)),
                converged = FALSE, iterations = 0L,
                reason = "non-finite score matrix")
    return(if (return_details) out else out$statistic)
  }

  r <- ncol(G)
  lambda <- rep(0, r)
  score <- colSums(G)

  if (max(abs(score)) <= tol) {
    out <- list(statistic = 0, lambda = lambda, converged = TRUE,
                iterations = 0L, reason = "score already zero")
    return(if (return_details) out else out$statistic)
  }

  converged <- FALSE
  reason <- "maximum iterations reached"
  for (iter in seq_len(max_iter)) {
    den <- as.vector(1 + G %*% lambda)
    if (any(!is.finite(den)) || any(den <= 0)) {
      reason <- "left empirical-likelihood domain"
      break
    }
    W <- G / den
    score <- colSums(W)
    if (max(abs(score)) <= tol * (1 + sqrt(sum(G^2)))) {
      converged <- TRUE
      reason <- "converged"
      break
    }

    H <- crossprod(W)
    diag(H) <- diag(H) + ridge
    step <- tryCatch(
      solve(H, score),
      error = function(e) qr.solve(H, score, tol = 1e-12)
    )
    if (any(!is.finite(step))) {
      reason <- "singular Newton system"
      break
    }

    old_norm <- sqrt(sum(score^2))
    accepted <- FALSE
    scale <- 1
    for (ls in 0:40) {
      candidate <- lambda + scale * step
      den_new <- as.vector(1 + G %*% candidate)
      if (all(is.finite(den_new)) && min(den_new) > 1e-10) {
        score_new <- colSums(G / den_new)
        if (all(is.finite(score_new)) &&
            sqrt(sum(score_new^2)) <= old_norm * (1 - 1e-4 * scale)) {
          lambda <- candidate
          accepted <- TRUE
          break
        }
      }
      scale <- scale / 2
    }
    if (!accepted) {
      reason <- "line search failed"
      break
    }
  }

  den <- as.vector(1 + G %*% lambda)
  stat <- if (converged && all(den > 0)) 2 * sum(log(den)) else NA_real_
  if (is.finite(stat) && stat < 0 && stat > -1e-8) stat <- 0

  out <- list(
    statistic = stat,
    lambda = lambda,
    converged = converged,
    iterations = if (exists("iter")) iter else 0L,
    reason = reason,
    min_denominator = if (all(is.finite(den))) min(den) else NA_real_
  )
  if (return_details) out else out$statistic
}

jel_global_test <- function(theta_hat, pseudo_values, alpha = 0.05) {
  K <- length(theta_hat)
  C <- global_contrast_matrix(K)
  G <- pseudo_values %*% t(C)
  stat <- el_mean_statistic(G)
  list(
    statistic = stat,
    df = K - 1L,
    p_value = if (is.finite(stat)) stats::pchisq(stat, K - 1L, lower.tail = FALSE) else NA_real_,
    reject = is.finite(stat) && stat > stats::qchisq(1 - alpha, K - 1L),
    contrast = C
  )
}

wald_global_test <- function(theta_hat, pseudo_values, alpha = 0.05) {
  K <- length(theta_hat)
  C <- global_contrast_matrix(K)
  delta <- as.vector(C %*% theta_hat)
  G <- pseudo_values %*% t(C)
  S <- crossprod(scale(G, center = TRUE, scale = FALSE)) / nrow(G)
  eig <- eigen((S + t(S)) / 2, symmetric = TRUE)
  keep <- eig$values > max(eig$values) * 1e-10
  rank <- sum(keep)
  if (rank == 0L) {
    return(list(statistic = NA_real_, df = 0L, p_value = NA_real_, reject = NA))
  }
  Sinv <- eig$vectors[, keep, drop = FALSE] %*%
    diag(1 / eig$values[keep], rank, rank) %*%
    t(eig$vectors[, keep, drop = FALSE])
  stat <- nrow(G) * drop(t(delta) %*% Sinv %*% delta)
  list(
    statistic = stat,
    df = rank,
    p_value = stats::pchisq(stat, rank, lower.tail = FALSE),
    reject = stat > stats::qchisq(1 - alpha, rank)
  )
}

jel_pairwise_tests <- function(theta_hat, pseudo_values) {
  K <- length(theta_hat)
  pm <- pair_contrast_matrix(K)
  q <- nrow(pm$P)
  out <- vector("list", q)
  for (j in seq_len(q)) {
    cvec <- pm$P[j, ]
    delta <- drop(cvec %*% theta_hat)
    stat <- el_mean_statistic(matrix(pseudo_values %*% cvec, ncol = 1L))
    out[[j]] <- data.frame(
      pair_index = j,
      classifier_a = pm$pairs[j, 1L],
      classifier_b = pm$pairs[j, 2L],
      difference = delta,
      statistic = stat,
      signed_root = sign(delta) * sqrt(stat),
      p_value = stats::pchisq(stat, 1L, lower.tail = FALSE),
      stringsAsFactors = FALSE
    )
  }
  ans <- do.call(rbind, out)
  ans$holm_p <- stats::p.adjust(ans$p_value, method = "holm")
  ans$bonferroni_p <- stats::p.adjust(ans$p_value, method = "bonferroni")
  ans
}

max_jel_calibration <- function(theta_hat, pseudo_values, ns,
                                alpha = 0.05, B = 1999L, seed = 20260826L,
                                multiplier_matrix = NULL) {
  K <- length(theta_hat)
  pm <- pair_contrast_matrix(K)
  P <- pm$P
  delta <- as.vector(P %*% theta_hat)
  q <- length(delta)
  roots <- numeric(q)

  for (j in seq_len(q)) {
    stat <- el_mean_statistic(matrix(pseudo_values %*% P[j, ], ncol = 1L))
    roots[j] <- sign(delta[j]) * sqrt(stat)
  }

  a <- pooled_centering_vector(ns)
  G <- pseudo_values %*% t(P) - a %o% delta
  scale_j <- sqrt(colMeans(G^2))
  if (any(!is.finite(scale_j)) || any(scale_j <= 0)) {
    stop("At least one pairwise contrast has zero or undefined first-order variance.")
  }

  if (is.null(multiplier_matrix)) {
    set.seed(seed)
    Xi <- matrix(stats::rnorm(B * nrow(G)), nrow = B, ncol = nrow(G))
  } else {
    Xi <- as.matrix(multiplier_matrix)
    if (ncol(Xi) != nrow(G)) stop("multiplier_matrix has the wrong number of columns.")
    B <- nrow(Xi)
  }
  Zstar <- (Xi %*% G) / sqrt(nrow(G))
  Zstar <- sweep(Zstar, 2L, scale_j, "/")
  max_star <- apply(abs(Zstar), 1L, max)
  critical <- as.numeric(stats::quantile(max_star, 1 - alpha, type = 1L))
  adjusted_p <- vapply(abs(roots), function(t0) mean(max_star >= t0), numeric(1L))

  data.frame(
    pair_index = seq_len(q),
    classifier_a = pm$pairs[, 1L],
    classifier_b = pm$pairs[, 2L],
    difference = delta,
    signed_root = roots,
    critical_value = critical,
    adjusted_p = adjusted_p,
    reject = abs(roots) > critical,
    stringsAsFactors = FALSE
  )
}

jel_profile_interval <- function(vcontrast, a, estimate, critical,
                                 lower_bound = -1, upper_bound = 1,
                                 grid_size = 401L, tol = 1e-8) {
  vcontrast <- as.numeric(vcontrast)
  a <- as.numeric(a)
  if (length(vcontrast) != length(a)) stop("vcontrast and a must have equal length.")
  target <- critical^2
  stat_at <- function(delta0) {
    el_mean_statistic(matrix(vcontrast - a * delta0, ncol = 1L))
  }
  f <- function(delta0) stat_at(delta0) - target

  grid <- sort(unique(c(
    seq(lower_bound, upper_bound, length.out = grid_size),
    max(lower_bound, min(upper_bound, estimate))
  )))
  vals <- vapply(grid, f, numeric(1L))
  inside <- is.finite(vals) & vals <= 0
  if (!any(inside)) return(c(lower = NA_real_, upper = NA_real_))

  i0 <- which.min(abs(grid - estimate))
  if (!inside[i0]) {
    i0 <- which(inside)[which.min(abs(grid[inside] - estimate))]
  }

  left_idx <- i0
  while (left_idx > 1L && inside[left_idx - 1L]) left_idx <- left_idx - 1L
  right_idx <- i0
  while (right_idx < length(grid) && inside[right_idx + 1L]) right_idx <- right_idx + 1L

  lower <- grid[left_idx]
  if (left_idx > 1L && is.finite(vals[left_idx - 1L])) {
    lower <- tryCatch(
      stats::uniroot(f, c(grid[left_idx - 1L], grid[left_idx]), tol = tol)$root,
      error = function(e) grid[left_idx]
    )
  }
  upper <- grid[right_idx]
  if (right_idx < length(grid) && is.finite(vals[right_idx + 1L])) {
    upper <- tryCatch(
      stats::uniroot(f, c(grid[right_idx], grid[right_idx + 1L]), tol = tol)$root,
      error = function(e) grid[right_idx]
    )
  }
  c(lower = lower, upper = upper)
}

jel_simultaneous_intervals <- function(theta_hat, pseudo_values, ns,
                                       alpha = 0.05, B = 19999L,
                                       seed = 20260826L,
                                       method = c("maxJEL", "bonferroni")) {
  method <- match.arg(method)
  K <- length(theta_hat)
  pm <- pair_contrast_matrix(K)
  P <- pm$P
  q <- nrow(P)
  delta <- as.vector(P %*% theta_hat)
  a <- pooled_centering_vector(ns)

  if (method == "maxJEL") {
    cal <- max_jel_calibration(theta_hat, pseudo_values, ns, alpha, B, seed)
    critical <- unique(cal$critical_value)
  } else {
    critical <- sqrt(stats::qchisq(1 - alpha / q, df = 1L))
  }

  intervals <- matrix(NA_real_, q, 2L)
  for (j in seq_len(q)) {
    intervals[j, ] <- jel_profile_interval(
      vcontrast = as.vector(pseudo_values %*% P[j, ]),
      a = a,
      estimate = delta[j],
      critical = critical
    )
  }

  data.frame(
    pair_index = seq_len(q),
    classifier_a = pm$pairs[, 1L],
    classifier_b = pm$pairs[, 2L],
    difference = delta,
    lower = intervals[, 1L],
    upper = intervals[, 2L],
    method = method,
    critical_value = critical,
    stringsAsFactors = FALSE
  )
}

jel_self_check <- function() {
  X <- matrix(c(0, 1, 2, 3, 4, 5), ncol = 2)
  Y <- matrix(c(1, 2, 3, 4, 5, 6), ncol = 2)
  Z <- matrix(c(2, 3, 4, 5, 6, 7), ncol = 2)
  fit <- vus_and_pseudovalues(list(X, Y, Z), ties = TRUE)
  err <- max(abs(colMeans(fit$pseudo_values) - fit$theta_hat))
  if (err > 5e-8) stop("Self-check failed.")
  invisible(TRUE)
}
