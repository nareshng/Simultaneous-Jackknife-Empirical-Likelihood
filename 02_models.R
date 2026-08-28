# 02_models.R
# Deterministic base-R classifiers used in the urinary-biomarker application.
# No test observation is used to estimate transformations, tune hyperparameters,
# or fit a classifier.

stable_softmax <- function(eta) {
  eta <- as.matrix(eta)
  eta <- eta - apply(eta, 1L, max)
  e <- exp(eta)
  e / rowSums(e)
}

multiclass_log_loss <- function(prob, y, eps = 1e-15) {
  prob <- as.matrix(prob)
  y <- as.integer(y)
  p <- pmin(1 - eps, pmax(eps, prob[cbind(seq_along(y), y + 1L)]))
  -mean(log(p))
}

classification_accuracy <- function(prob, y) {
  pred <- max.col(prob, ties.method = "first") - 1L
  mean(pred == as.integer(y))
}

fit_preprocessor <- function(df) {
  required <- c("age", "creatinine", "LYVE1", "REG1B", "TFF1")
  if (!all(required %in% names(df))) {
    stop("Missing required predictors: ", paste(setdiff(required, names(df)), collapse = ", "))
  }
  raw <- cbind(
    age = as.numeric(df$age),
    log_creatinine = log(pmax(as.numeric(df$creatinine), .Machine$double.xmin)),
    log_LYVE1 = log(pmax(as.numeric(df$LYVE1), .Machine$double.xmin)),
    log_REG1B = log(pmax(as.numeric(df$REG1B), .Machine$double.xmin)),
    log_TFF1 = log(pmax(as.numeric(df$TFF1), .Machine$double.xmin))
  )
  center <- colMeans(raw)
  scale <- apply(raw, 2L, stats::sd)
  scale[!is.finite(scale) | scale <= 0] <- 1
  list(center = center, scale = scale, variables = required)
}

apply_preprocessor <- function(df, prep) {
  raw <- cbind(
    age = as.numeric(df$age),
    log_creatinine = log(pmax(as.numeric(df$creatinine), .Machine$double.xmin)),
    log_LYVE1 = log(pmax(as.numeric(df$LYVE1), .Machine$double.xmin)),
    log_REG1B = log(pmax(as.numeric(df$REG1B), .Machine$double.xmin)),
    log_TFF1 = log(pmax(as.numeric(df$TFF1), .Machine$double.xmin))
  )
  X <- sweep(raw, 2L, prep$center, "-")
  X <- sweep(X, 2L, prep$scale, "/")
  storage.mode(X) <- "double"
  X
}

fit_ridge_multinomial <- function(X, y, lambda = 1, maxit = 1000L) {
  X <- as.matrix(X); y <- as.integer(y)
  classes <- sort(unique(y))
  if (!identical(classes, 0:(length(classes) - 1L))) {
    stop("Classes must be coded as consecutive integers starting at zero.")
  }
  K <- length(classes)
  D <- cbind(`(Intercept)` = 1, X)
  p <- ncol(D)
  n <- nrow(D)
  Y <- matrix(0, nrow = n, ncol = K)
  Y[cbind(seq_len(n), y + 1L)] <- 1

  unpack <- function(par) matrix(par, nrow = p, ncol = K - 1L)
  objective <- function(par) {
    B <- unpack(par)
    eta <- cbind(D %*% B, 0)
    prob <- stable_softmax(eta)
    nll <- -sum(log(pmax(prob[cbind(seq_len(n), y + 1L)], 1e-15)))
    penalty <- 0.5 * lambda * sum(B[-1L, , drop = FALSE]^2)
    nll + penalty
  }
  gradient <- function(par) {
    B <- unpack(par)
    eta <- cbind(D %*% B, 0)
    prob <- stable_softmax(eta)
    G <- crossprod(D, prob[, seq_len(K - 1L), drop = FALSE] -
                     Y[, seq_len(K - 1L), drop = FALSE])
    G[-1L, ] <- G[-1L, ] + lambda * B[-1L, ]
    as.vector(G)
  }

  fit <- stats::optim(
    par = rep(0, p * (K - 1L)),
    fn = objective,
    gr = gradient,
    method = "BFGS",
    control = list(maxit = maxit, reltol = 1e-10)
  )
  if (fit$convergence != 0L) {
    warning("Ridge multinomial optimizer did not report convergence: code ", fit$convergence)
  }
  list(type = "ridge_multinomial", coefficients = unpack(fit$par),
       lambda = lambda, K = K, convergence = fit$convergence,
       objective = fit$value)
}

predict_ridge_multinomial <- function(model, X) {
  D <- cbind(`(Intercept)` = 1, as.matrix(X))
  stable_softmax(cbind(D %*% model$coefficients, 0))
}

.safe_cov <- function(X) {
  X <- as.matrix(X)
  if (nrow(X) <= 1L) return(diag(0, ncol(X)))
  S <- stats::cov(X)
  if (is.null(dim(S))) S <- matrix(S, 1L, 1L)
  (S + t(S)) / 2
}

fit_regularized_qda <- function(X, y, lambda = 1) {
  X <- as.matrix(X); y <- as.integer(y)
  K <- max(y) + 1L
  p <- ncol(X)
  means <- matrix(NA_real_, nrow = K, ncol = p)
  covariances <- vector("list", K)
  priors <- numeric(K)
  for (k in 0:(K - 1L)) {
    Xk <- X[y == k, , drop = FALSE]
    means[k + 1L, ] <- colMeans(Xk)
    S <- .safe_cov(Xk) + lambda * diag(p)
    eig <- eigen((S + t(S)) / 2, symmetric = TRUE, only.values = TRUE)$values
    if (min(eig) <= 1e-10) S <- S + (1e-10 - min(eig) + 1e-8) * diag(p)
    covariances[[k + 1L]] <- S
    priors[k + 1L] <- nrow(Xk) / nrow(X)
  }
  list(type = "regularized_qda", means = means, covariances = covariances,
       priors = priors, lambda = lambda, K = K)
}

predict_regularized_qda <- function(model, X) {
  X <- as.matrix(X)
  n <- nrow(X); K <- model$K; p <- ncol(X)
  logscore <- matrix(NA_real_, n, K)
  constant <- p * log(2 * pi)
  for (k in seq_len(K)) {
    S <- model$covariances[[k]]
    R <- chol(S)
    centered <- sweep(X, 2L, model$means[k, ], "-")
    z <- forwardsolve(t(R), t(centered))
    mahal <- colSums(z^2)
    logdet <- 2 * sum(log(diag(R)))
    logscore[, k] <- log(pmax(model$priors[k], 1e-15)) -
      0.5 * (constant + logdet + mahal)
  }
  stable_softmax(logscore)
}

rbf_kernel <- function(X1, X2, gamma) {
  X1 <- as.matrix(X1); X2 <- as.matrix(X2)
  d2 <- outer(rowSums(X1^2), rowSums(X2^2), "+") - 2 * tcrossprod(X1, X2)
  d2[d2 < 0 & d2 > -1e-10] <- 0
  exp(-gamma * pmax(d2, 0))
}

fit_rbf_krr <- function(X, y, gamma = 0.08, lambda = 0.1) {
  X <- as.matrix(X); y <- as.integer(y)
  n <- nrow(X); Kclass <- max(y) + 1L
  Y <- matrix(0, n, Kclass)
  Y[cbind(seq_len(n), y + 1L)] <- 1
  Kmat <- rbf_kernel(X, X, gamma)
  A <- solve(Kmat + lambda * diag(n), Y)
  list(type = "rbf_krr", X = X, alpha = A, gamma = gamma,
       lambda = lambda, K = Kclass)
}

predict_rbf_krr <- function(model, X) {
  scores <- rbf_kernel(as.matrix(X), model$X, model$gamma) %*% model$alpha
  stable_softmax(scores)
}

fit_weighted_knn <- function(X, y, k = 25L) {
  X <- as.matrix(X); y <- as.integer(y)
  if (k < 1L || k > nrow(X)) stop("k must be between 1 and the training sample size.")
  list(type = "weighted_knn", X = X, y = y, k = as.integer(k),
       K = max(y) + 1L)
}

predict_weighted_knn <- function(model, X) {
  X <- as.matrix(X)
  d2 <- outer(rowSums(X^2), rowSums(model$X^2), "+") - 2 * tcrossprod(X, model$X)
  d <- sqrt(pmax(d2, 0))
  out <- matrix(0, nrow(X), model$K)
  for (i in seq_len(nrow(X))) {
    ord <- order(d[i, ], seq_along(d[i, ]))[seq_len(model$k)]
    w <- 1 / (d[i, ord] + 1e-8)
    for (kk in 0:(model$K - 1L)) {
      out[i, kk + 1L] <- sum(w[model$y[ord] == kk])
    }
    out[i, ] <- out[i, ] / sum(out[i, ])
  }
  out
}

fit_classifier <- function(df, y, method, hyperparameters) {
  prep <- fit_preprocessor(df)
  X <- apply_preprocessor(df, prep)
  model <- switch(
    method,
    ridge_multinomial = fit_ridge_multinomial(
      X, y, lambda = as.numeric(hyperparameters$lambda)
    ),
    regularized_qda = fit_regularized_qda(
      X, y, lambda = as.numeric(hyperparameters$lambda)
    ),
    rbf_krr = fit_rbf_krr(
      X, y, gamma = as.numeric(hyperparameters$gamma),
      lambda = as.numeric(hyperparameters$lambda)
    ),
    weighted_knn = fit_weighted_knn(
      X, y, k = as.integer(hyperparameters$k)
    ),
    stop("Unknown method: ", method)
  )
  list(method = method, preprocessing = prep, model = model,
       hyperparameters = hyperparameters)
}

predict_classifier <- function(fitted, new_df) {
  X <- apply_preprocessor(new_df, fitted$preprocessing)
  switch(
    fitted$method,
    ridge_multinomial = predict_ridge_multinomial(fitted$model, X),
    regularized_qda = predict_regularized_qda(fitted$model, X),
    rbf_krr = predict_rbf_krr(fitted$model, X),
    weighted_knn = predict_weighted_knn(fitted$model, X),
    stop("Unknown fitted method: ", fitted$method)
  )
}

make_stratified_folds <- function(y, v = 5L, seed = 20260826L) {
  y <- as.integer(y)
  set.seed(seed)
  fold <- integer(length(y))
  for (cls in sort(unique(y))) {
    ii <- which(y == cls)
    ii <- sample(ii, length(ii), replace = FALSE)
    fold[ii] <- rep(seq_len(v), length.out = length(ii))
  }
  fold
}

default_parameter_grid <- function(method) {
  switch(
    method,
    ridge_multinomial = data.frame(lambda = c(0.01, 0.1, 1, 10)),
    regularized_qda = data.frame(lambda = c(0.01, 0.1, 0.5, 1, 2)),
    rbf_krr = expand.grid(
      gamma = c(0.02, 0.05, 0.08, 0.12, 0.20),
      lambda = c(0.01, 0.1, 1),
      KEEP.OUT.ATTRS = FALSE
    ),
    weighted_knn = data.frame(k = c(5L, 10L, 15L, 25L, 35L)),
    stop("Unknown method: ", method)
  )
}

cross_validate_classifier <- function(df, y, method,
                                      grid = default_parameter_grid(method),
                                      v = 5L, seed = 20260826L,
                                      verbose = TRUE) {
  fold <- make_stratified_folds(y, v = v, seed = seed)
  scores <- rep(NA_real_, nrow(grid))
  failures <- integer(nrow(grid))

  for (g in seq_len(nrow(grid))) {
    hyper <- as.list(grid[g, , drop = FALSE])
    fold_loss <- rep(NA_real_, v)
    for (ff in seq_len(v)) {
      train <- fold != ff
      valid <- !train
      ans <- tryCatch({
        fitted <- fit_classifier(df[train, , drop = FALSE], y[train], method, hyper)
        prob <- predict_classifier(fitted, df[valid, , drop = FALSE])
        multiclass_log_loss(prob, y[valid])
      }, error = function(e) {
        if (verbose) {
          message("CV failure for ", method, ", grid row ", g,
                  ", fold ", ff, ": ", conditionMessage(e))
        }
        NA_real_
      })
      fold_loss[ff] <- ans
    }
    failures[g] <- sum(!is.finite(fold_loss))
    scores[g] <- if (failures[g] == 0L) mean(fold_loss) else Inf
    if (verbose) {
      message(sprintf("%s grid %d/%d: mean log loss %.6f",
                      method, g, nrow(grid), scores[g]))
    }
  }

  result <- cbind(grid, mean_log_loss = scores, failures = failures)
  best <- which.min(result$mean_log_loss)
  if (!is.finite(result$mean_log_loss[best])) {
    stop("All hyperparameter settings failed for method ", method)
  }
  list(
    method = method,
    best_hyperparameters = as.list(grid[best, , drop = FALSE]),
    best_cv_log_loss = result$mean_log_loss[best],
    table = result,
    folds = fold
  )
}

ordinal_score_from_probabilities <- function(prob) {
  prob <- as.matrix(prob)
  drop(prob %*% (0:(ncol(prob) - 1L)))
}
