# 05_interval_simulation.R
# Simultaneous confidence-interval coverage and length.
#
# Manuscript run:
#   B_MULT=2000 


script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
if (is.na(script_file) || !nzchar(script_file)) {
  package_root <- normalizePath(getwd())
} else {
  package_root <- normalizePath(file.path(dirname(script_file), ".."))
}
source(file.path(package_root, "R", "01_jel_core.R"))

results_dir <- file.path(package_root, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

Bmult <- as.integer(Sys.getenv("B_MULT", "2000"))
cores <- as.integer(Sys.getenv("MC_CORES", "1"))
master_seed <- as.integer(Sys.getenv("INT_SEED", "2026"))
suffix <- Sys.getenv("OUTPUT_SUFFIX", "_R")
override_reps <- Sys.getenv("INT_REPS", "")

scenarios <- list(
  list(id="N3", ns=c(40,40,40), rho=.8, dist="gaussian",
       d=rep(.80,4), reps=1000L),
  list(id="N4", ns=c(102,119,37), rho=.5, dist="gaussian",
       d=rep(.80,4), reps=500L),
  list(id="N6", ns=c(40,40,40), rho=.5, dist="t3",
       d=rep(.80,4), reps=500L),
  list(id="N9", ns=c(40,40,40), rho=.5, dist="rounded",
       d=rep(.80,4), reps=500L),
  list(id="A3", ns=c(40,40,40), rho=.5, dist="gaussian",
       d=c(.80,.80,.80,1.00), reps=500L)
)
if (nzchar(override_reps)) {
  for (j in seq_along(scenarios)) scenarios[[j]]$reps <- as.integer(override_reps)
}

generate_correlated <- function(n, K, rho) {
  S <- matrix(rho, K, K); diag(S) <- 1
  matrix(stats::rnorm(n*K), n, K) %*% chol(S)
}
generate_scores <- function(s) {
  K <- length(s$d)
  lapply(0:2, function(cc) {
    n <- s$ns[cc+1L]
    E <- generate_correlated(n, K, s$rho)
    if (s$dist == "t3") E <- E / sqrt(stats::rchisq(n, 3) / 3)
    X <- sweep(E, 2L, cc * s$d, "+")
    if (s$dist == "rounded") X <- round(X / .5) * .5
    X
  })
}

gaussian_vus <- function(d) {
  stats::integrate(
    function(y) stats::dnorm(y, mean=d, sd=1) *
      stats::pnorm(y, mean=0, sd=1) *
      stats::pnorm(y, mean=2*d, sd=1, lower.tail=FALSE),
    lower = -Inf, upper = Inf,
    rel.tol = 1e-11, subdivisions = 1000L
  )$value
}

scenario_truth <- function(s) {
  K <- length(s$d)
  pm <- pair_contrast_matrix(K)
  if (length(unique(s$d)) == 1L) return(rep(0, nrow(pm$P)))
  if (s$dist != "gaussian") {
    stop("Non-null interval truth is implemented analytically only for Gaussian margins.")
  }
  theta <- vapply(s$d, gaussian_vus, numeric(1L))
  as.vector(pm$P %*% theta)
}

run_one <- function(s, ss) {
  truth <- scenario_truth(s)
  seeds <- master_seed + ss * 1000000L + seq_len(s$reps)
  worker <- function(rr) {
    set.seed(seeds[rr])
    scores <- generate_scores(s)
    tryCatch({
      fit <- vus_and_pseudovalues(scores, ties=TRUE)
      max_ci <- jel_simultaneous_intervals(
        fit$theta_hat, fit$pseudo_values, fit$ns,
        alpha=.05, B=Bmult, seed=seeds[rr]+500000L, method="maxJEL"
      )
      bonf_ci <- jel_simultaneous_intervals(
        fit$theta_hat, fit$pseudo_values, fit$ns,
        alpha=.05, B=Bmult, seed=seeds[rr]+500000L, method="bonferroni"
      )
      c(
        replicate=rr, seed=seeds[rr], failure=0,
        max_coverage=as.numeric(all(truth >= max_ci$lower & truth <= max_ci$upper)),
        bonf_coverage=as.numeric(all(truth >= bonf_ci$lower & truth <= bonf_ci$upper)),
        max_mean_length=mean(max_ci$upper-max_ci$lower),
        bonf_mean_length=mean(bonf_ci$upper-bonf_ci$lower)
      )
    }, error=function(e) {
      c(replicate=rr, seed=seeds[rr], failure=1,
        max_coverage=NA, bonf_coverage=NA,
        max_mean_length=NA, bonf_mean_length=NA)
    })
  }
  if (cores > 1L && .Platform$OS.type != "windows") {
    z <- parallel::mclapply(seq_len(s$reps), worker, mc.cores=cores)
  } else {
    z <- lapply(seq_len(s$reps), worker)
  }
  M <- as.data.frame(do.call(rbind,z))
  M$ScenarioID <- s$id
  M <- M[,c("ScenarioID",setdiff(names(M),"ScenarioID"))]
  utils::write.csv(
    M, file.path(results_dir,paste0("interval_replicates_",s$id,suffix,".csv")),
    row.names=FALSE
  )
  nvalid <- sum(M$failure==0)
  pmax <- mean(M$max_coverage,na.rm=TRUE)
  pbonf <- mean(M$bonf_coverage,na.rm=TRUE)
  data.frame(
    ScenarioID=s$id, Replications=s$reps, Valid_replicates=nvalid,
    Numerical_failures=sum(M$failure==1),
    Max_JEL_simultaneous_coverage=pmax,
    Max_JEL_coverage_MCSE=sqrt(pmax*(1-pmax)/nvalid),
    Bonferroni_JEL_simultaneous_coverage=pbonf,
    Bonferroni_JEL_coverage_MCSE=sqrt(pbonf*(1-pbonf)/nvalid),
    Max_JEL_mean_length=mean(M$max_mean_length,na.rm=TRUE),
    Bonferroni_JEL_mean_length=mean(M$bonf_mean_length,na.rm=TRUE),
    stringsAsFactors=FALSE
  )
}

out <- do.call(rbind,lapply(seq_along(scenarios),function(j) run_one(scenarios[[j]],j)))
utils::write.csv(
  out, file.path(results_dir,paste0("simulation_interval_results",suffix,".csv")),
  row.names=FALSE
)
message("Interval simulation completed.")
