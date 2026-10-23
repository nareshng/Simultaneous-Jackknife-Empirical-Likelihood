# 06_make_figures.R
# Figures from CSV result files using base R.

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
results_dir <- file.path(package_root, "results")
figures_dir <- file.path(package_root, "figures")
dir.create(figures_dir, recursive=TRUE, showWarnings=FALSE)
suffix <- Sys.getenv("OUTPUT_SUFFIX", "_R")
alpha <- as.numeric(Sys.getenv("ALPHA", "0.05"))

use_file <- function(stem) {
  # Exact suffix only: never mix an old unsuffixed run with current results.
  path <- file.path(results_dir, paste0(stem, suffix, ".csv"))
  if (!file.exists(path)) stop("Missing result file for the selected OUTPUT_SUFFIX: ", path,
    ". Rerun the matching original/full/quick workflow first.")
  path
}

# Derive the wide plotting tables from the long-format file written by 04.
# MCSE bars describe Monte Carlo uncertainty of each rate, not uncertainty of
# a paired difference between two procedures.
simulation_table <- function(kind) {
  map <- if (kind == "null") c(
    Global_JEL="global_jel_reject", Global_Wald="global_wald_reject",
    Unadjusted_any="unadjusted_any", Holm_JEL="holm_any",
    Max_JEL="maxjel_any", Wald_maxT="waldmax_any"
  ) else c(
    Holm_any_true="holm_any_true_discovery",
    Max_JEL_any_true="maxjel_any_true_discovery",
    Wald_maxT_any_true="waldmax_any_true_discovery"
  )
  aggregate_path <- use_file("simulation_aggregate_results")
  aggregate <- utils::read.csv(aggregate_path, check.names=FALSE)
  required <- c("ScenarioID","Scenario","Metric","Estimate","MCSE")
  if (!all(required %in% names(aggregate))) stop("Incomplete simulation aggregate table.")
  prefix <- if (kind == "null") "N" else "A"
  selected <- aggregate[startsWith(aggregate$ScenarioID,prefix),,drop=FALSE]
  ids <- unique(selected$ScenarioID)
  if (!length(ids)) stop("No ",kind," scenarios in aggregate table.")
  out <- data.frame(ScenarioID=ids,stringsAsFactors=FALSE)
  for (nm in names(map)) {
    z <- selected[selected$Metric == map[[nm]],,drop=FALSE]
    if (anyDuplicated(z$ScenarioID) || !all(ids %in% z$ScenarioID)) {
      stop("Missing or duplicated aggregate metric: ",map[[nm]])
    }
    ix <- match(ids,z$ScenarioID)
    out[[nm]] <- z$Estimate[ix]
    out[[paste0(nm,"_MCSE")]] <- z$MCSE[ix]
    if ("Valid_replicates" %in% names(z)) out[[paste0(nm,"_valid_replicates")]] <- z$Valid_replicates[ix]
    if ("Numerical_failures" %in% names(z)) out[[paste0(nm,"_numerical_failures")]] <- z$Numerical_failures[ix]
  }
  utils::write.csv(out,file.path(results_dir,paste0("simulation_",kind,
                    "_results",suffix,".csv")),row.names=FALSE)
  out
}

rate_plot <- function(dat, methods, labels, title, ylab, output, null=FALSE) {
  estimate <- as.matrix(dat[,methods,drop=FALSE])
  se <- as.matrix(dat[,paste0(methods,"_MCSE"),drop=FALSE])
  lower <- pmax(estimate-1.96*se,0)
  upper <- pmin(estimate+1.96*se,1)
  nx <- nrow(estimate); nm <- ncol(estimate)
  xx <- outer(seq_len(nx),seq(-.18,.18,length.out=nm),"+")
  png(file.path(figures_dir,output),width=2400,height=1450,res=180)
  op <- par(mar=c(5,5,4,1))
  ylim <- if (null) c(0,max(.26,estimate,upper,na.rm=TRUE)) else c(0,1)
  matplot(xx,estimate,type="b",lty=1,pch=seq_len(nm),col=seq_len(nm),
          xaxt="n",xlab=if (null) "Null scenario" else "Alternative scenario",
          ylab=ylab,ylim=ylim,xlim=c(.5,nx+.5),main=title)
  for (j in seq_len(nm)) {
    good <- is.finite(lower[,j]) & is.finite(upper[,j])
    if (any(good)) arrows(xx[good,j],lower[good,j],xx[good,j],upper[good,j],
                         angle=90,code=3,length=.025,col=j)
  }
  axis(1,at=seq_len(nx),labels=dat$ScenarioID)
  if (null) abline(h=alpha,lty=2,lwd=2)
  legend("topleft",legend=labels,lty=1,pch=seq_len(nm),col=seq_len(nm),
         cex=.82,ncol=if (null) 2 else 1,bty="n")
  mtext("Bars: pointwise estimate +/- 1.96 Monte Carlo SE",side=3,line=.2,cex=.75)
  par(op); dev.off()
}

# Fig. 1: null type-I error and FWER, with pointwise MC uncertainty.
null <- simulation_table("null")
methods <- c("Global_JEL","Global_Wald","Unadjusted_any","Holm_JEL","Max_JEL","Wald_maxT")
labels <- c("Global JEL","Global Wald","Unadjusted pairwise","Holm-JEL","max-JEL","Wald maxT")
rate_plot(null,methods,labels,
          "Type-I error and familywise error under the global null",
          "Rejection probability",paste0("fig1_null_fwer",suffix,".png"),null=TRUE)

# Fig. 2: power with the same pointwise MC uncertainty convention.
alt <- simulation_table("alternative")
methods2 <- c("Holm_any_true","Max_JEL_any_true","Wald_maxT_any_true")
labels2 <- c("Holm-JEL","max-JEL","Wald maxT")
rate_plot(alt,methods2,labels2,"Multiplicity-adjusted power under alternatives",
          "Probability of at least one true discovery",
          paste0("fig2_alternative_power",suffix,".png"))

# Fig. 3: real-data VUS.
perf <- utils::read.csv(use_file("real_data_model_performance"),check.names=FALSE)
png(file.path(figures_dir,paste0("fig3_real_vus",suffix,".png")),width=1800,height=1100,res=180)
op <- par(mar=c(8,5,3,1))
plot(seq_len(nrow(perf)),perf$Heldout_VUS,pch=19,cex=1.5,xaxt="n",
     ylim=range(c(1/6,perf$Heldout_VUS))+c(-.03,.03),
     xlab="",ylab="Held-out VUS",
     main="Classifier VUS on the independent test cohort")
axis(1,at=seq_len(nrow(perf)),labels=perf$Classifier,las=2)
abline(h=1/6,lty=2)
par(op)
dev.off()

# Fig. 4: simultaneous pairwise intervals.
pair <- utils::read.csv(use_file("real_data_pairwise_inference"),check.names=FALSE)
lo_name <- if ("Max_JEL_CI_lower" %in% names(pair)) "Max_JEL_CI_lower" else "max_JEL_CI_lower"
hi_name <- if ("Max_JEL_CI_upper" %in% names(pair)) "Max_JEL_CI_upper" else "max_JEL_CI_upper"
diff_name <- if ("VUS_difference_A_minus_B" %in% names(pair)) "VUS_difference_A_minus_B" else "difference"
labels4 <- paste(pair$Classifier_A,"-",pair$Classifier_B)
ord <- rev(seq_len(nrow(pair)))
png(file.path(figures_dir,paste0("fig4_pairwise_intervals",suffix,".png")),width=2100,height=1350,res=180)
op <- par(mar=c(5,16,3,1))
xlim <- range(c(pair[[lo_name]],pair[[hi_name]],0), na.rm=TRUE)
plot(pair[[diff_name]][ord],seq_along(ord),xlim=xlim+c(-.02,.02),
     ylim=c(.5,length(ord)+.5),yaxt="n",ylab="",
     xlab="VUS difference (A - B)",pch=19,
     main=sprintf("%.1f%% simultaneous max-JEL intervals", 100*(1-alpha)))
segments(pair[[lo_name]][ord],seq_along(ord),
         pair[[hi_name]][ord],seq_along(ord),lwd=2)
axis(2,at=seq_along(ord),labels=labels4[ord],las=2,cex.axis=.75)
abline(v=0,lty=2)
par(op)
dev.off()
message("Figures written to ",figures_dir)

