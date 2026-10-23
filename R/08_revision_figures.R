# 08_revision_figures.R
# Publication figures and a TeX calibration table from R-generated revision results.
# Base R only. Run after 07_revision_simulations.R, from any working directory:
#   OUTPUT_SUFFIX=_R Rscript R/08_revision_figures.R
# No Python result files or manuscript point estimates are read by this script.

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[1L]) else ""
script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
source_files <- vapply(sys.frames(), function(frame) {
  if (exists("ofile", envir = frame, inherits = FALSE)) {
    as.character(get("ofile", envir = frame))[1L]
  } else ""
}, character(1L))
source_files <- source_files[nzchar(source_files)]
explicit_root <- Sys.getenv("JEL_PROJECT_ROOT", "")
if (nzchar(explicit_root)) {
  package_root <- normalizePath(explicit_root, mustWork = TRUE)
} else if (length(source_files)) {
  package_root <- normalizePath(file.path(dirname(tail(source_files, 1L)), ".."), mustWork = TRUE)
} else if (nzchar(script_file)) {
  package_root <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
} else if (file.exists(file.path(getwd(), "R", "01_jel_core.R"))) {
  package_root <- normalizePath(getwd(), mustWork = TRUE)
} else if (file.exists(file.path(getwd(), "01_jel_core.R"))) {
  package_root <- normalizePath(file.path(getwd(), ".."), mustWork = TRUE)
} else {
  stop("Run this file with Rscript or source it from the R package root.")
}
if (!file.exists(file.path(package_root, "R", "01_jel_core.R"))) {
  stop("JEL project root does not contain R/01_jel_core.R: ", package_root)
}
source(file.path(package_root, "R", "01_jel_core.R"), local = TRUE)
suffix <- Sys.getenv("OUTPUT_SUFFIX", "_R")
if (grepl("[/\\\\]", suffix)) stop("OUTPUT_SUFFIX cannot contain path separators.")
results_dir <- file.path(package_root, "results")
figures_dir <- file.path(package_root, "figures")
input <- file.path(results_dir, paste0("revision_summary", suffix, ".csv"))
if (!file.exists(input)) {
  stop("R revision summary does not exist: ", input,
       ". Run 07_revision_simulations.R with the same OUTPUT_SUFFIX first.")
}
dat <- utils::read.csv(input, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c("scenario", "K", "ns", "dependence", "alternative", "metric",
                      "estimate", "mcse", "low", "high", "reps", "B",
                      "valid_replicates", "total_replicates", "numerical_failures")
missing_columns <- setdiff(required_columns, names(dat))
if (length(missing_columns)) {
  stop("Incomplete R revision summary; missing columns: ", paste(missing_columns, collapse = ", "))
}
if (!nrow(dat)) stop("R revision summary contains no results.")
for (nm in c("K", "estimate", "mcse", "low", "high", "reps", "B",
             "valid_replicates", "total_replicates", "numerical_failures")) {
  dat[[nm]] <- suppressWarnings(as.numeric(dat[[nm]]))
}
key <- paste(dat$scenario, dat$metric, sep = ":")
if (anyDuplicated(key)) stop("Duplicate scenario/metric rows in R revision summary.")
small_ids <- c("S10", "S20", "S40", "S80")
dependence_ids <- c("HN", "HA", "LN", "LA", "K8N", "K8A", "K12N", "K12A")
power_ids <- c("HA", "LA", "K8A", "K12A")
all_ids <- c(small_ids, dependence_ids)
fwer_metrics <- c("holm_fwer", "maxjel_fwer", "waldmax_fwer")
coverage_metrics <- c("maxjel_coverage", "bonf_coverage")
power_metrics <- c("holm_discovery", "maxjel_discovery", "waldmax_discovery")
required_keys <- c(as.vector(outer(all_ids, c(fwer_metrics, coverage_metrics), paste, sep = ":")),
                   as.vector(outer(power_ids, power_metrics, paste, sep = ":")))
missing_keys <- setdiff(required_keys, key)
if (length(missing_keys)) {
  stop("Incomplete revision metrics; no figures or tables were produced. Missing: ",
       paste(missing_keys, collapse = ", "))
}
used <- dat[match(required_keys, key), , drop = FALSE]
finite_columns <- c("estimate", "mcse", "low", "high", "valid_replicates",
                    "total_replicates", "numerical_failures")
if (any(!is.finite(as.matrix(used[, finite_columns])))) {
  stop("Nonfinite estimates, intervals, or replicate counts in required revision metrics.")
}
probability_tolerance <- 1e-12
if (any(used$estimate < -probability_tolerance | used$estimate > 1 + probability_tolerance |
        used$low < -probability_tolerance | used$high > 1 + probability_tolerance |
        used$low > used$estimate + probability_tolerance |
        used$high < used$estimate - probability_tolerance | used$mcse < 0 |
        used$valid_replicates < 1 | used$total_replicates < used$valid_replicates |
        used$numerical_failures < 0)) {
  stop("Invalid probabilities, Monte Carlo intervals, or replicate counts in revision summary.")
}
# At a zero or unit estimate, Wilson arithmetic can differ by machine epsilon.
# Only rounding within the above tolerance is repaired for plotting/formatting.
dat$low <- pmax(0, pmin(dat$low, dat$estimate))
dat$high <- pmin(1, pmax(dat$high, dat$estimate))
if (any(used$numerical_failures > 0)) {
  warning("Numerical failures are reported. Figures use valid replicate denominators; ",
          "the TeX table includes the valid/total counts.")
}
extract <- function(ids, metric) {
  dat[match(paste(ids, metric, sep = ":"), key), , drop = FALSE]
}
scenario_row <- function(id) dat[match(id, dat$scenario), , drop = FALSE]
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)
blue <- c("#08306B", "#2171B5", "#6BAED6")
method_labels <- c("Holm-JEL", "max-JEL", "Wald max-T")
coverage_labels <- c("max-JEL", "Bonferroni-JEL")

# Each figure is drawn independently on both devices so the PDF retains vector text.
save_figure <- function(stem, width, height, draw) {
  for (ext in c("pdf", "png")) {
    destination <- file.path(figures_dir, paste0(stem, suffix, ".", ext))
    if (ext == "pdf") {
      grDevices::pdf(destination, width = width, height = height,
                     family = "sans", useDingbats = FALSE, title = stem)
    } else {
      grDevices::png(destination, width = round(width * 300), height = round(height * 300),
                     res = 300, type = if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
    }
    tryCatch(draw(), finally = grDevices::dev.off())
    message("Wrote ", destination)
  }
}

panel <- function(ids, metrics, labels, ylab, title, nominal = NULL,
                  xlabels = ids, xlab = "Scenario", connect = FALSE,
                  coverage = FALSE, legend_columns = length(metrics), fixed_ylim = NULL) {
  rr <- lapply(metrics, function(m) extract(ids, m))
  values <- unlist(lapply(rr, function(x) c(x$low, x$high)))
  if (!is.null(fixed_ylim)) {
    ylim <- fixed_ylim
  } else if (coverage) {
    ylim <- c(max(0, min(c(values, .95)) - .04), 1.025)
  } else {
    upper <- min(1.025, max(c(values, .05)) * 1.23 + .015)
    ylim <- c(0, upper)
  }
  xx <- seq_along(ids)
  plot(xx, rr[[1L]]$estimate, type = "n", xlim = c(.55, length(ids) + .45), ylim = ylim,
       xaxt = "n", xlab = xlab, ylab = ylab, main = "", bty = "l",
       cex.main = 1.03, font.main = 2)
  mtext(title, side = 3, line = 2.6, font = 2, cex = 1.03)
  axis(1, at = xx, labels = xlabels)
  abline(h = pretty(ylim), col = "#E5EAF0", lwd = .6)
  if (!is.null(nominal)) abline(h = nominal, lty = 2, lwd = 1.1, col = "#4C5663")
  offsets <- seq(-.085, .085, length.out = length(metrics))
  for (j in seq_along(metrics)) {
    r <- rr[[j]]
    xpos <- xx + offsets[j]
    if (connect) lines(xpos, r$estimate, col = blue[j], lwd = 1.2)
    arrows(xpos, r$low, xpos, r$high, angle = 90, code = 3,
           length = .035, col = blue[j], lwd = 1.1)
    points(xpos, r$estimate, pch = c(16, 17, 15)[j], col = blue[j], cex = .85)
  }
  legend("top", inset = c(0, -.10), legend = labels, col = blue[seq_along(metrics)],
         pch = c(16, 17, 15)[seq_along(metrics)], lty = if (connect) 1 else 0,
         lwd = 1.1, bty = "n", xpd = NA, ncol = legend_columns, cex = .79,
         x.intersp = .65, y.intersp = .7)
  box(bty = "l")
}

save_figure("supp_smallest_class", 10.4, 4.15, function() {
  par(mfrow = c(1, 2), mar = c(4.6, 4.5, 4.0, .8), oma = c(1.0, 0, 0, 0),
      mgp = c(2.7, .7, 0), las = 1, cex = .9, col.axis = "#303846", col.lab = "#172033")
  class_sizes <- vapply(small_ids, function(id) {
    ns <- strsplit(as.character(scenario_row(id)$ns), "[/, ]+")[[1L]]
    ns[1L]
  }, character(1L))
  panel(small_ids, fwer_metrics, method_labels, "Complete-null FWER",
        "a  Familywise error", .05, class_sizes,
        "Class 1 size (class 2 = 40; class 3 = 160)", connect = TRUE)
  panel(small_ids, coverage_metrics, coverage_labels, "Joint 95% interval coverage",
        "b  Simultaneous coverage", .95, class_sizes,
        "Class 1 size (class 2 = 40; class 3 = 160)", connect = TRUE, coverage = TRUE)
  mtext("Bars show 95% Monte Carlo intervals. Dashed lines mark the nominal levels.",
        outer = TRUE, side = 1, line = -.1, cex = .77, col = "#4C5663")
})

save_figure("supp_dependence_scalability", 9.0, 6.8, function() {
  par(mfrow = c(2, 1), mar = c(3.6, 4.8, 3.9, .8), oma = c(2.3, 0, 0, 0),
      mgp = c(2.8, .65, 0), las = 1, cex = .9, col.axis = "#303846", col.lab = "#172033")
  panel(dependence_ids, fwer_metrics, method_labels, "True-null FWER",
        "a  Dependence and number of classifiers", .05, xlab = "")
  panel(dependence_ids, coverage_metrics, coverage_labels, "Joint 95% interval coverage",
        "b  Simultaneous coverage", .95, xlab = "Scenario", coverage = TRUE)
  mtext("H: heterogeneous Gaussian; L: nonlinear copula; K8/K12: 8/12 classifiers.",
        side = 1, outer = TRUE, line = .1, cex = .78, col = "#4C5663")
  mtext("N: complete null; A: sparse alternative. Bars: 95% Monte Carlo intervals.",
        side = 1, outer = TRUE, line = 1.15, cex = .78, col = "#4C5663")
})

save_figure("supp_power", 7.2, 4.2, function() {
  par(mar = c(4.4, 4.8, 4.2, .8), oma = c(1.15, 0, 0, 0),
      mgp = c(2.9, .7, 0), las = 1, cex = .95, col.axis = "#303846", col.lab = "#172033")
  panel(power_ids, power_metrics, method_labels, "Probability of any true discovery",
        "Multiplicity-adjusted power", fixed_ylim = c(0, 1.02))
  mtext("Bars show 95% Monte Carlo intervals; power is reported with true-null FWER in the companion figure.",
        side = 1, outer = TRUE, line = 0, cex = .72, col = "#4C5663")
})

# TeX table is supplementary source, deliberately independent of manuscript numbering.
tex_probability <- function(id, metric) {
  r <- extract(id, metric)
  sprintf("%.3f [%.3f, %.3f]", r$estimate, r$low, r$high)
}
tex_escape <- function(x) {
  x <- gsub("\\", "\\textbackslash{}", x, fixed = TRUE)
  for (ch in c("_", "%", "&", "#", "$")) x <- gsub(ch, paste0("\\", ch), x, fixed = TRUE)
  x
}
row_end <- paste0(" ", strrep("\\", 2L))
tex_lines <- c(
  "% Generated by R/08_revision_figures.R from revision_summary; base R only.",
  "\\begin{table*}[t]", "\\centering", "\\small",
  "\\caption{R-generated revision simulations. Entries are probabilities with 95\\% Monte Carlo confidence intervals. Coverage refers to simultaneous coverage of all pairwise contrasts. Valid/total denotes the number of usable simulation replications; discarded numerical failures, if any, must also be reported.}",
  "\\label{tab:R-revision-calibration}",
  "\\begin{tabular}{llrrrr}", "\\hline",
  paste0("Scenario & Class sizes & Valid/total & max-JEL FWER & max-JEL coverage & Bonferroni-JEL coverage", row_end),
  "\\hline"
)
for (id in all_ids) {
  r <- extract(id, "maxjel_fwer")
  tex_lines <- c(tex_lines, paste0(
    tex_escape(id), " & ", tex_escape(as.character(r$ns)), " & ",
    format(r$valid_replicates, scientific = FALSE, trim = TRUE), "/",
    format(r$total_replicates, scientific = FALSE, trim = TRUE), " & ",
    tex_probability(id, "maxjel_fwer"), " & ", tex_probability(id, "maxjel_coverage"),
    " & ", tex_probability(id, "bonf_coverage"), row_end))
}
tex_lines <- c(tex_lines, "\\hline", "\\end{tabular}", "\\end{table*}")

# Optional paired gains are computed by 07 from the same replicate-level indicators.
# They are never reconstructed by treating the methods as independent samples.
paired_file <- file.path(results_dir, paste0("revision_paired_power_gains", suffix, ".csv"))
if (file.exists(paired_file)) {
  paired <- utils::read.csv(paired_file, stringsAsFactors = FALSE)
  needed <- c("scenario", "gain", "low", "high", "reps", "failed_replicates")
  if (!all(needed %in% names(paired)) || anyDuplicated(paired$scenario) ||
      !all(power_ids %in% paired$scenario)) {
    stop("Incomplete paired-power file: ", paired_file)
  }
  paired <- paired[match(power_ids, paired$scenario), , drop = FALSE]
  if (any(!is.finite(as.matrix(paired[, setdiff(needed, "scenario")])))) {
    stop("Nonfinite paired-power intervals: ", paired_file)
  }
  tex_lines <- c(tex_lines, "", "\\begin{table*}[t]", "\\centering", "\\small",
    "\\caption{Multiplicity-adjusted power and paired power differences. Power is the probability of at least one true discovery. The max-JEL minus Holm-JEL difference uses the paired replicate indicators and a 95\\% Monte Carlo interval. These power values must be interpreted with the true-null FWER results.}",
    "\\label{tab:R-revision-power}", "\\begin{tabular}{lrrrr}", "\\hline",
    paste0("Scenario & Holm-JEL & max-JEL & Wald max-T & Paired difference", row_end), "\\hline")
  for (i in seq_along(power_ids)) {
    id <- power_ids[i]
    gain <- sprintf("%.3f [%.3f, %.3f]", paired$gain[i], paired$low[i], paired$high[i])
    tex_lines <- c(tex_lines, paste0(tex_escape(id), " & ",
      tex_probability(id, "holm_discovery"), " & ", tex_probability(id, "maxjel_discovery"),
      " & ", tex_probability(id, "waldmax_discovery"), " & ", gain, row_end))
  }
  tex_lines <- c(tex_lines, "\\hline", "\\end{tabular}", "\\end{table*}")
}

runtime_file <- file.path(results_dir, paste0("revision_controlled_runtime", suffix, ".csv"))
if (file.exists(runtime_file)) {
  timing <- utils::read.csv(runtime_file, stringsAsFactors = FALSE)
  needed <- c("K", "q", "n", "B", "repetitions", "median_ms", "p95_ms", "numerical_failures")
  if (!all(c(needed, "alternative") %in% names(timing)) || !nrow(timing) ||
      any(!is.finite(as.matrix(timing[, needed]))) ||
      any(!tolower(as.character(timing$alternative)) %in% c("true", "false", "1", "0")) ||
      any(timing$median_ms < 0 | timing$p95_ms < timing$median_ms | timing$repetitions < 1)) {
    stop("Incomplete or invalid controlled-runtime file: ", runtime_file)
  }
  tex_lines <- c(tex_lines, "", "\\begin{table*}[t]", "\\centering", "\\small",
    "\\caption{Controlled R timing, in milliseconds per simulated replicate. The reported median and 95th percentile apply to the recorded execution host and the stated number of timing repetitions.}",
    "\\label{tab:R-revision-runtime}", "\\begin{tabular}{lrrrrrrrr}", "\\hline",
    paste0("Setting & $K$ & $q$ & $n$ & $B$ & Repeats & Median ms & 95th pct. ms & Failures", row_end), "\\hline")
  for (i in seq_len(nrow(timing))) {
    r <- timing[i, ]
    setting <- if (tolower(as.character(r$alternative)) %in% c("true", "1")) "$H_1$" else "$H_0$"
    tex_lines <- c(tex_lines, paste0(paste(c(setting, r$K, r$q, r$n, r$B, r$repetitions,
      sprintf("%.2f", r$median_ms), sprintf("%.2f", r$p95_ms), r$numerical_failures),
      collapse = " & "), row_end))
  }
  tex_lines <- c(tex_lines, "\\hline", "\\end{tabular}", "\\end{table*}")
}
table_file <- file.path(results_dir, paste0("revision_tables", suffix, ".tex"))
writeLines(tex_lines, table_file)
message("Wrote ", table_file)
message("All supplementary figures were constructed exclusively from ", input)
