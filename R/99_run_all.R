# 99_run_all.R
# Base-R workflow driver. From a terminal:
#   Rscript R/99_run_all.R check
#   Rscript R/99_run_all.R quick
#   Rscript R/99_run_all.R full
# From RStudio (source defines the function without starting a long run):
#   source("R/99_run_all.R")
#   run_jel_workflow("check")
#   run_jel_workflow("full")
# Existing environment variables override all mode defaults.

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
.jel_driver_root <- .jel_find_root()

run_jel_workflow <- function(mode = "real", project_root = .jel_driver_root) {
  mode <- tolower(as.character(mode))
  allowed <- c("check", "quick", "real", "original", "revision", "full")
  if (length(mode) != 1L || !mode %in% allowed) {
    stop("Mode must be one of: ", paste(allowed, collapse = ", "))
  }
  project_root <- normalizePath(project_root, mustWork = TRUE)
  if (getRversion() < "3.6.0") stop("R >= 3.6.0 is required.")
  rscript <- file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript")
  if (!file.exists(rscript)) stop("Cannot find the Rscript executable belonging to this R installation: ", rscript)
  defaults <- c(
    RETUNE = "true", REAL_SEED = "2026", REAL_B_MULT = "20000",
    R_REPS = "5000", B_MULT = "2000", MC_CORES = "1", ALPHA = "0.05",
    SIM_SEED = "2026", INT_SEED = "2026", REV_REPS = "5000",
    REV_B_MULT = "1999", REV_SEED = "20261016", REV_RUNTIME_REPS = "500",
    REV_RESUME = "true", REV_CHECKPOINT_EVERY = "100", OUTPUT_SUFFIX = "_R"
  )
  if (mode == "check") defaults["OUTPUT_SUFFIX"] <- "_check"
  if (mode == "quick") {
    defaults[c("REAL_B_MULT", "R_REPS", "B_MULT", "REV_REPS", "REV_B_MULT",
               "REV_RUNTIME_REPS", "REV_CHECKPOINT_EVERY", "OUTPUT_SUFFIX")] <-
      c("199", "10", "199", "5", "199", "2", "5", "_quick")
    defaults <- c(defaults, INT_REPS = "5")
  }
  inherited <- Sys.getenv(names(defaults), unset = NA_character_)
  effective <- defaults
  effective[!is.na(inherited)] <- inherited[!is.na(inherited)]
  # INT_REPS is deliberately absent from original/full defaults: script 05's
  # original scenario-specific 500/1000 replication counts are preserved.
  if ("INT_REPS" %in% names(effective) && !nzchar(effective[["INT_REPS"]])) {
    effective <- effective[names(effective) != "INT_REPS"]
  }
  if (mode == "quick" && !tolower(effective[["RETUNE"]]) %in% c("1", "true", "yes")) {
    message("QUICK with RETUNE=false is a fixed-hyperparameter smoke check; CV losses will be NA.")
  }
  selected <- switch(mode,
    check = "00_self_check.R",
    real = c("00_self_check.R", "03_real_data_application.R"),
    original = c("00_self_check.R", "03_real_data_application.R",
                 "04_simulation_study.R", "05_interval_simulation.R", "06_make_figures.R"),
    revision = c("00_self_check.R", "07_revision_simulations.R", "08_revision_figures.R"),
    full = c("00_self_check.R", "03_real_data_application.R", "04_simulation_study.R",
             "05_interval_simulation.R", "06_make_figures.R",
             "07_revision_simulations.R", "08_revision_figures.R"),
    quick = c("00_self_check.R", "03_real_data_application.R", "04_simulation_study.R",
              "05_interval_simulation.R", "06_make_figures.R",
              "07_revision_simulations.R", "08_revision_figures.R")
  )
  # Independent deterministic numerical validation runs once per workflow.
  selected <- append(selected, "09_validate_revision.R", after = 1L)
  scripts <- file.path(project_root, "R", selected)
  if (any(!file.exists(scripts))) stop("Missing R script(s): ", paste(scripts[!file.exists(scripts)], collapse = ", "))
  # Check mode always parses every R file and loads 07 only under its no-run guard.
  run_id <- paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_", Sys.getpid())
  results_dir <- file.path(project_root, "results")
  log_dir <- file.path(results_dir, "logs", paste0(mode, "_", run_id))
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  restore_names <- unique(c(names(effective), "JEL_PROJECT_ROOT", "JEL_RUN_ID", "REVISION_HELPERS_ONLY"))
  old <- Sys.getenv(restore_names, unset = NA_character_)
  on.exit({
    if (any(is.na(old))) Sys.unsetenv(restore_names[is.na(old)])
    if (any(!is.na(old))) do.call(Sys.setenv, as.list(stats::setNames(old[!is.na(old)], restore_names[!is.na(old)])))
  }, add = TRUE)
  do.call(Sys.setenv, as.list(effective))
  Sys.setenv(JEL_PROJECT_ROOT = project_root, JEL_RUN_ID = run_id,
             REVISION_HELPERS_ONLY = "false")
  configuration <- c(effective, JEL_PROJECT_ROOT = project_root,
    JEL_RUN_ID = run_id, RUN_MODE = mode,
    INT_REPS_EFFECTIVE = Sys.getenv("INT_REPS", "scenario-specific 500/1000"),
    REV_SCENARIOS_EFFECTIVE = Sys.getenv("REV_SCENARIOS", "all 12"))
  utils::write.table(data.frame(variable = names(configuration), value = unname(configuration)),
    file.path(log_dir, "configuration.tsv"), sep = "\t", quote = TRUE, row.names = FALSE)
  writeLines(capture.output(sessionInfo()), file.path(log_dir, "session_info.txt"))
  source_files <- list.files(file.path(project_root, "R"), pattern = "[.]R$", full.names = TRUE)
  hashes <- tools::md5sum(source_files)
  utils::write.csv(data.frame(file = basename(names(hashes)), md5 = unname(hashes)),
    file.path(log_dir, "source_md5.csv"), row.names = FALSE)
  record <- data.frame(stage = character(), script = character(), start = character(),
    elapsed_seconds = numeric(), exit_status = integer(), log = character(), stringsAsFactors = FALSE)
  message("Mode: ", mode, "; serial workers by default; MC_CORES=", effective[["MC_CORES"]])
  message("Output suffix: '", effective[["OUTPUT_SUFFIX"]], "'. Logs: ", log_dir)
  message("Effective counts: R_REPS=", effective[["R_REPS"]],
    "; INT_REPS=", Sys.getenv("INT_REPS", "scenario-specific 500/1000"),
    "; B_MULT=", effective[["B_MULT"]], "; REAL_B_MULT=", effective[["REAL_B_MULT"]],
    "; REV_REPS=", effective[["REV_REPS"]], "; REV_B_MULT=", effective[["REV_B_MULT"]],
    "; REV_RUNTIME_REPS=", effective[["REV_RUNTIME_REPS"]])
  if (mode == "quick") message("QUICK is a low-replication pipeline check; its simulation rates are not scientific estimates.")
  if (mode %in% c("full", "original", "revision")) {
    message("This mode can take substantial time. Each child writes a stage log; inspect those files for progress.")
  }
  for (i in seq_along(selected)) {
    stage_log <- file.path(log_dir, sprintf("%02d_%s.log", i, sub("[.]R$", "", selected[i])))
    started <- Sys.time()
    message("[", i, "/", length(selected), "] Running ", selected[i], "\n  Log: ", stage_log)
    # system2 quotes its executable; explicitly quote the script argument so
    # project paths containing spaces work on macOS, Linux and Windows.
    status <- system2(command = rscript, args = c("--vanilla", shQuote(scripts[i])),
                      stdout = stage_log, stderr = stage_log, wait = TRUE)
    elapsed <- as.numeric(difftime(Sys.time(), started, units = "secs"))
    if (is.null(status)) status <- 0L
    record <- rbind(record, data.frame(stage = as.character(i), script = selected[i],
      start = format(started, "%Y-%m-%dT%H:%M:%S%z"), elapsed_seconds = elapsed,
      exit_status = as.integer(status), log = stage_log, stringsAsFactors = FALSE))
    utils::write.csv(record, file.path(log_dir, "stage_manifest.csv"), row.names = FALSE)
    output <- if (file.exists(stage_log)) readLines(stage_log, warn = FALSE) else character()
    if (status != 0L) {
      if (length(output)) cat(paste(tail(output, 40L), collapse = "\n"), "\n")
      stop(selected[i], " failed with exit status ", status, ". Full log: ", stage_log)
    }
    if (length(output)) cat(paste(tail(output, 5L), collapse = "\n"), "\n")
    message("  Completed in ", sprintf("%.2f", elapsed), " seconds.")
  }
  writeLines(c(sprintf("Mode %s completed successfully.", mode),
    sprintf("Run ID: %s", run_id), sprintf("Completed: %s", format(Sys.time())),
    "New R runs must be used to regenerate manuscript tables; they need not match prior Python/reference results."),
    file.path(log_dir, "COMPLETED.txt"))
  message("Workflow completed. Results: ", results_dir, "; logs: ", log_dir)
  invisible(list(mode = mode, run_id = run_id, configuration = configuration,
                 stages = record, log_directory = log_dir))
}

# Sourcing is safe in RStudio: it defines run_jel_workflow without launching it.
if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  selected_mode <- if (length(args)) args[1L] else Sys.getenv("RUN_MODE", "real")
  run_jel_workflow(selected_mode)
}
