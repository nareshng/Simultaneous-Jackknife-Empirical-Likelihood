# 99_run_all.R
# Reproduction driver.
#
#   Rscript R/99_run_all.R real   # real-data application only
#   Rscript R/99_run_all.R full   # manuscript-scale simulation (HPC recommended)

args <- commandArgs(trailingOnly=TRUE)
mode <- if (length(args)) tolower(args[1L]) else "real"
if (!mode %in% c("real","quick","full")) stop("Mode must be real, quick, or full.")

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
package_root <- normalizePath(file.path(dirname(script_file), ".."))
rscript <- file.path(R.home("bin"),"Rscript")

run <- function(script,env=character()) {
  message("\n--- Running ",script," ---")
  status <- system2(
    rscript,
    args=file.path(package_root,"R",script),
    env=env,
    stdout="",stderr=""
  )
  if (!identical(status,0L)) stop(script," failed with exit status ",status)
}

run("00_self_check.R")

if (mode=="real") {
  run("03_real_data_application.R",
      env=c("RETUNE=true","REAL_B_MULT=19999","OUTPUT_SUFFIX=_R"))
} else if (mode=="quick") {
  common <- c(
    "RETUNE=false","REAL_B_MULT=399","R_REPS=200","B_MULT=399",
    "INT_REPS=50","MC_CORES=1","OUTPUT_SUFFIX=_quick"
  )
  run("03_real_data_application.R",env=common)
  run("04_simulation_study.R",env=common)
  run("05_interval_simulation.R",env=common)
  run("06_make_figures.R",env=common)
} else {
  common <- c(
    "RETUNE=true","REAL_B_MULT=19999","R_REPS=5000","B_MULT=1999",
    paste0("MC_CORES=",Sys.getenv("MC_CORES","1")),"OUTPUT_SUFFIX=_R"
  )
  run("03_real_data_application.R",env=common)
  run("04_simulation_study.R",env=common)
  run("05_interval_simulation.R",env=common)
  run("06_make_figures.R",env=common)
}
message("\nReproduction workflow completed in mode: ",mode)
