# 00_self_check.R
# Lightweight implementation and data-integrity checks.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
if (is.na(script_file) || !nzchar(script_file)) {
  package_root <- normalizePath(getwd())
} else {
  package_root <- normalizePath(file.path(dirname(script_file), ".."))
}
source(file.path(package_root,"R","01_jel_core.R"))
jel_self_check()

dat <- utils::read.csv(
  file.path(package_root,"data","urinary_biomarkers_pancreatic_cancer.csv"),
  stringsAsFactors=FALSE,check.names=FALSE
)
stopifnot(nrow(dat)==590L)
stopifnot(sum(dat$patient_cohort==1)==332L)
stopifnot(sum(dat$patient_cohort==2)==258L)
tab <- table(dat$diagnosis[dat$patient_cohort==2])
stopifnot(identical(as.integer(tab),c(102L,119L,37L)))
stopifnot(!anyNA(dat[,c("age","creatinine","LYVE1","REG1B","TFF1")]))
cat("Core pseudo-value and data-integrity checks passed.\n")
