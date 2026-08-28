# 06_make_figures.R
# Figures from CSV result files using base R.

script_file <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1L])
if (is.na(script_file) || !nzchar(script_file)) {
  package_root <- normalizePath(getwd())
} else {
  package_root <- normalizePath(file.path(dirname(script_file), ".."))
}
results_dir <- file.path(package_root, "results")
figures_dir <- file.path(package_root, "figures")
dir.create(figures_dir, recursive=TRUE, showWarnings=FALSE)

use_file <- function(stem) {
  candidates <- c(
    file.path(results_dir,paste0(stem,Sys.getenv("OUTPUT_SUFFIX","_R"),".csv")),
    file.path(results_dir,paste0(stem,".csv"))
  )
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("No result file found for ",stem)
  hit[1L]
}

# Fig. 1: null type-I error and FWER.
null <- utils::read.csv(use_file("simulation_null_results"),check.names=FALSE)
methods <- c("Global_JEL","Global_Wald","Unadjusted_any","Holm_JEL","Max_JEL","Wald_maxT")
labels <- c("Global JEL","Global Wald","Unadjusted pairwise","Holm-JEL","max-JEL","Wald maxT")
png(file.path(figures_dir,"fig1_null_fwer_R.png"),width=2200,height=1250,res=180)
matplot(
  seq_len(nrow(null)),as.matrix(null[,methods]),type="b",lty=1,pch=seq_along(methods),
  xaxt="n",xlab="Null scenario",ylab="Rejection probability",
  ylim=c(0,max(.26,null[,methods],na.rm=TRUE)),
  main="Type-I error and familywise error under the global null"
)
axis(1,at=seq_len(nrow(null)),labels=null$ScenarioID)
abline(h=.05,lty=2,lwd=2)
legend("topleft",legend=labels,lty=1,pch=seq_along(methods),cex=.82,ncol=2,bty="n")
dev.off()

# Fig. 2: power for alternative scenarios.
alt <- utils::read.csv(use_file("simulation_alternative_results"),check.names=FALSE)
methods2 <- c("Holm_any_true","Max_JEL_any_true","Wald_maxT_any_true")
labels2 <- c("Holm-JEL","max-JEL","Wald maxT")
png(file.path(figures_dir,"fig2_alternative_power_R.png"),width=2200,height=1250,res=180)
matplot(
  seq_len(nrow(alt)),as.matrix(alt[,methods2]),type="b",lty=1,pch=1:3,
  xaxt="n",xlab="Alternative scenario",ylab="Probability of at least one true discovery",
  ylim=c(0,1),main="Multiplicity-adjusted power under alternatives"
)
axis(1,at=seq_len(nrow(alt)),labels=alt$ScenarioID)
legend("topleft",legend=labels2,lty=1,pch=1:3,bty="n")
dev.off()

# Fig. 3: real-data VUS.
perf <- utils::read.csv(use_file("real_data_model_performance"),check.names=FALSE)
png(file.path(figures_dir,"fig3_real_vus_R.png"),width=1800,height=1100,res=180)
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
png(file.path(figures_dir,"fig4_pairwise_intervals_R.png"),width=2100,height=1350,res=180)
op <- par(mar=c(5,16,3,1))
xlim <- range(c(pair[[lo_name]],pair[[hi_name]],0))
plot(pair[[diff_name]][ord],seq_along(ord),xlim=xlim+c(-.02,.02),
     ylim=c(.5,length(ord)+.5),yaxt="n",ylab="",
     xlab="VUS difference (A - B)",pch=19,
     main="95% simultaneous max-JEL intervals")
segments(pair[[lo_name]][ord],seq_along(ord),
         pair[[hi_name]][ord],seq_along(ord),lwd=2)
axis(2,at=seq_along(ord),labels=labels4[ord],las=2,cex.axis=.75)
abline(v=0,lty=2)
par(op)
dev.off()
message("Figures written to ",figures_dir)
