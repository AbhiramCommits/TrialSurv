#!/usr/bin/env Rscript
# Run the complete TrialSurv pipeline in order and record reproducibility
# metadata. Each stage is a separate R process, exactly as a user would run it.

root <- getwd()
if (!file.exists(file.path(root, "DESCRIPTION"))) {
  stop("run_all.R must be run from the package root (the directory with DESCRIPTION).")
}

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

run_stage <- function(label, command, args = character()) {
  cat(sprintf("\n===== %s =====\n", label))
  t0 <- Sys.time()
  status <- system2(command, args = args)
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  if (!identical(status, 0L)) {
    stop(sprintf("Stage '%s' failed with status %s.", label, status), call. = FALSE)
  }
  cat(sprintf("Stage '%s' finished in %.1f seconds.\n", label, elapsed))
  data.frame(stage = label, seconds = elapsed, stringsAsFactors = FALSE)
}

timings <- list()
timings[[length(timings) + 1L]] <- run_stage(
  "01_frequentist", "Rscript", "scripts/01_frequentist.R"
)
timings[[length(timings) + 1L]] <- run_stage(
  "02_bayesian", "Rscript", "scripts/02_bayesian.R"
)
timings[[length(timings) + 1L]] <- run_stage(
  "03_simulation", "Rscript", "scripts/03_simulation.R"
)
timings[[length(timings) + 1L]] <- run_stage(
  "04_report", "Rscript",
  c("-e", shQuote("pkgload::load_all('.', quiet = TRUE); trialsurv::render_report()"))
)

timing_tbl <- do.call(rbind, timings)
timing_tbl$total <- sum(timing_tbl$seconds)

packages <- c(
  "survival", "survminer", "ggsurvfit", "brms", "rstan", "posterior",
  "dplyr", "tidyr", "ggplot2", "purrr", "testthat", "roxygen2", "quarto",
  "covr"
)
versions <- vapply(packages, function(p) {
  if (requireNamespace(p, quietly = TRUE)) {
    as.character(utils::packageVersion(p))
  } else {
    "not installed"
  }
}, character(1))

sink("results/session_info.txt")
cat("TrialSurv pipeline session information\n")
cat("=======================================\n\n")
cat("R version: ", R.version.string, "\n", sep = "")
cat("Platform: ", R.version$platform, "\n", sep = "")
cat("Run finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n", sep = "")
cat("Stage runtimes (seconds):\n")
print(timing_tbl, row.names = FALSE)
cat("\nTotal runtime (seconds): ", round(sum(timing_tbl$seconds), 1), "\n\n",
  sep = ""
)
cat("Package versions:\n")
for (p in names(versions)) cat(sprintf("  %-12s %s\n", p, versions[[p]]))
cat("\nFull sessionInfo():\n\n")
print(utils::sessionInfo())
sink()

timing_csv <- timing_tbl[, c("stage", "seconds")]
utils::write.csv(timing_csv, "results/stage_runtimes.csv", row.names = FALSE)

cat("\n===== Pipeline complete =====\n")
print(timing_tbl, row.names = FALSE)
cat(sprintf("Total runtime: %.1f seconds (%.1f minutes)\n",
  sum(timing_tbl$seconds), sum(timing_tbl$seconds) / 60
))
cat("Wrote results/session_info.txt and results/stage_runtimes.csv\n")
