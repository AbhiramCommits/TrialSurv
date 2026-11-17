#!/usr/bin/env Rscript
# Step 1 analysis: Kaplan-Meier, Cox proportional hazards, PH diagnostics and
# held-out concordance on the colon cancer trial data.

suppressPackageStartupMessages({
  if (file.exists("DESCRIPTION") && requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(".", quiet = TRUE)
  } else {
    library(trialsurv)
  }
})
suppressPackageStartupMessages(library(dplyr))

set.seed(20261117)
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

dat <- load_trial_data()
cohort <- attr(dat, "cohort")

cat("\n=== TrialSurv: frequentist survival analysis ===\n")
cat("Cohort:", length(unique(dat$id)), "patients\n")
print(as.data.frame(cohort), row.names = FALSE)

write_result(
  dplyr::mutate(cohort, dataset = "survival::colon", time_unit = "months"),
  "cohort_summary"
)

# ---- Reproducible stratified 70/30 split of patients -----------------------
death <- dat[dat$endpoint == "death", , drop = FALSE]
set.seed(20261117)
key <- interaction(death$rx, death$status, drop = TRUE)
train_rows <- unlist(lapply(split(seq_len(nrow(death)), key), function(i) {
  sample(i, size = floor(0.7 * length(i)))
}), use.names = FALSE)
train_ids <- death$id[train_rows]
train <- dat[dat$id %in% train_ids, , drop = FALSE]
test <- dat[!dat$id %in% train_ids, , drop = FALSE]
test_death <- test[test$endpoint == "death", , drop = FALSE]

split_summary <- tibble::tibble(
  set = c("train", "test"),
  n_patients = c(length(unique(train$id)), length(unique(test$id))),
  n_death_events = c(
    sum(train$status[train$endpoint == "death"]),
    sum(test_death$status)
  )
)
write_result(split_summary, "split_summary")
cat("\nSplit:\n")
print(as.data.frame(split_summary), row.names = FALSE)

# ---- Kaplan-Meier -----------------------------------------------------------
km <- km_fit(dat, by = "rx", endpoint = "death")
write_result(km$summary, "km_summary")

cat("\nKaplan-Meier, death endpoint:\n")
print(as.data.frame(km$summary[, c(
  "arm", "n", "n_events", "median", "median_lower", "median_upper",
  "surv_12m", "surv_36m", "surv_60m", "logrank_p"
)]), row.names = FALSE, digits = 4)
cat(sprintf(
  "Log-rank: chi-square = %.3f on %d df, p = %.4g\n",
  km$logrank$chisq, km$logrank$df, km$logrank$p_value
))

plot_km(dat, by = "rx", endpoint = "death", name = "km_by_arm.png")

# ---- Cox models -------------------------------------------------------------
cx <- cox_fit(train, endpoint = "death")
write_result(cx$hr, "cox_hazard_ratios")

cat("\nCox hazard ratios (trained on the 70% split):\n")
print(as.data.frame(cx$hr[, c(
  "model", "label", "hr", "conf.low", "conf.high", "p.value"
)]), row.names = FALSE, digits = 4)

# ---- Proportional hazards diagnostics --------------------------------------
diag <- cox_diagnostics(cx)
write_result(diag, "ph_diagnostics")

cat("\nProportional hazards diagnostics (cox.zph, KM time transform):\n")
print(as.data.frame(diag[, c(
  "term", "statistic", "df", "p_value", "ph_violation", "slope", "slope_p"
)]), row.names = FALSE, digits = 4)
violators <- attr(diag, "violators")
cat(sprintf(
  "Covariates violating PH at alpha = 0.05: %s\n",
  if (length(violators) == 0L) "none" else paste(violators, collapse = ", ")
))
cat(sprintf(
  "Global cox.zph test: chi-square = %.3f, p = %.4g\n",
  diag$statistic[diag$term == "GLOBAL"], diag$p_value[diag$term == "GLOBAL"]
))
cat(sprintf(
  "Deviance residual outliers (|r| > 3): %d of %d (%.1f%%)\n",
  diag$deviance_outliers[1], cx$n, diag$deviance_outlier_pct[1]
))

# ---- Held-out concordance ---------------------------------------------------
cc_un <- concordance_eval(cx$unadjusted, test_death)
cc_adj <- concordance_eval(cx$adjusted, test_death)
cc <- dplyr::bind_rows(
  dplyr::mutate(cc_un, model = "unadjusted"),
  dplyr::mutate(cc_adj, model = "adjusted")
)
write_result(cc, "concordance_frequentist")

cat("\nHeld-out concordance (30% test set):\n")
print(as.data.frame(cc[, c(
  "model", "c_index", "std_error", "conf.low", "conf.high", "n", "n_events"
)]), row.names = FALSE, digits = 4)

plot_schoenfeld(cx$adjusted, name = "schoenfeld.png")

cat("\nWrote results/km_summary.csv, results/cox_hazard_ratios.csv,")
cat(" results/ph_diagnostics.csv, results/concordance_frequentist.csv\n")
cat("Wrote figures/km_by_arm.png and figures/schoenfeld.png\n")
cat("=== Done ===\n")
