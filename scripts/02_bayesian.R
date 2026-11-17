#!/usr/bin/env Rscript
# Step 2 analysis: hierarchical Bayesian Weibull model with partial borrowing
# across nodal-risk strata, compared with the frequentist Cox model.

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
if (requireNamespace("rstan", quietly = TRUE)) {
  rstan::rstan_options(auto_write = TRUE)
}

cat("\n=== TrialSurv: Bayesian hierarchical survival analysis ===\n")
dat <- load_trial_data()
sp <- trial_split(dat, prop = 0.7, seed = 20261117, endpoint = "death")
train <- sp$train
test <- sp$test
cat(sprintf(
  "Train: %d patients, %d deaths; test: %d patients, %d deaths\n",
  length(unique(train$id)),
  sum(train$status[train$endpoint == "death"]),
  length(unique(test$id)),
  sum(test$status[test$endpoint == "death"])
))

# ---- Fit the hierarchical model --------------------------------------------
cat("\nFitting the brms/rstan hierarchical Weibull model (4000 iterations, 4 chains)...\n")
t0 <- Sys.time()
fit <- bayes_surv_fit(
  train,
  iter = 4000, chains = 4, seed = 20261117, backend = "brms",
  prior_scale = 0.5, adapt_delta = 0.95
)
fit_seconds <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
bs <- bayes_summary(fit)

max_rhat <- suppressWarnings(max(bs$rhat, na.rm = TRUE))
min_ess <- suppressWarnings(min(bs$ess_bulk, na.rm = TRUE))
min_ess_tail <- suppressWarnings(min(bs$ess_tail, na.rm = TRUE))

cat(sprintf(
  "First fit: %.1f s, max R-hat = %.4f, min bulk ESS = %.1f, divergences = %d\n",
  fit_seconds, max_rhat, min_ess, fit$n_divergent
))

# One more conservative re-fit if the sampler still shows any problem. brms
# already uses a non-centred parameterisation; this second pass raises the
# target acceptance rate and doubles the post-warm-up sample.
if (max_rhat > 1.01 || fit$n_divergent > 0) {
  message("Convergence check failed; refitting with adapt_delta = 0.97 and 6000 iterations.")
  t0 <- Sys.time()
  fit <- bayes_surv_fit(
    train,
    iter = 6000, chains = 4, seed = 20261117, backend = "brms",
    prior_scale = 0.5, adapt_delta = 0.97
  )
  fit_seconds <- fit_seconds + as.numeric(difftime(Sys.time(), t0, units = "secs"))
  bs <- bayes_summary(fit)
  max_rhat <- suppressWarnings(max(bs$rhat, na.rm = TRUE))
  min_ess <- suppressWarnings(min(bs$ess_bulk, na.rm = TRUE))
  min_ess_tail <- suppressWarnings(min(bs$ess_tail, na.rm = TRUE))
  cat(sprintf(
    "Second fit: max R-hat = %.4f, min bulk ESS = %.1f, divergences = %d\n",
    max_rhat, min_ess, fit$n_divergent
  ))
}

if (max_rhat > 1.01 || fit$n_divergent > 0) {
  stop(sprintf(
    paste0(
      "Bayesian fit did not converge after the conservative re-fit: ",
      "max R-hat = %.4f (limit 1.01), divergences = %d. ",
      "Do not report these posterior summaries."
    ),
    max_rhat, fit$n_divergent
  ), call. = FALSE)
}
cat("Convergence checks passed: max R-hat <= 1.01 and zero divergent transitions.\n")

# ---- Posterior summaries and shrinkage -------------------------------------
write_result(bs, "bayes_posterior_summary")

cx <- cox_fit(train, endpoint = "death")
sh <- shrinkage_report(fit, cx)
write_result(sh, "shrinkage_report")

cat("\nGlobal posterior hazard ratios (hazard scale):\n")
hr_global <- bs[bs$type == "hazard_ratio" & grepl("_global$", bs$parameter), ,
  drop = FALSE
]
print(as.data.frame(hr_global[, c(
  "parameter", "mean", "q2.5", "q97.5", "rhat", "ess_bulk"
)]), row.names = FALSE, digits = 3)

cat("\nSubgroup hazard ratios after partial borrowing:\n")
hr_sub <- bs[bs$type == "hazard_ratio" & !grepl("_global$", bs$parameter), ,
  drop = FALSE
]
print(as.data.frame(hr_sub[, c(
  "parameter", "mean", "q2.5", "q97.5", "rhat", "ess_bulk"
)]), row.names = FALSE, digits = 3)

cat("\nShrinkage report:\n")
print(as.data.frame(sh[, c(
  "arm", "stratum", "standalone_log_hr", "pooled_log_hr", "global_log_hr",
  "pct_toward_global", "pct_width_reduction"
)]), row.names = FALSE, digits = 3)

# ---- Model comparison and held-out concordance -----------------------------
cmp <- compare_models(cx, fit, test)
write_result(cmp, "model_comparison")

cat("\nFrequentist vs Bayesian comparison:\n")
print(as.data.frame(cmp[, c(
  "model", "parameter", "estimate", "conf.low", "conf.high",
  "interval_width", "c_index", "c_index_low", "c_index_high", "waic"
)]), row.names = FALSE, digits = 3)

# ---- Convergence artefacts and figures -------------------------------------
conv <- tibble::tibble(
  backend = fit$backend,
  formula = paste(deparse(fit$formula), collapse = " "),
  family = fit$family,
  iterations_per_chain = fit$iter,
  chains = fit$chains,
  prior_scale = fit$prior_scale,
  max_rhat = max_rhat,
  min_ess_bulk = min_ess,
  min_ess_tail = min_ess_tail,
  n_divergent = fit$n_divergent,
  fit_seconds = fit_seconds,
  seed = fit$seed
)
write_result(conv, "bayes_convergence")

meta <- list(
  backend = fit$backend,
  family = fit$family,
  formula = paste(deparse(fit$formula), collapse = " "),
  subgroup = fit$subgroup,
  stratum = fit$stratum,
  prior_scale = fit$prior_scale,
  iter = fit$iter,
  chains = fit$chains,
  seed = fit$seed,
  adapt_delta_attempts = as.list(fit$attempts),
  max_rhat = max_rhat,
  min_ess_bulk = min_ess,
  n_divergent = fit$n_divergent
)
write_result(meta, "bayes_fit_meta", format = "json")

plot_posterior_hr(fit, name = "posterior_hr.png")
plot_trace(fit, name = "trace.png")

cat("\nWrote results/bayes_posterior_summary.csv, results/shrinkage_report.csv,")
cat(" results/model_comparison.csv, results/bayes_convergence.csv\n")
cat("Wrote figures/posterior_hr.png and figures/trace.png\n")
cat("=== Done ===\n")
