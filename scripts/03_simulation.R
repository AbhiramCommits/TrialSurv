#!/usr/bin/env Rscript
# Step 3 analysis: Monte Carlo power curves and type I error for a two-arm
# time-to-event trial, plus simulated versus Schoenfeld sample sizes.

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
cores <- max(1L, min(8L, parallel::detectCores() - 1L))

cat("\n=== TrialSurv: trial design simulation ===\n")
cat(sprintf("Using %d cores; 1000 simulations per design cell.\n", cores))

# ---- Main power grid --------------------------------------------------------
grid <- expand.grid(
  hr = c(1.0, 0.9, 0.8, 0.7, 0.6),
  n_per_arm = c(150, 250, 400, 600),
  dropout_rate = c(0, 0.10, 0.20),
  KEEP.OUT.ATTRS = FALSE
)
grid$accrual_months <- 24
grid$followup_months <- 60
grid$median_control <- 60

t0 <- Sys.time()
power <- trial_power(
  grid, n_sims = 1000, test = "cox", cores = cores, seed = 20261117
)
grid_runtime <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
write_result(power, "power_curve")

cat(sprintf(
  "Power grid: %d cells x 1000 simulations = %s trials in %.1f s\n",
  nrow(power), format(nrow(power) * 1000, big.mark = ","), grid_runtime
))

# ---- Empirical type I error at HR = 1 --------------------------------------
type1_cells <- power[power$hr == 1, , drop = FALSE]
type1_cells$scope <- "cell"
type1_cells$within_mc_error <- abs(type1_cells$rejection_rate - 0.05) <=
  2 * type1_cells$mc_se

n_pooled <- sum(type1_cells$n_sims)
x_pooled <- sum(type1_cells$rejection_rate * type1_cells$n_sims)
p_pooled <- x_pooled / n_pooled
se_pooled <- sqrt(p_pooled * (1 - p_pooled) / n_pooled)
z <- stats::qnorm(0.975)
pooled_row <- tibble::tibble(
  n_per_arm = NA_integer_,
  hr = 1,
  accrual_months = 24,
  followup_months = 60,
  dropout_rate = NA_real_,
  median_control = 60,
  test = "cox",
  alpha = 0.05,
  n_sims = as.integer(n_pooled),
  rejection_rate = p_pooled,
  mc_se = se_pooled,
  mc_lower = max(0, p_pooled - z * se_pooled),
  mc_upper = min(1, p_pooled + z * se_pooled),
  mean_events = NA_real_,
  mean_hr = NA_real_,
  mean_log_hr = NA_real_,
  scope = "pooled",
  within_mc_error = abs(p_pooled - 0.05) <= 2 * se_pooled
)
type1 <- dplyr::bind_rows(pooled_row, type1_cells)
write_result(type1, "type1_error")

cat("\nEmpirical type I error at HR = 1 (alpha = 0.05):\n")
print(as.data.frame(type1[, c(
  "scope", "n_per_arm", "dropout_rate", "n_sims", "rejection_rate",
  "mc_se", "mc_lower", "mc_upper", "within_mc_error"
)]), row.names = FALSE, digits = 4)

if (!isTRUE(pooled_row$within_mc_error)) {
  stop(sprintf(
    "Pooled empirical type I error %.4f is outside two Monte Carlo SEs of 0.05 (SE = %.4f).",
    p_pooled, se_pooled
  ), call. = FALSE)
}
cat("Type I error check passed: pooled empirical alpha is within two MC SEs of 0.05.\n")

# ---- Required sample size versus Schoenfeld --------------------------------
scenarios <- expand.grid(
  hr = c(0.8, 0.7),
  dropout_rate = c(0, 0.20),
  KEEP.OUT.ATTRS = FALSE
)
n_range <- c(100, 150, 200, 250, 300, 400, 500, 600, 800)

t1 <- Sys.time()
req <- purrr::map_dfr(seq_len(nrow(scenarios)), function(i) {
  required_n(
    target_power = 0.8,
    hr = scenarios$hr[i],
    n_range = n_range,
    n_sims = 1000,
    test = "cox",
    accrual_months = 24,
    followup_months = 60,
    dropout_rate = scenarios$dropout_rate[i],
    median_control = 60,
    alpha = 0.05,
    cores = cores,
    seed = 20261117
  )$summary
})
req_runtime <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
write_result(req, "sample_size_vs_schoenfeld")

cat("\nRequired per-arm n for 80% power (simulation vs Schoenfeld):\n")
print(as.data.frame(req[, c(
  "hr", "dropout_rate", "n_simulated", "power_at_simulated",
  "n_schoenfeld", "ratio_sim_over_schoenfeld"
)]), row.names = FALSE, digits = 3)

# ---- Power curves -----------------------------------------------------------
power_plot <- power |>
  dplyr::mutate(
    dropout_label = factor(
      sprintf("Dropout %d%%", round(100 * .data$dropout_rate)),
      levels = sprintf("Dropout %d%%", round(100 * sort(unique(.data$dropout_rate))))
    ),
    hr_label = factor(sprintf("HR = %.1f", .data$hr))
  )

p <- ggplot2::ggplot(
  power_plot,
  ggplot2::aes(
    x = .data$n_per_arm, y = .data$rejection_rate,
    colour = .data$hr_label, group = .data$hr_label
  )
) +
  ggplot2::geom_hline(yintercept = 0.8, linetype = 2, colour = "grey40") +
  ggplot2::geom_errorbar(
    ggplot2::aes(ymin = .data$mc_lower, ymax = .data$mc_upper),
    width = 8, alpha = 0.35
  ) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_point(size = 2.2) +
  ggplot2::facet_wrap(~dropout_label) +
  ggplot2::scale_y_continuous(
    limits = c(0, 1), breaks = seq(0, 1, 0.2),
    labels = function(x) sprintf("%.0f%%", 100 * x)
  ) +
  ggplot2::labs(
    title = "Power versus per-arm sample size",
    subtitle = "Two-arm trial, 24-month accrual, 60-month follow-up; dashed line is 80% power",
    x = "Patients per arm",
    y = "Empirical rejection rate (power)",
    colour = NULL,
    caption = "Bars are 95% Wilson intervals for the Monte Carlo rejection rate."
  ) +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(legend.position = "top")

fig_save(p, "power_curves.png", width = 10, height = 6)

# ---- Simulation metadata ----------------------------------------------------
total_sims <- nrow(power) * 1000 + nrow(scenarios) * length(n_range) * 1000
meta <- tibble::tibble(
  n_cells_power = nrow(power),
  n_sims_per_cell = 1000L,
  n_cells_required_n = nrow(req),
  total_simulations = as.integer(total_sims),
  cores = cores,
  runtime_power_sec = grid_runtime,
  runtime_required_n_sec = req_runtime,
  runtime_total_sec = grid_runtime + req_runtime,
  seed = 20261117L,
  test = "cox",
  alpha = 0.05,
  pooled_type1 = p_pooled,
  pooled_type1_mc_se = se_pooled
)
write_result(meta, "simulation_meta")

cat(sprintf(
  "\nTotal simulations: %s (%.1f s total runtime)\n",
  format(total_sims, big.mark = ","), grid_runtime + req_runtime
))
cat("Wrote results/power_curve.csv, results/type1_error.csv,")
cat(" results/sample_size_vs_schoenfeld.csv, results/simulation_meta.csv\n")
cat("Wrote figures/power_curves.png\n")
cat("=== Done ===\n")
