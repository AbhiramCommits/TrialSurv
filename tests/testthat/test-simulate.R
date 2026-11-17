test_that("simulate_trial() respects the design and produces censoring", {
  d <- simulate_trial(
    n_per_arm = 100, hr = 0.7, accrual_months = 12,
    followup_months = 6, dropout_rate = 0, median_control = 24, seed = 1
  )
  expect_s3_class(d, "tbl_df")
  expect_equal(nrow(d), 200L)
  expect_equal(sum(d$arm == "Control"), 100L)
  expect_equal(sum(d$arm == "Treatment"), 100L)
  expect_true(all(d$time > 0))
  expect_true(any(d$status == 0L))
  expect_true(all(d$status %in% c(0L, 1L)))
  expect_true(all(d$time <= d$admin_time + 1e-8))
})

test_that("event rate increases as the hazard ratio decreases", {
  d1 <- simulate_trial(2000,
    hr = 1.0, accrual_months = 24,
    followup_months = 24, dropout_rate = 0, median_control = 24, seed = 2
  )
  d2 <- simulate_trial(2000,
    hr = 0.5, accrual_months = 24,
    followup_months = 24, dropout_rate = 0, median_control = 24, seed = 2
  )
  expect_gt(mean(d1$status), mean(d2$status))
})

test_that("trial_power() returns a type I error consistent with alpha", {
  grid <- data.frame(
    n_per_arm = 300, hr = 1, dropout_rate = 0,
    accrual_months = 24, followup_months = 60, median_control = 60
  )
  p <- trial_power(grid, n_sims = 1000, test = "logrank", cores = 2, seed = 1)

  expect_s3_class(p, "tbl_df")
  expect_equal(nrow(p), 1L)
  expect_true(all(c(
    "rejection_rate", "mc_se", "mc_lower", "mc_upper",
    "mean_events", "mean_hr"
  ) %in% names(p)))
  expect_lte(abs(p$rejection_rate - 0.05), 3 * p$mc_se + 0.005)
  expect_gt(p$mean_events, 0)
  expect_gt(p$mean_hr, 0)
})

test_that("power is monotone in n and decreasing in HR", {
  grid <- expand.grid(
    n_per_arm = c(150, 400), hr = c(1.0, 0.6),
    dropout_rate = 0, KEEP.OUT.ATTRS = FALSE
  )
  grid$accrual_months <- 24
  grid$followup_months <- 60
  grid$median_control <- 60
  p <- trial_power(grid, n_sims = 600, test = "cox", cores = 2, seed = 3)

  pow_150 <- p$rejection_rate[p$n_per_arm == 150 & p$hr == 0.6]
  pow_400 <- p$rejection_rate[p$n_per_arm == 400 & p$hr == 0.6]
  pow_null <- p$rejection_rate[p$n_per_arm == 400 & p$hr == 1.0]
  expect_gt(pow_400, pow_150)
  expect_gt(pow_400, pow_null)
})

test_that("required_n() compares the simulation with Schoenfeld", {
  r <- required_n(
    target_power = 0.8, hr = 0.7,
    n_range = c(100, 200, 300), n_sims = 200,
    test = "logrank", cores = 2, seed = 4
  )
  expect_s3_class(r, "trialsurv_required_n")
  expect_s3_class(r$summary, "tbl_df")
  expect_true(all(c(
    "n_simulated", "n_schoenfeld", "power_at_simulated",
    "ratio_sim_over_schoenfeld"
  ) %in% names(r$summary)))
  expect_gt(r$summary$n_schoenfeld, 0)
  expect_equal(nrow(r$search), 3L)
})
