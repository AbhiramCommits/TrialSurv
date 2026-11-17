test_that("a Cox model recovers a known hazard ratio", {
  set.seed(20261117)
  n <- 4000
  trt <- rep(0:1, length.out = n)
  hr_true <- 0.5
  rate <- 0.08
  event_time <- stats::rexp(n, rate = rate * ifelse(trt == 1, hr_true, 1))
  admin <- stats::quantile(event_time, 0.95)
  dat <- data.frame(
    id = seq_len(n),
    time = pmin(event_time, admin),
    status = as.integer(event_time <= admin),
    endpoint = factor("death", levels = c("recurrence", "death")),
    rx = factor(ifelse(trt == 1, "Lev+5FU", "Obs"), levels = c("Obs", "Lev", "Lev+5FU")),
    age = 60,
    sex = factor("Female", levels = c("Female", "Male")),
    nodes = 5,
    extent = factor("Serosa", levels = c("Submucosa", "Muscle", "Serosa", "Contiguous")),
    differ = factor("Moderate", levels = c("Well", "Moderate", "Poor")),
    obstruct = 0L
  )

  cx <- cox_fit(dat, covariates = character(0))
  hr <- cx$hr$hr[cx$hr$model == "unadjusted"][1]
  ci <- cx$hr[cx$hr$model == "unadjusted", ]

  expect_true(abs(hr - hr_true) <= 0.08)
  expect_true(ci$conf.low <= hr_true && ci$conf.high >= hr_true)
})

test_that("cox_diagnostics() has one row per coefficient plus global", {
  d <- load_trial_data()
  cx <- cox_fit(d, endpoint = "death")
  diag <- cox_diagnostics(cx)

  expect_s3_class(diag, "tbl_df")
  expect_equal(nrow(diag), length(stats::coef(cx$adjusted)) + 1L)
  expect_equal(sum(diag$term == "GLOBAL"), 1L)
  expect_true(all(c("statistic", "p_value", "ph_violation", "slope", "deviance_outliers") %in% names(diag)))
  expect_type(attr(diag, "violators"), "character")
  expect_true(all(diag$p_value >= 0 & diag$p_value <= 1))
})

test_that("concordance_eval() is in a sensible range for informative covariates", {
  set.seed(42)
  n <- 800
  x <- stats::rnorm(n)
  event_time <- stats::rexp(n, rate = 0.05 * exp(0.8 * x))
  admin <- stats::quantile(event_time, 0.9)
  dat <- data.frame(
    time = pmin(event_time, admin),
    status = as.integer(event_time <= admin),
    endpoint = factor("death", levels = c("recurrence", "death")),
    rx = factor("Obs", levels = c("Obs", "Lev", "Lev+5FU")),
    x = x
  )
  fit <- survival::coxph(survival::Surv(time, status) ~ x,
    data = dat, x = TRUE, y = TRUE
  )
  cc <- concordance_eval(fit, dat)

  expect_s3_class(cc, "tbl_df")
  expect_gt(cc$c_index, 0.5)
  expect_lt(cc$c_index, 1)
  expect_true(cc$conf.low <= cc$c_index && cc$c_index <= cc$conf.high)
})
