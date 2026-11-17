# Helper: simulate a tidy colon-like trial with subgroup-varying treatment effects.
.sim_bayes_data <- function(n = 240, seed = 1, tau = 0.5, beta_lev = 0.15,
                            beta_f5 = -0.30) {
  set.seed(seed)
  rx <- factor(sample(c("Obs", "Lev", "Lev+5FU"), n, TRUE),
    levels = c("Obs", "Lev", "Lev+5FU")
  )
  nodes <- sample(0:12, n, TRUE)
  stratum <- cut(nodes,
    breaks = c(-Inf, 1, 3, 7, Inf),
    labels = c("n01", "n23", "n47", "n8p"), include.lowest = TRUE
  )
  u <- matrix(stats::rnorm(4 * 2, 0, tau), nrow = 4, ncol = 2)
  eta <- -5 + 0.12 * nodes +
    (rx == "Lev") * (beta_lev + u[as.integer(stratum), 1]) +
    (rx == "Lev+5FU") * (beta_f5 + u[as.integer(stratum), 2])
  k <- 1.2
  tt <- ((-log(stats::runif(n))) / exp(eta))^(1 / k)
  cens <- stats::rexp(n, 1 / 200)
  data.frame(
    id = seq_len(n),
    time = pmin(tt, cens),
    status = as.integer(tt <= cens),
    endpoint = factor("death", levels = c("recurrence", "death")),
    rx = rx,
    sex = factor(sample(c("Female", "Male"), n, TRUE),
      levels = c("Female", "Male")
    ),
    age = stats::rnorm(n, 60, 8),
    obstruct = rbinom(n, 1, 0.3),
    perfor = rbinom(n, 1, 0.02),
    adhere = rbinom(n, 1, 0.1),
    nodes = nodes,
    differ = factor(sample(c("Well", "Moderate", "Poor"), n, TRUE,
      prob = c(0.1, 0.7, 0.2)
    ), levels = c("Well", "Moderate", "Poor")),
    extent = factor(sample(c("Submucosa", "Muscle", "Serosa", "Contiguous"),
      n, TRUE,
      prob = c(0.05, 0.15, 0.7, 0.1)
    ), levels = c("Submucosa", "Muscle", "Serosa", "Contiguous")),
    surg = stats::runif(n, 0, 2),
    node4 = as.integer(nodes > 4)
  )
}

test_that("the fallback sampler fits the hierarchical model", {
  skip_if_not_installed("posterior")
  d <- .sim_bayes_data(n = 240, seed = 2, tau = 0.5)

  fit <- bayes_surv_fit(
    d,
    iter = 1000, chains = 2, seed = 11, backend = "mcmc",
    cores = 2, prior_scale = 0.5
  )
  expect_s3_class(fit, "trialsurv_bayes")
  expect_identical(fit$backend, "mcmc")
  expect_equal(fit$n_divergent, 0L)

  bs <- bayes_summary(fit)
  expect_s3_class(bs, "tbl_df")
  expect_true(all(c(
    "parameter", "type", "mean", "q2.5", "q97.5", "rhat",
    "ess_bulk", "ess_tail", "hazard_ratio"
  ) %in% names(bs)))
  expect_true(any(bs$type == "hazard_ratio"))
  expect_true(all(is.finite(bs$rhat[bs$parameter %in% c("b_lev", "b_f5", "shape")])))
  expect_true(all(bs$ess_bulk[bs$parameter %in% c("b_lev", "b_f5", "shape")] > 0))
  expect_true(is.finite(attr(bs, "max_rhat")))

  p_hr <- plot_posterior_hr(fit, save = FALSE)
  expect_s3_class(p_hr, "ggplot")
  p_trace <- plot_trace(fit, save = FALSE)
  expect_s3_class(p_trace, "ggplot")
})

test_that("shrinkage_report() quantifies movement toward the global effect", {
  skip_if_not_installed("posterior")
  d <- .sim_bayes_data(n = 400, seed = 3, tau = 0.7)
  fit <- bayes_surv_fit(
    d,
    iter = 1500, chains = 2, seed = 12, backend = "mcmc",
    cores = 2, prior_scale = 0.5
  )
  cx <- cox_fit(d, endpoint = "death")
  sh <- shrinkage_report(fit, cx)

  expect_s3_class(sh, "tbl_df")
  expect_equal(nrow(sh), 4L)
  expect_true(all(c(
    "standalone_log_hr", "pooled_log_hr", "global_log_hr",
    "shrinkage_abs", "pct_toward_global", "pct_width_reduction"
  ) %in% names(sh)))
  expect_true(any(is.finite(sh$pooled_log_hr)))
  expect_true(any(sh$pct_width_reduction > 0, na.rm = TRUE))

  # Borrowing should reduce the average distance to the global effect.
  dist_standalone <- mean(abs(sh$standalone_log_hr - sh$global_log_hr),
    na.rm = TRUE
  )
  dist_pooled <- mean(abs(sh$pooled_log_hr - sh$global_log_hr), na.rm = TRUE)
  expect_lte(dist_pooled, dist_standalone + 0.15)
})

test_that("compare_models() aligns the frequentist and Bayesian estimates", {
  skip_if_not_installed("posterior")
  d <- .sim_bayes_data(n = 260, seed = 4, tau = 0.4)
  fit <- bayes_surv_fit(
    d,
    iter = 800, chains = 2, seed = 13, backend = "mcmc",
    cores = 2, prior_scale = 0.5
  )
  cx <- cox_fit(d, endpoint = "death")
  cmp <- compare_models(cx, fit, d)

  expect_s3_class(cmp, "tbl_df")
  expect_equal(nrow(cmp), 2L)
  expect_setequal(cmp$model, c("frequentist", "bayesian"))
  expect_true(all(c("estimate", "conf.low", "conf.high", "c_index") %in% names(cmp)))
  expect_true(all(cmp$c_index > 0.5 & cmp$c_index < 1))
})

test_that("the brms backend fits the same hierarchical model", {
  skip_on_ci()
  skip_if_not_installed("brms")
  skip_if_not_installed("rstan")
  skip_if(
    Sys.getenv("TRIALSURV_TEST_BRMS") != "true",
    "Stan smoke test is opt-in; set TRIALSURV_TEST_BRMS=true"
  )
  d <- .sim_bayes_data(n = 180, seed = 5, tau = 0.4)
  fit <- bayes_surv_fit(
    d,
    iter = 1000, chains = 2, seed = 14, backend = "brms",
    adapt_delta = 0.95
  )
  expect_s3_class(fit, "trialsurv_bayes")
  expect_identical(fit$backend, "brms")
  bs <- bayes_summary(fit)
  expect_true(any(bs$type == "hazard_ratio"))
  expect_true(all(is.finite(bs$rhat[bs$parameter %in% c("b_lev", "b_f5", "shape")])))

  cmp <- compare_models(cox_fit(d, endpoint = "death"), fit, d)
  expect_equal(nrow(cmp), 2L)
})

test_that("bayes_surv_fit() validates arguments and supports pooling all arms", {
  skip_if_not_installed("posterior")
  d <- .sim_bayes_data(n = 150, seed = 7, tau = 0.4)

  expect_error(bayes_surv_fit(d, prior_scale = 0),
    class = "trialsurv_input_error"
  )
  expect_error(bayes_surv_fit(d, iter = 50), class = "trialsurv_input_error")
  expect_error(bayes_surv_fit(d, chains = 0), class = "trialsurv_input_error")
  expect_error(
    bayes_surv_fit(d, formula = time ~ 1, backend = "mcmc"),
    class = "trialsurv_input_error"
  )
  expect_error(bayes_surv_fit(d[0, ]), class = "trialsurv_input_error")
  expect_error(bayes_surv_fit(d, stratum = "not_a_column"),
    class = "trialsurv_input_error"
  )
  expect_error(bayes_summary(list()), class = "trialsurv_input_error")

  fake <- structure(list(hr_draws = NULL), class = c("trialsurv_bayes", "list"))
  expect_error(shrinkage_report(fake, NULL), class = "trialsurv_input_error")
  expect_error(plot_posterior_hr(fake), class = "trialsurv_input_error")
  expect_error(compare_models(1, fake, d), class = "trialsurv_input_error")

  fit_all <- bayes_surv_fit(
    d,
    iter = 600, chains = 2, seed = 8, backend = "mcmc",
    cores = 2, prior_scale = 0.5, pool = "all"
  )
  expect_true(all(fit_all$prep$arms$pooled))
  bs_all <- bayes_summary(fit_all)
  expect_true(any(grepl("^hr_lev_", bs_all$parameter)))
  expect_true(any(grepl("^hr_f5_", bs_all$parameter)))
})

test_that("compare_models() accepts a plain coxph fit", {
  skip_if_not_installed("posterior")
  d <- .sim_bayes_data(n = 200, seed = 9, tau = 0.4)
  fit <- bayes_surv_fit(
    d,
    iter = 600, chains = 2, seed = 10, backend = "mcmc",
    cores = 2, prior_scale = 0.5
  )
  cx <- cox_fit(d, endpoint = "death")
  cmp <- compare_models(cx$adjusted, fit, d)
  expect_equal(nrow(cmp), 2L)
  expect_setequal(cmp$model, c("frequentist", "bayesian"))
})
