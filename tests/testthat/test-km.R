test_that("km_fit() returns sane medians and stable arm ordering", {
  d <- load_trial_data()
  km <- km_fit(d, by = "rx", endpoint = "death")

  expect_s3_class(km, "trialsurv_km")
  expect_s3_class(km$fit, "survfit")
  expect_equal(nrow(km$median), 3L)
  expect_true(all(c("arm", "median", "median_lower", "median_upper") %in% names(km$median)))

  # Medians that are reached must lie inside the observed follow-up range.
  # NA medians mean survival never fell to 50% during follow-up.
  obs_range <- range(d$time[d$endpoint == "death"])
  reached <- !is.na(km$median$median)
  expect_gt(sum(reached), 0L)
  expect_true(all(km$median$median[reached] >= obs_range[1]))
  expect_true(all(km$median$median[reached] <= obs_range[2]))

  # Landmark survival is non-increasing in time within each arm.
  wide <- km$summary
  expect_true(all(wide$surv_12m >= wide$surv_36m))
  expect_true(all(wide$surv_36m >= wide$surv_60m))
  expect_true(all(wide$surv_60m >= 0 & wide$surv_60m <= 1))

  # The 5-FU arm should not do worse than observation on the median scale,
  # and the ordering is deterministic across repeated calls.
  med <- stats::setNames(km$median$median, km$median$arm)
  finite_med <- med[is.finite(med)]
  expect_true("Obs" %in% names(finite_med))
  # If the 5-FU median is reached it should not be shorter than observation.
  if ("Lev+5FU" %in% names(finite_med)) {
    expect_lte(finite_med[["Obs"]], finite_med[["Lev+5FU"]] + 1e-8)
  }

  km2 <- km_fit(d, by = "rx", endpoint = "death")
  expect_equal(km$median, km2$median)
  expect_equal(km$logrank$p_value, km2$logrank$p_value)
})

test_that("km_fit() computes a log-rank test and supports the overall curve", {
  d <- load_trial_data()
  km <- km_fit(d, by = "rx", endpoint = "recurrence")
  expect_true(is.finite(km$logrank$chisq))
  expect_true(km$logrank$p_value >= 0 && km$logrank$p_value <= 1)

  overall <- km_fit(d, by = NULL, endpoint = "death")
  expect_equal(nrow(overall$median), 1L)
  expect_equal(overall$median$arm, "Overall")
  expect_true(is.na(overall$logrank$p_value))
})

test_that("plot_km() returns a plot object", {
  skip_if_not_installed("ggsurvfit")
  d <- load_trial_data()
  p <- plot_km(d, by = "rx", endpoint = "death", save = FALSE)
  expect_true(inherits(p, "ggplot") || inherits(p, "patchwork"))
})

test_that("plot_km() saves a figure when requested", {
  skip_if_not_installed("ggsurvfit")
  d <- load_trial_data()
  tmp <- tempfile("trialsurv-km-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)
  plot_km(d,
    by = "rx", endpoint = "death", save = TRUE,
    dir = tmp, name = "km.png"
  )
  expect_true(file.exists(file.path(tmp, "km.png")))
})
