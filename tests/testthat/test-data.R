test_that("load_trial_data() returns the documented tidy columns", {
  d <- load_trial_data()

  expected <- c(
    "id", "time", "status", "endpoint", "rx", "sex", "age", "obstruct",
    "perfor", "adhere", "nodes", "differ", "extent", "surg", "node4"
  )
  expect_true(all(expected %in% names(d)))
  expect_s3_class(d, "tbl_df")
  expect_gt(nrow(d), 0)
})

test_that("status is 0/1 and both endpoints are present", {
  d <- load_trial_data()

  expect_type(d$status, "integer")
  expect_true(all(d$status %in% c(0L, 1L)))
  expect_setequal(unique(as.character(d$endpoint)), c("recurrence", "death"))
  expect_setequal(levels(d$rx), c("Obs", "Lev", "Lev+5FU"))

  # One row per patient per endpoint, and the two endpoints have equal n.
  counts <- table(d$id)
  expect_true(all(counts == 2L))
  expect_equal(
    sum(d$endpoint == "recurrence"),
    sum(d$endpoint == "death")
  )
})

test_that("analysis variables are complete after the default drop", {
  d <- load_trial_data()
  vars <- c(
    "time", "status", "rx", "sex", "age", "obstruct", "perfor",
    "adhere", "nodes", "differ", "extent", "surg", "node4"
  )
  expect_false(anyNA(d[, vars]))
  expect_true(all(d$time > 0))
})

test_that("the cohort attribute records n, events and follow-up", {
  d <- load_trial_data()
  cohort <- attr(d, "cohort")

  expect_s3_class(cohort, "tbl_df")
  expect_named(cohort, c("endpoint", "n", "n_events", "median_followup"))
  expect_equal(sort(cohort$endpoint), c("death", "recurrence"))
  expect_true(all(cohort$n > 0))
  expect_true(all(cohort$n_events > 0))
  expect_true(all(cohort$n_events <= cohort$n))
  expect_true(all(cohort$median_followup > 0))

  miss <- attr(d, "missing")
  expect_s3_class(miss, "tbl_df")
  expect_true(all(c("variable", "n_missing") %in% names(miss)))
})

test_that("na_action = 'keep' retains flagged missing values", {
  kept <- load_trial_data(na_action = "keep")
  dropped <- load_trial_data()

  expect_gt(nrow(kept), nrow(dropped))
  expect_true(anyNA(kept$differ))
  expect_equal(attr(attr(kept, "missing"), "na_action"), "keep")
  expect_equal(attr(attr(dropped, "missing"), "na_action"), "drop")
  expect_error(load_trial_data("nope"), class = "trialsurv_input_error")
})
