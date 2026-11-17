test_that("trial_split() is deterministic given a seed", {
  d <- load_trial_data()
  a <- trial_split(d, prop = 0.7, seed = 42)
  b <- trial_split(d, prop = 0.7, seed = 42)
  c <- trial_split(d, prop = 0.7, seed = 43)

  expect_s3_class(a, "trialsurv_split")
  expect_identical(a$train_ids, b$train_ids)
  expect_false(identical(a$train_ids, c$train_ids))
  expect_identical(a$train, b$train)
})

test_that("train and test are disjoint and cover every patient", {
  d <- load_trial_data()
  sp <- trial_split(d, prop = 0.7, seed = 20261117)

  expect_length(intersect(sp$train_ids, sp$test_ids), 0L)
  expect_setequal(c(sp$train_ids, sp$test_ids), unique(d$id))
  expect_setequal(unique(sp$train$id), sp$train_ids)
  expect_setequal(unique(sp$test$id), sp$test_ids)

  # Both endpoints travel with each patient.
  expect_equal(
    sum(sp$train$endpoint == "death"),
    length(sp$train_ids)
  )
  expect_equal(
    sum(sp$test$endpoint == "recurrence"),
    length(sp$test_ids)
  )
})

test_that("the split is stratified by arm and event status", {
  d <- load_trial_data()
  sp <- trial_split(d, prop = 0.7, seed = 20261117)
  tr <- sp$train[sp$train$endpoint == "death", , drop = FALSE]
  te <- sp$test[sp$test$endpoint == "death", , drop = FALSE]

  for (arm in levels(d$rx)) {
    for (ev in c(0, 1)) {
      n_tr <- sum(tr$rx == arm & tr$status == ev)
      n_te <- sum(te$rx == arm & te$status == ev)
      if (n_tr + n_te >= 20) {
        expect_equal(n_tr / (n_tr + n_te), 0.7, tolerance = 0.12)
      }
    }
  }
  expect_equal(sum(sp$summary$n_patients), length(unique(d$id)))
})
