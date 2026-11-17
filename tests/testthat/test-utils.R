test_that("assert_surv_inputs() validates the survival data contract", {
  ok <- data.frame(time = c(1, 2, 3), status = c(1, 0, 1))
  expect_true(assert_surv_inputs(ok))

  expect_error(assert_surv_inputs(1:3), class = "trialsurv_input_error")
  expect_error(
    assert_surv_inputs(data.frame(time = 1, status = 1), vars = "age"),
    class = "trialsurv_input_error"
  )
  expect_error(
    assert_surv_inputs(data.frame(time = "a", status = 1)),
    class = "trialsurv_input_error"
  )
  expect_error(
    assert_surv_inputs(data.frame(time = 0, status = 1)),
    class = "trialsurv_input_error"
  )
  expect_error(
    assert_surv_inputs(data.frame(time = c(1, NA), status = c(1, 0))),
    class = "trialsurv_input_error"
  )
  expect_error(
    assert_surv_inputs(data.frame(time = 1, status = 2)),
    class = "trialsurv_input_error"
  )
})

test_that("write_result() writes CSV and JSON and returns its input", {
  tmp <- tempfile("trialsurv-results-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  df <- data.frame(a = 1:3, b = c("x", "y", "z"))
  expect_invisible(write_result(df, "table", dir = tmp))
  csv <- file.path(tmp, "table.csv")
  expect_true(file.exists(csv))
  expect_equal(read.csv(csv)$a, 1:3)

  lst <- list(alpha = 1, beta = "two", gamma = c(3, 4))
  expect_invisible(write_result(lst, "meta", dir = tmp))
  json <- file.path(tmp, "meta.json")
  expect_true(file.exists(json))
  parsed <- jsonlite::read_json(json, simplifyVector = TRUE)
  expect_equal(parsed$alpha, 1)
  expect_equal(parsed$beta, "two")

  write_result(df, "table.json", dir = tmp)
  expect_true(file.exists(file.path(tmp, "table.json")))
  expect_error(
    write_result(list(a = 1:2, b = 3:4), "bad", dir = tmp, format = "csv"),
    class = "trialsurv_input_error"
  )
})

test_that("fig_save() writes a non-empty PNG", {
  skip_if_not_installed("ggplot2")
  tmp <- tempfile("trialsurv-figures-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  p <- ggplot2::ggplot(
    data.frame(x = 1:3, y = 3:1),
    ggplot2::aes(x = .data$x, y = .data$y)
  ) +
    ggplot2::geom_point()
  path <- fig_save(p, "smoke", dir = tmp, width = 3, height = 2, dpi = 72)

  expect_true(file.exists(path))
  expect_gt(file.size(path), 0)
  expect_match(path, "smoke\\.png$")
  expect_invisible(fig_save(p, "smoke2.png", dir = tmp, width = 3, height = 2, dpi = 72))
})
