test_that(".check_report_inputs() lists every missing artifact", {
  tmp <- tempfile("trialsurv-report-")
  dir.create(file.path(tmp, "results"), recursive = TRUE)
  dir.create(file.path(tmp, "figures"), recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  err <- tryCatch(
    .check_report_inputs(file.path(tmp, "results"), file.path(tmp, "figures")),
    error = identity
  )
  expect_s3_class(err, "trialsurv_report_error")
  expect_match(conditionMessage(err), "cohort_summary.csv")

  results <- c(
    "cohort_summary.csv", "km_summary.csv", "cox_hazard_ratios.csv",
    "ph_diagnostics.csv", "concordance_frequentist.csv",
    "bayes_posterior_summary.csv", "shrinkage_report.csv",
    "model_comparison.csv", "bayes_convergence.csv", "power_curve.csv",
    "type1_error.csv", "sample_size_vs_schoenfeld.csv", "simulation_meta.csv"
  )
  figures <- c(
    "km_by_arm.png", "schoenfeld.png", "posterior_hr.png", "power_curves.png"
  )
  file.create(file.path(tmp, "results", results))
  file.create(file.path(tmp, "figures", figures))
  expect_true(.check_report_inputs(file.path(tmp, "results"), file.path(tmp, "figures")))

  unlink(file.path(tmp, "results", "km_summary.csv"))
  err2 <- tryCatch(
    .check_report_inputs(file.path(tmp, "results"), file.path(tmp, "figures")),
    error = identity
  )
  expect_s3_class(err2, "trialsurv_report_error")
  expect_match(conditionMessage(err2), "km_summary.csv")
})

test_that("render_report() errors clearly when the template is absent", {
  skip_if_not_installed("quarto")
  tmp <- tempfile("trialsurv-report-template-")
  dir.create(tmp)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  err <- tryCatch(render_report(output_dir = tmp), error = identity)
  expect_s3_class(err, "trialsurv_report_error")
  expect_match(conditionMessage(err), "Report template not found")
})

test_that("render_report() renders a minimal report when inputs exist", {
  skip_if_not_installed("quarto")
  tmp <- tempfile("trialsurv-report-render-")
  dir.create(file.path(tmp, "results"), recursive = TRUE)
  dir.create(file.path(tmp, "figures"), recursive = TRUE)
  on.exit(unlink(tmp, recursive = TRUE), add = TRUE)

  results <- c(
    "cohort_summary.csv", "km_summary.csv", "cox_hazard_ratios.csv",
    "ph_diagnostics.csv", "concordance_frequentist.csv",
    "bayes_posterior_summary.csv", "shrinkage_report.csv",
    "model_comparison.csv", "bayes_convergence.csv", "power_curve.csv",
    "type1_error.csv", "sample_size_vs_schoenfeld.csv", "simulation_meta.csv"
  )
  figures <- c(
    "km_by_arm.png", "schoenfeld.png", "posterior_hr.png", "power_curves.png"
  )
  file.create(file.path(tmp, "results", results))
  file.create(file.path(tmp, "figures", figures))
  writeLines(
    c("---", "title: smoke", "format: html", "---", "", "Minimal report."),
    file.path(tmp, "trial_report.qmd")
  )

  html <- render_report(
    output_dir = tmp,
    results_dir = file.path(tmp, "results"),
    figures_dir = file.path(tmp, "figures"),
    quiet = TRUE
  )
  expect_true(file.exists(html))
  expect_match(html, "trial_report\\.html$")
})
