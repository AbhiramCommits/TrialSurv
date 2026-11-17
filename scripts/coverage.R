#!/usr/bin/env Rscript
# Measure testthat coverage and record the test-count summary. This is a
# development helper, not part of the analysis pipeline. The optional Stan
# smoke test is enabled when brms/rstan are available; it is still skipped on
# CI by testthat::skip_on_ci().

Sys.setenv(TRIALSURV_TEST_BRMS = "true")
suppressPackageStartupMessages(library(testthat))

res <- test_local(reporter = "silent")
df <- as.data.frame(res)
test_summary <- data.frame(
  test_files = length(unique(df$file)),
  test_blocks = nrow(df),
  expectations = sum(df$nb),
  failures = sum(df$failed),
  skipped = sum(df$skipped)
)
dir.create("results", showWarnings = FALSE)
utils::write.csv(test_summary, "results/test_summary.csv", row.names = FALSE)
cat("Test summary:\n")
print(test_summary)

cov <- covr::package_coverage()
pct <- covr::percent_coverage(cov)
writeLines(
  c(
    sprintf("Overall coverage: %.2f%%", pct),
    "",
    capture.output(print(cov))
  ),
  "results/coverage.txt"
)
cat(sprintf("Coverage: %.2f%%\n", pct))
