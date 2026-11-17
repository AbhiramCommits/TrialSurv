#' Kaplan-Meier survival estimation
#'
#' Fits Kaplan-Meier curves for one endpoint of the tidied trial data and
#' returns the fitted `survfit` object together with tidy summaries: median
#' survival with a confidence interval, landmark survival probabilities at
#' 1, 3 and 5 years, and the log-rank test comparing groups.
#'
#' @param data A tidy trial data frame, usually from [load_trial_data()].
#' @param by Character scalar naming the grouping column (default `"rx"`),
#'   or `NULL` for a single overall curve.
#' @param endpoint One of `"death"` or `"recurrence"`.
#' @param times Numeric vector of landmark times in months. Defaults to
#'   `c(12, 36, 60)` (1, 3 and 5 years).
#' @param conf.int Confidence level for the survival intervals.
#'
#' @return A list of class `trialsurv_km` with elements:
#'   * `fit`: the [survival::survfit()] object;
#'   * `summary`: a wide tibble with one row per group containing n, events,
#'     median survival with CI, landmark survival with CIs and the log-rank
#'     statistics;
#'   * `median`: a tidy tibble of median survival;
#'   * `landmark`: a long tibble of landmark survival probabilities;
#'   * `logrank`: a one-row tibble with the log-rank chi-square, df and p-value.
#'
#' @examples
#' d <- load_trial_data()
#' km <- km_fit(d, by = "rx", endpoint = "death")
#' km$summary
#' km$logrank
#'
#' @export
km_fit <- function(data, by = "rx", endpoint = "death",
                   times = c(12, 36, 60), conf.int = 0.95) {
  endpoint <- match.arg(endpoint, c("recurrence", "death"))
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.", class = "trialsurv_input_error")
  }
  if (!is.null(by) && (!is.character(by) || length(by) != 1L)) {
    rlang::abort("`by` must be NULL or a single column name.",
      class = "trialsurv_input_error"
    )
  }
  if (!all(c("time", "status", "endpoint") %in% names(data))) {
    rlang::abort(
      "`data` must contain time, status and endpoint columns.",
      class = "trialsurv_input_error"
    )
  }

  d <- data[data$endpoint == endpoint, , drop = FALSE]
  if (nrow(d) == 0L) {
    rlang::abort(sprintf("No rows for endpoint %s.", encodeString(endpoint, quote = "\"")),
      class = "trialsurv_input_error"
    )
  }
  assert_surv_inputs(d, vars = by)

  rhs <- if (is.null(by)) "1" else by
  form <- stats::reformulate(rhs, response = "survival::Surv(time, status)")
  fit <- survival::survfit(form, data = d, conf.int = conf.int)

  tab <- summary(fit)$table
  if (is.null(dim(tab))) {
    tab <- matrix(tab, nrow = 1L, dimnames = list("Overall", names(tab)))
  }
  arm <- sub("^.*=", "", rownames(tab))
  if (is.null(by)) arm <- "Overall"

  lcl <- grep("LCL$", colnames(tab))[1]
  ucl <- grep("UCL$", colnames(tab))[1]
  median_tbl <- tibble::tibble(
    arm = arm,
    n = as.integer(tab[, "records"]),
    n_events = as.integer(tab[, "events"]),
    median = as.numeric(tab[, "median"]),
    median_lower = as.numeric(tab[, lcl]),
    median_upper = as.numeric(tab[, ucl])
  )

  s <- summary(fit, times = times, extend = TRUE)
  s_arm <- if (is.null(s$strata)) {
    rep("Overall", length(s$time))
  } else {
    sub("^.*=", "", as.character(s$strata))
  }
  landmark <- tibble::tibble(
    arm = s_arm,
    time = as.numeric(s$time),
    surv = as.numeric(s$surv),
    lower = as.numeric(s$lower),
    upper = as.numeric(s$upper)
  )

  landmark_wide <- landmark |>
    dplyr::mutate(label = paste0(time, "m")) |>
    dplyr::select("arm", "label", "surv", "lower", "upper") |>
    tidyr::pivot_wider(
      names_from = "label",
      values_from = c("surv", "lower", "upper"),
      names_glue = "{.value}_{label}"
    )

  summary_tbl <- dplyr::left_join(median_tbl, landmark_wide, by = "arm")

  if (is.null(by)) {
    logrank <- tibble::tibble(
      test = "log-rank", chisq = NA_real_, df = NA_integer_, p_value = NA_real_
    )
  } else {
    sd <- survival::survdiff(form, data = d)
    df <- length(sd$n) - 1L
    logrank <- tibble::tibble(
      test = "log-rank",
      chisq = unname(sd$chisq),
      df = as.integer(df),
      p_value = stats::pchisq(sd$chisq, df = df, lower.tail = FALSE)
    )
  }

  summary_tbl <- dplyr::mutate(
    summary_tbl,
    logrank_chisq = logrank$chisq,
    logrank_df = logrank$df,
    logrank_p = logrank$p_value
  )

  out <- list(
    fit = fit,
    summary = summary_tbl,
    median = median_tbl,
    landmark = landmark,
    logrank = logrank,
    endpoint = endpoint,
    by = by
  )
  class(out) <- c("trialsurv_km", "list")
  out
}
