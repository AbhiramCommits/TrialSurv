#' Compare frequentist and Bayesian survival models
#'
#' Puts the primary treatment effect, interval width and held-out concordance
#' of a frequentist Cox model and a Bayesian hierarchical model side by side.
#' When the Bayesian fit used the brms backend, WAIC is added; the fallback
#' sampler has no WAIC because it does not produce a pointwise log likelihood.
#'
#' @param freq_fit A `trialsurv_cox` object (or a [survival::coxph()] fit).
#' @param bayes_fit A `trialsurv_bayes` fit.
#' @param test_data Held-out data used for both concordance calculations.
#'
#' @return A tibble with one row per model containing the treatment hazard
#'   ratio, its interval, interval width, held-out C-index and (for brms)
#'   WAIC.
#'
#' @examples
#' \dontrun{
#' d <- load_trial_data()
#' sp <- trial_split(d)
#' cx <- cox_fit(sp$train)
#' bf <- bayes_surv_fit(sp$train, iter = 1000, chains = 2)
#' compare_models(cx, bf, sp$test)
#' }
#'
#' @export
compare_models <- function(freq_fit, bayes_fit, test_data) {
  if (inherits(freq_fit, "trialsurv_cox")) {
    hr_tbl <- freq_fit$hr
  } else if (inherits(freq_fit, "coxph")) {
    hr_tbl <- .tidy_cox(freq_fit, "adjusted")
  } else {
    rlang::abort("`freq_fit` must be a trialsurv_cox or coxph object.",
      class = "trialsurv_input_error"
    )
  }
  if (!inherits(bayes_fit, "trialsurv_bayes")) {
    rlang::abort("`bayes_fit` must be a trialsurv_bayes object.",
      class = "trialsurv_input_error"
    )
  }

  freq_row <- hr_tbl[hr_tbl$model == "adjusted" & grepl("^rx", hr_tbl$term), ,
    drop = FALSE
  ]
  if (nrow(freq_row) == 0L) {
    freq_row <- hr_tbl[hr_tbl$model == "adjusted", , drop = FALSE][1, , drop = FALSE]
  }
  primary <- freq_row[grepl("5FU", freq_row$term), , drop = FALSE]
  if (nrow(primary) == 0L) primary <- freq_row[1, , drop = FALSE]

  bs <- bayes_summary(bayes_fit)
  hr_rows <- bs[bs$type == "hazard_ratio" &
    grepl("_global$", bs$parameter), , drop = FALSE]
  bprimary <- hr_rows[grepl("_f5_global$", hr_rows$parameter), , drop = FALSE]
  if (nrow(bprimary) == 0L && nrow(hr_rows) > 0L) {
    bprimary <- hr_rows[nrow(hr_rows), , drop = FALSE]
  }

  cc_f <- concordance_eval(freq_fit, test_data)
  cc_b <- bayes_concordance(bayes_fit, test_data)

  waic <- elpd_waic <- p_waic <- NA_real_
  if (identical(bayes_fit$backend, "brms") && !isTRUE(bayes_fit$custom_formula)) {
    w <- tryCatch(brms::waic(bayes_fit$fit), error = function(e) NULL)
    if (!is.null(w)) {
      waic <- as.numeric(w$estimates["waic", "Estimate"])
      elpd_waic <- as.numeric(w$estimates["elpd_waic", "Estimate"])
      p_waic <- as.numeric(w$estimates["p_waic", "Estimate"])
    }
  }

  tibble::tibble(
    model = c("frequentist", "bayesian"),
    parameter = c(primary$label[1], bprimary$parameter[1]),
    estimate = c(primary$hr[1], bprimary$mean[1]),
    conf.low = c(primary$conf.low[1], bprimary$q2.5[1]),
    conf.high = c(primary$conf.high[1], bprimary$q97.5[1]),
    interval_width = c(
      primary$conf.high[1] - primary$conf.low[1],
      bprimary$q97.5[1] - bprimary$q2.5[1]
    ),
    c_index = c(cc_f$c_index[1], cc_b$c_index[1]),
    c_index_low = c(cc_f$conf.low[1], cc_b$conf.low[1]),
    c_index_high = c(cc_f$conf.high[1], cc_b$conf.high[1]),
    waic = c(NA_real_, waic),
    elpd_waic = c(NA_real_, elpd_waic),
    p_waic = c(NA_real_, p_waic),
    n_test = c(cc_f$n[1], cc_b$n[1]),
    n_events_test = c(cc_f$n_events[1], cc_b$n_events[1])
  )
}
