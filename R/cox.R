#' Cox proportional hazards models
#'
#' Fits an unadjusted treatment-only Cox model and a covariate-adjusted Cox
#' model with Efron handling of ties, and returns both fits together with a
#' tidy table of hazard ratios, confidence intervals and p-values.
#'
#' @param data A tidy trial data frame, usually from [load_trial_data()].
#' @param covariates Character vector of adjustment covariates. The treatment
#'   variable given by `by` is always included.
#' @param by Character scalar naming the treatment variable (default `"rx"`).
#' @param endpoint One of `"death"` or `"recurrence"`.
#' @param conf.int Confidence level for the hazard ratio intervals.
#'
#' @return A list of class `trialsurv_cox` with elements `unadjusted` and
#'   `adjusted` (both [survival::coxph()] fits), `hr` (a tidy tibble with
#'   columns `term`, `label`, `model`, `hr`, `conf.low`, `conf.high`,
#'   `p.value`, `n`, `n_events`), and the model metadata.
#'
#' @examples
#' d <- load_trial_data()
#' cx <- cox_fit(d, endpoint = "death")
#' cx$hr
#'
#' @export
cox_fit <- function(data,
                    covariates = c("age", "sex", "nodes", "extent", "differ", "obstruct"),
                    by = "rx", endpoint = "death", conf.int = 0.95) {
  endpoint <- match.arg(endpoint, c("recurrence", "death"))
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.", class = "trialsurv_input_error")
  }
  if (!is.character(by) || length(by) != 1L) {
    rlang::abort("`by` must be a single column name.", class = "trialsurv_input_error")
  }
  if (!is.character(covariates)) {
    rlang::abort("`covariates` must be a character vector.",
      class = "trialsurv_input_error"
    )
  }

  d <- data[data$endpoint == endpoint, , drop = FALSE]
  if (nrow(d) == 0L) {
    rlang::abort(sprintf("No rows for endpoint %s.", encodeString(endpoint, quote = "\"")),
      class = "trialsurv_input_error"
    )
  }
  assert_surv_inputs(d, vars = unique(c(by, covariates)))

  for (v in unique(c(by, covariates))) {
    if (is.character(d[[v]])) d[[v]] <- factor(d[[v]])
    if (is.factor(d[[v]])) d[[v]] <- droplevels(d[[v]])
  }

  f_un <- stats::reformulate(by, response = "survival::Surv(time, status)")
  f_adj <- stats::reformulate(unique(c(by, covariates)),
    response = "survival::Surv(time, status)"
  )

  fit_un <- survival::coxph(f_un,
    data = d, ties = "efron",
    x = TRUE, y = TRUE, model = TRUE
  )
  fit_adj <- survival::coxph(f_adj,
    data = d, ties = "efron",
    x = TRUE, y = TRUE, model = TRUE
  )

  hr <- dplyr::bind_rows(
    .tidy_cox(fit_un, "unadjusted", conf.int),
    .tidy_cox(fit_adj, "adjusted", conf.int)
  )

  out <- list(
    unadjusted = fit_un,
    adjusted = fit_adj,
    hr = hr,
    n = nrow(d),
    n_events = sum(d$status == 1L),
    endpoint = endpoint,
    covariates = covariates,
    by = by,
    formula_unadjusted = f_un,
    formula_adjusted = f_adj
  )
  class(out) <- c("trialsurv_cox", "list")
  out
}

#' Tidy a coxph fit into a hazard-ratio table
#' @noRd
.tidy_cox <- function(fit, model, conf.int = 0.95) {
  s <- summary(fit, conf.int = conf.int)
  ci <- s$conf.int
  cf <- s$coefficients
  term <- rownames(cf)
  tibble::tibble(
    term = term,
    label = .pretty_term(term),
    model = model,
    hr = as.numeric(ci[, 1L]),
    conf.low = as.numeric(ci[, 3L]),
    conf.high = as.numeric(ci[, 4L]),
    p.value = as.numeric(cf[, ncol(cf)]),
    n = as.integer(fit$n),
    n_events = as.integer(fit$nevent)
  )
}

#' Readable labels for common model terms
#' @noRd
.pretty_term <- function(term) {
  out <- term
  out <- sub("^rxLev\\+5FU$", "Lev+5FU vs Obs", out)
  out <- sub("^rxLev$", "Lev vs Obs", out)
  out <- sub("^sexMale$", "Male vs Female", out)
  out <- sub("^extentMuscle$", "Muscle vs Submucosa", out)
  out <- sub("^extentSerosa$", "Serosa vs Submucosa", out)
  out <- sub("^extentContiguous$", "Contiguous vs Submucosa", out)
  out <- sub("^differModerate$", "Moderate vs Well", out)
  out <- sub("^differPoor$", "Poor vs Well", out)
  out <- sub("^surg$", "Time since surgery", out)
  out <- sub("^age$", "Age (per year)", out)
  out <- sub("^nodes$", "Nodes (per node)", out)
  out <- sub("^obstruct$", "Obstruction (yes vs no)", out)
  out <- sub("^perfor$", "Perforation (yes vs no)", out)
  out <- sub("^adhere$", "Adherence (yes vs no)", out)
  out
}

#' Proportional hazards diagnostics for a Cox model
#'
#' Runs [survival::cox.zph()] with a Kaplan-Meier time transform, extracts the
#' per-coefficient Schoenfeld residual slopes, counts deviance-residual
#' outliers, and flags covariates that violate proportional hazards at
#' `alpha`.
#'
#' @param fit A [survival::coxph()] fit or a `trialsurv_cox` object (the
#'   adjusted model is used).
#' @param alpha Significance level for the proportional hazards flag.
#' @param outlier_threshold Absolute deviance residual above which an
#'   observation is counted as an outlier.
#'
#' @return A tibble with one row per model coefficient plus a `"GLOBAL"` row,
#'   containing the cox.zph statistic, df and p-value, the Schoenfeld residual
#'   slope, and cohort-level deviance outlier counts. The names of covariates
#'   violating proportional hazards are attached as `attr(x, "violators")`.
#'
#' @examples
#' d <- load_trial_data()
#' cx <- cox_fit(d)
#' cox_diagnostics(cx)
#'
#' @export
cox_diagnostics <- function(fit, alpha = 0.05, outlier_threshold = 3) {
  if (inherits(fit, "trialsurv_cox")) fit <- fit$adjusted
  if (!inherits(fit, "coxph")) {
    rlang::abort("`fit` must be a coxph or trialsurv_cox object.",
      class = "trialsurv_input_error"
    )
  }
  if (!is.numeric(alpha) || length(alpha) != 1L || alpha <= 0 || alpha >= 1) {
    rlang::abort("`alpha` must be a single number in (0, 1).",
      class = "trialsurv_input_error"
    )
  }

  zph <- survival::cox.zph(fit, transform = "km", terms = FALSE)
  tab <- as.data.frame(zph$table)
  term <- rownames(tab)

  res <- stats::residuals(fit, type = "schoenfeld")
  res <- as.matrix(res)
  tt <- attr(res, "time")
  if (is.null(tt)) tt <- suppressWarnings(as.numeric(rownames(res)))
  tt <- as.numeric(tt)
  if (length(tt) != nrow(res)) {
    rlang::abort("Could not align Schoenfeld residuals with their times.",
      class = "trialsurv_input_error"
    )
  }

  slope <- rep(NA_real_, length(term))
  slope_se <- rep(NA_real_, length(term))
  slope_p <- rep(NA_real_, length(term))
  for (j in seq_len(ncol(res))) {
    nm <- colnames(res)[j]
    k <- match(nm, term)
    if (is.na(k)) next
    ok <- is.finite(res[, j]) & is.finite(tt)
    if (sum(ok) < 3L) next
    m <- stats::lm(res[ok, j] ~ tt[ok])
    co <- summary(m)$coefficients
    slope[k] <- co[2L, 1L]
    slope_se[k] <- co[2L, 2L]
    slope_p[k] <- co[2L, 4L]
  }

  dv <- stats::residuals(fit, type = "deviance")
  n_out <- sum(abs(dv) > outlier_threshold, na.rm = TRUE)
  n_dev <- sum(is.finite(dv))

  out <- tibble::tibble(
    term = term,
    statistic = as.numeric(tab$chisq),
    df = as.numeric(tab$df),
    p_value = as.numeric(tab$p),
    ph_violation = ifelse(term == "GLOBAL", NA, as.numeric(tab$p) < alpha),
    slope = slope,
    slope_se = slope_se,
    slope_p = slope_p,
    n = as.integer(fit$n),
    n_events = as.integer(fit$nevent),
    deviance_outliers = as.integer(n_out),
    deviance_outlier_pct = 100 * n_out / n_dev
  )

  attr(out, "violators") <- out$term[which(out$ph_violation)]
  attr(out, "alpha") <- alpha
  out
}

#' Harrell's concordance for a Cox model on new data
#'
#' Computes Harrell's C-index on held-out data using the Cox linear predictor
#' as the risk score, with an approximate standard error and confidence
#' interval from [survival::concordance()].
#'
#' @param fit A [survival::coxph()] fit or a `trialsurv_cox` object (the
#'   adjusted model is used).
#' @param newdata A data frame with the same covariates and the `time` and
#'   `status` columns.
#' @param conf.int Confidence level for the C-index interval.
#'
#' @return A one-row tibble with `model`, `c_index`, `std_error`,
#'   `conf.low`, `conf.high`, `n` and `n_events`.
#'
#' @examples
#' d <- load_trial_data()
#' ids <- unique(d$id)
#' set.seed(1)
#' train_ids <- sample(ids, size = floor(0.7 * length(ids)))
#' train <- d[d$id %in% train_ids, ]
#' test <- d[!d$id %in% train_ids, ]
#' cx <- cox_fit(train)
#' concordance_eval(cx, test)
#'
#' @export
concordance_eval <- function(fit, newdata, conf.int = 0.95) {
  if (inherits(fit, "trialsurv_cox")) fit <- fit$adjusted
  if (!inherits(fit, "coxph")) {
    rlang::abort("`fit` must be a coxph or trialsurv_cox object.",
      class = "trialsurv_input_error"
    )
  }
  assert_surv_inputs(newdata)

  lp <- stats::predict(fit, newdata = newdata, type = "lp")
  d <- newdata
  d$.risk_score <- as.numeric(lp)
  d <- d[is.finite(d$.risk_score), , drop = FALSE]
  if (nrow(d) == 0L) {
    rlang::abort("No finite linear predictors on `newdata`.",
      class = "trialsurv_input_error"
    )
  }

  cc <- survival::concordance(
    survival::Surv(time, status) ~ .risk_score,
    data = d, reverse = TRUE
  )
  est <- as.numeric(cc$concordance)
  se <- if (!is.null(cc$std.err)) {
    as.numeric(cc$std.err)
  } else if (!is.null(cc$var)) {
    sqrt(as.numeric(cc$var))
  } else {
    NA_real_
  }
  z <- stats::qnorm(1 - (1 - conf.int) / 2)
  tibble::tibble(
    model = "cox",
    c_index = est,
    std_error = se,
    conf.low = max(0, est - z * se),
    conf.high = min(1, est + z * se),
    n = nrow(d),
    n_events = sum(d$status == 1L)
  )
}
