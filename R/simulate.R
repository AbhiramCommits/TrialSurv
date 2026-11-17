#' Simulate a two-arm time-to-event trial
#'
#' Simulates a randomised two-arm trial with uniform staggered accrual,
#' exponential survival in the control arm, proportional hazards in the
#' experimental arm, administrative censoring at the analysis date and
#' independent random dropout. Event times are in months.
#'
#' @param n_per_arm Number of patients randomised to each arm.
#' @param hr Hazard ratio of the experimental arm versus control.
#' @param accrual_months Length of the uniform accrual window.
#' @param followup_months Additional follow-up after the last patient is
#'   randomised. The analysis occurs at `accrual_months + followup_months`.
#' @param dropout_rate Proportion of patients expected to drop out within
#'   `followup_months`; converted to an exponential dropout hazard. Set to 0
#'   for no dropout.
#' @param median_control Median survival in the control arm, in months.
#' @param alpha Significance level stored with the simulated dataset (it does
#'   not change the simulation).
#' @param seed Optional integer seed.
#'
#' @return A tibble with one row per patient and columns `id`, `arm`, `time`
#'   (observed follow-up), `status` (1 event, 0 censored), `event_time`,
#'   `dropout_time`, `admin_time` and `accrual_time`.
#'
#' @examples
#' d <- simulate_trial(150, hr = 0.7, accrual_months = 12,
#'   followup_months = 24, dropout_rate = 0.1, median_control = 24, seed = 1)
#' mean(d$status)
#'
#' @export
simulate_trial <- function(n_per_arm, hr, accrual_months, followup_months,
                           dropout_rate, median_control, alpha = 0.05,
                           seed = NULL) {
  if (!is.numeric(n_per_arm) || length(n_per_arm) != 1L || n_per_arm < 1) {
    rlang::abort("`n_per_arm` must be a positive integer.",
      class = "trialsurv_input_error"
    )
  }
  if (!is.numeric(hr) || length(hr) != 1L || hr <= 0) {
    rlang::abort("`hr` must be a positive number.", class = "trialsurv_input_error")
  }
  if (!is.numeric(accrual_months) || accrual_months <= 0) {
    rlang::abort("`accrual_months` must be positive.",
      class = "trialsurv_input_error"
    )
  }
  if (!is.numeric(followup_months) || followup_months <= 0) {
    rlang::abort("`followup_months` must be positive.",
      class = "trialsurv_input_error"
    )
  }
  if (!is.numeric(dropout_rate) || dropout_rate < 0 || dropout_rate >= 1) {
    rlang::abort("`dropout_rate` must be in [0, 1).",
      class = "trialsurv_input_error"
    )
  }
  if (!is.numeric(median_control) || median_control <= 0) {
    rlang::abort("`median_control` must be positive.",
      class = "trialsurv_input_error"
    )
  }
  if (!is.null(seed)) set.seed(as.integer(seed))

  n_per_arm <- as.integer(n_per_arm)
  n_total <- 2L * n_per_arm
  arm <- factor(rep(c("Control", "Treatment"), each = n_per_arm),
    levels = c("Control", "Treatment")
  )
  accrual_time <- stats::runif(n_total, 0, accrual_months)
  analysis_time <- accrual_months + followup_months
  admin_time <- analysis_time - accrual_time

  rate_control <- log(2) / median_control
  rate <- ifelse(arm == "Treatment", rate_control * hr, rate_control)
  event_time <- stats::rexp(n_total, rate = rate)

  if (dropout_rate > 0) {
    dropout_hazard <- -log(1 - dropout_rate) / followup_months
    dropout_time <- stats::rexp(n_total, rate = dropout_hazard)
  } else {
    dropout_time <- rep(Inf, n_total)
  }

  censor_time <- pmin(dropout_time, admin_time)
  time <- pmin(event_time, censor_time)
  status <- as.integer(event_time <= censor_time)

  out <- tibble::tibble(
    id = seq_len(n_total),
    arm = arm,
    time = time,
    status = status,
    event_time = event_time,
    dropout_time = dropout_time,
    admin_time = admin_time,
    accrual_time = accrual_time
  )
  attr(out, "hr") <- hr
  attr(out, "alpha") <- alpha
  attr(out, "median_control") <- median_control
  out
}

#' Monte Carlo power and type I error for a two-arm trial
#'
#' Runs [simulate_trial()] across a grid of design scenarios and returns the
#' empirical rejection rate, Monte Carlo standard error and interval, mean
#' observed events and mean estimated hazard ratio for each cell. At
#' `hr = 1` the rejection rate is the empirical type I error rate.
#'
#' @param grid A data frame with one row per scenario. Required columns are
#'   `n_per_arm`, `hr` and `dropout_rate`; `accrual_months`,
#'   `followup_months` and `median_control` default to 24, 60 and 60 months.
#' @param n_sims Number of simulated trials per cell.
#' @param test `"logrank"` or `"cox"` for the primary p-value. A Cox model is
#'   always fitted to report the estimated hazard ratio.
#' @param cores Number of cores for `parallel::mclapply()`; defaults to one
#'   fewer than the detected cores.
#' @param seed Integer seed used to derive a deterministic seed per
#'   simulation.
#' @param alpha Significance level.
#'
#' @return A tidy tibble with one row per scenario. The total runtime in
#'   seconds is attached as `attr(x, "runtime")`.
#'
#' @examples
#' grid <- data.frame(n_per_arm = c(100, 200), hr = c(1, 0.6),
#'   dropout_rate = 0)
#' trial_power(grid, n_sims = 100, test = "logrank", cores = 1, seed = 1)
#'
#' @export
trial_power <- function(grid, n_sims = 1000,
                        test = c("logrank", "cox"), cores = NULL,
                        seed = 20261117, alpha = 0.05) {
  test <- match.arg(test)
  if (!is.data.frame(grid)) {
    rlang::abort("`grid` must be a data frame.", class = "trialsurv_input_error")
  }
  if (!is.numeric(n_sims) || n_sims < 2) {
    rlang::abort("`n_sims` must be at least 2.", class = "trialsurv_input_error")
  }
  defaults <- list(
    n_per_arm = 300, hr = 0.7, accrual_months = 24,
    followup_months = 60, dropout_rate = 0, median_control = 60
  )
  for (nm in names(defaults)) {
    if (!nm %in% names(grid)) grid[[nm]] <- defaults[[nm]]
  }
  if (is.null(cores)) {
    cores <- max(1L, min(8L, parallel::detectCores() - 1L))
  }

  t0 <- Sys.time()
  rows <- vector("list", nrow(grid))
  for (i in seq_len(nrow(grid))) {
    g <- grid[i, , drop = FALSE]
    sim_one <- function(s) {
      d <- simulate_trial(
        n_per_arm = g$n_per_arm, hr = g$hr,
        accrual_months = g$accrual_months,
        followup_months = g$followup_months,
        dropout_rate = g$dropout_rate,
        median_control = g$median_control,
        alpha = alpha, seed = seed + i * 100000L + s
      )
      fit <- survival::coxph(
        survival::Surv(time, status) ~ arm,
        data = d, ties = "efron"
      )
      if (identical(test, "cox")) {
        pval <- summary(fit)$coefficients[1, "Pr(>|z|)"]
      } else {
        sd <- survival::survdiff(survival::Surv(time, status) ~ arm, data = d)
        pval <- stats::pchisq(sd$chisq, df = 1L, lower.tail = FALSE)
      }
      c(
        p_value = as.numeric(pval),
        events = sum(d$status),
        estimated_hr = exp(unname(stats::coef(fit)[1]))
      )
    }
    if (.Platform$OS.type == "unix" && cores > 1L && n_sims > 1L) {
      res <- parallel::mclapply(seq_len(n_sims), sim_one,
        mc.cores = min(cores, n_sims), mc.set.seed = FALSE
      )
    } else {
      res <- lapply(seq_len(n_sims), sim_one)
    }
    mat <- do.call(rbind, res)
    reject <- mean(mat[, "p_value"] < alpha, na.rm = TRUE)
    se <- sqrt(reject * (1 - reject) / n_sims)
    ci <- .wilson_ci(sum(mat[, "p_value"] < alpha, na.rm = TRUE), n_sims)
    rows[[i]] <- tibble::tibble(
      n_per_arm = as.integer(g$n_per_arm),
      hr = as.numeric(g$hr),
      accrual_months = as.numeric(g$accrual_months),
      followup_months = as.numeric(g$followup_months),
      dropout_rate = as.numeric(g$dropout_rate),
      median_control = as.numeric(g$median_control),
      test = test,
      alpha = alpha,
      n_sims = as.integer(n_sims),
      rejection_rate = reject,
      mc_se = se,
      mc_lower = ci[1],
      mc_upper = ci[2],
      mean_events = mean(mat[, "events"], na.rm = TRUE),
      mean_hr = mean(mat[, "estimated_hr"], na.rm = TRUE),
      mean_log_hr = mean(log(mat[, "estimated_hr"]), na.rm = TRUE)
    )
  }
  out <- dplyr::bind_rows(rows)
  attr(out, "runtime") <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  attr(out, "test") <- test
  attr(out, "seed") <- as.integer(seed)
  out
}

#' Wilson score interval for a binomial proportion
#' @noRd
.wilson_ci <- function(x, n, conf = 0.95) {
  p <- x / n
  z <- stats::qnorm(1 - (1 - conf) / 2)
  denom <- 1 + z^2 / n
  centre <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(max(0, centre - half), min(1, centre + half))
}

#' Analytic per-arm sample size using the Schoenfeld formula
#'
#' Computes the number of events required for a two-arm log-rank test with
#' equal allocation from the Schoenfeld formula, then converts events to
#' per-arm sample size using the probability of observing an event under
#' uniform accrual, exponential survival and random dropout. The event
#' probability is integrated numerically over the accrual window.
#'
#' @param hr Hazard ratio.
#' @param alpha Two-sided significance level.
#' @param power Target power.
#' @param median_control Median control survival in months.
#' @param accrual_months Accrual duration.
#' @param followup_months Additional follow-up.
#' @param dropout_rate Proportion dropping out within `followup_months`.
#'
#' @return A list with `events` (required events), `n_per_arm`, and the
#'   control and treatment event probabilities.
#' @keywords internal
schoenfeld_n <- function(hr, alpha = 0.05, power = 0.8, median_control = 60,
                         accrual_months = 24, followup_months = 60,
                         dropout_rate = 0) {
  if (!is.numeric(hr) || hr <= 0) {
    rlang::abort("`hr` must be positive.", class = "trialsurv_input_error")
  }
  z_a <- stats::qnorm(1 - alpha / 2)
  z_b <- stats::qnorm(power)
  events <- 4 * (z_a + z_b)^2 / (log(hr)^2)
  drop_hazard <- if (dropout_rate > 0) {
    -log(1 - dropout_rate) / followup_months
  } else {
    0
  }
  analysis_time <- accrual_months + followup_months
  event_prob <- function(rate) {
    f <- function(a) {
      mf <- analysis_time - a
      rate / (rate + drop_hazard) * (1 - exp(-(rate + drop_hazard) * mf))
    }
    stats::integrate(f, 0, accrual_months)$value / accrual_months
  }
  rate_c <- log(2) / median_control
  pc <- event_prob(rate_c)
  pt <- event_prob(rate_c * hr)
  n_per_arm <- ceiling(events / (pc + pt))
  list(
    events = events,
    n_per_arm = n_per_arm,
    p_event_control = pc,
    p_event_treatment = pt,
    median_control = median_control,
    accrual_months = accrual_months,
    followup_months = followup_months,
    dropout_rate = dropout_rate
  )
}

#' Find the per-arm sample size required for a target power
#'
#' Searches a grid of per-arm sample sizes by simulation, then compares the
#' smallest simulated size meeting the target with the analytic Schoenfeld
#' sample size.
#'
#' @param target_power Target power (default 0.8).
#' @param hr Hazard ratio.
#' @param n_range Candidate per-arm sample sizes.
#' @param n_sims Simulations per candidate size.
#' @param test `"logrank"` or `"cox"`.
#' @param accrual_months,followup_months,dropout_rate,median_control Design
#'   parameters passed to [simulate_trial()].
#' @param alpha Significance level.
#' @param cores Cores for [trial_power()].
#' @param seed Seed for the simulation search.
#'
#' @return A list of class `trialsurv_required_n` with `search` (the simulated
#'   power curve) and `summary` (a one-row tibble with the simulated and
#'   Schoenfeld per-arm sizes).
#'
#' @examples
#' \dontrun{
#' required_n(hr = 0.7, n_range = c(100, 200), n_sims = 100)
#' }
#'
#' @export
required_n <- function(target_power = 0.8, hr, n_range = c(
                         100, 150, 200, 250, 300, 400, 500, 600, 800
                       ),
                       n_sims = 1000, test = c("logrank", "cox"),
                       accrual_months = 24, followup_months = 60,
                       dropout_rate = 0, median_control = 60, alpha = 0.05,
                       cores = NULL, seed = 20261117) {
  test <- match.arg(test)
  if (!is.numeric(target_power) || target_power <= 0 || target_power >= 1) {
    rlang::abort("`target_power` must be in (0, 1).",
      class = "trialsurv_input_error"
    )
  }
  grid <- expand.grid(
    n_per_arm = n_range, hr = hr, dropout_rate = dropout_rate,
    KEEP.OUT.ATTRS = FALSE
  )
  grid$accrual_months <- accrual_months
  grid$followup_months <- followup_months
  grid$median_control <- median_control
  search <- trial_power(
    grid, n_sims = n_sims, test = test, cores = cores,
    seed = seed, alpha = alpha
  )
  search <- search[order(search$n_per_arm), , drop = FALSE]
  hit <- search[search$rejection_rate >= target_power, , drop = FALSE]
  n_sim <- if (nrow(hit) > 0L) hit$n_per_arm[1] else NA_integer_
  power_at <- if (nrow(hit) > 0L) hit$rejection_rate[1] else NA_real_
  sch <- schoenfeld_n(
    hr = hr, alpha = alpha, power = target_power,
    median_control = median_control, accrual_months = accrual_months,
    followup_months = followup_months, dropout_rate = dropout_rate
  )
  summary_tbl <- tibble::tibble(
    target_power = target_power,
    hr = hr,
    dropout_rate = dropout_rate,
    accrual_months = accrual_months,
    followup_months = followup_months,
    median_control = median_control,
    n_simulated = n_sim,
    power_at_simulated = power_at,
    n_schoenfeld = sch$n_per_arm,
    schoenfeld_events = sch$events,
    p_event_control = sch$p_event_control,
    p_event_treatment = sch$p_event_treatment,
    ratio_sim_over_schoenfeld = n_sim / sch$n_per_arm,
    n_sims = as.integer(n_sims),
    test = test
  )
  out <- list(search = search, summary = summary_tbl)
  class(out) <- c("trialsurv_required_n", "list")
  out
}
