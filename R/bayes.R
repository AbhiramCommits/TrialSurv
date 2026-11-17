#' Bayesian hierarchical survival model with partial borrowing
#'
#' Fits a parametric Weibull proportional hazards model in which the treatment
#' effects of the non-reference arms are allowed to vary across prognostic
#' strata, with a hierarchical shrinkage prior that partially pools each
#' stratum-specific effect toward the corresponding global (population-level)
#' treatment effect. This is the "partial borrowing" mechanism: small or
#' unstable subgroup estimates are pulled toward the overall effect instead of
#' being reported as if they were estimated independently.
#'
#' @details
#' Let \eqn{h(t)} be the hazard of death. The default model is
#' \deqn{h(t \mid arm, stratum) = k t^{k-1} \exp(\eta),}
#' \deqn{\eta = \beta_0 + \sum_a \beta_a 1[arm=a]
#'   + \gamma' z + \sum_a u_{a,stratum} 1[arm=a],}
#' where \eqn{k} is the Weibull shape, \eqn{z} are adjustment covariates,
#' \eqn{\beta_a} is the global log hazard ratio of arm \eqn{a} versus the
#' reference arm, and \eqn{u_{a,s} \sim N(0, \tau_a^2)} is the
#' stratum-specific deviation. The \eqn{\tau_a} are given half-normal
#' \eqn{N^+(0, prior\_scale)} priors, so a stratum with little information is
#' shrunk strongly toward \eqn{\beta_a}. The strata are nodal-risk groups
#' derived from the number of positive lymph nodes (0-1, 2-3, 4-7, 8+), a
#' strong prognostic factor in colon cancer.
#'
#' Priors are weakly informative on the log-hazard scale: normal(0, 2.5) for
#' covariate and treatment coefficients, normal(-5, 2) for the baseline
#' log-hazard parameter in the fallback sampler (normal(4, 2.5) for the
#' log-scale intercept in the brms AFT parameterisation), gamma(2, 2) for the
#' Weibull shape, and half-normal(0, prior_scale) for the random-effect
#' standard deviations. Right censoring enters through the likelihood (the
#' `Surv`/`cens` response); no imputation is performed.
#'
#' The `"brms"` backend uses brms + rstan (or cmdstanr if configured) and
#' brms's non-centred parameterisation. If brms/rstan cannot be installed, the
#' self-contained `"mcmc"` backend fits the same model with a
#' Metropolis-within-Gibbs sampler written in base R: adaptive
#' random-walk Metropolis blocks for the fixed effects, the non-centred
#' random effects, the random-effect standard deviations and the log shape,
#' over the exact Weibull likelihood above. Both backends return posterior
#' draws, R-hat and bulk/tail ESS.
#'
#' @param data A tidy trial data frame with the `death` endpoint, usually from
#'   [load_trial_data()].
#' @param formula Optional custom brms formula. If `NULL` (recommended), the
#'   hierarchical model described above is constructed. Custom formulas are
#'   supported for the `"brms"` backend but disable the automatic hazard-ratio
#'   and shrinkage extraction, because the subgroup structure is unknown.
#' @param family Character; only `"weibull"` is currently implemented.
#' @param subgroup Character naming the treatment factor whose effects are
#'   partially pooled. Defaults to `"rx"`.
#' @param stratum Optional character naming a factor column that defines the
#'   borrowing strata. If `NULL`, nodal-risk strata are derived from `nodes`.
#' @param covariates Character vector of adjustment covariates. Set to
#'   `character(0)` for an unadjusted model.
#' @param pool Which treatment arms receive the hierarchical subgroup random
#'   effect. `"primary"` (default) pools the last level of `subgroup`, which is
#'   the experimental arm of interest; `"all"` pools every non-reference arm
#'   (a more weakly identified model).
#' @param prior_scale Scale of the half-normal prior on the random-effect
#'   standard deviations. Smaller values borrow more aggressively.
#' @param iter Total MCMC iterations per chain, including warm-up.
#' @param chains Number of MCMC chains.
#' @param seed Integer seed.
#' @param backend `"brms"` (default) or `"mcmc"` for the self-contained
#'   fallback sampler.
#' @param adapt_delta Target acceptance probability for the brms/Stan sampler.
#' @param max_treedepth Maximum tree depth for the brms/Stan sampler.
#' @param cores Number of cores for parallel chains (defaults to `chains`).
#' @param ... Further arguments passed to [brms::brm()].
#'
#' @return An object of class `trialsurv_bayes` containing the backend fit,
#'   posterior draws (a `posterior` draws array), derived subgroup hazard-ratio
#'   draws, the analysis data, and convergence metadata.
#'
#' @references
#' Bürkner P-C (2017). brms: An R Package for Bayesian Multilevel Models Using
#' Stan. *Journal of Statistical Software*, 80(1), 1-28.
#'
#' Gelman A, Carlin JB, Stern HS, Dunson DB, Vehtari A, Rubin DB (2013).
#' *Bayesian Data Analysis*, 3rd edition. CRC Press.
#'
#' @examples
#' \dontrun{
#' d <- load_trial_data()
#' sp <- trial_split(d, seed = 20261117)
#' fit <- bayes_surv_fit(sp$train, iter = 1000, chains = 2, seed = 1)
#' bayes_summary(fit)
#' }
#'
#' @export
bayes_surv_fit <- function(data,
                           formula = NULL,
                           family = "weibull",
                           subgroup = "rx",
                           stratum = NULL,
                           covariates = c("age", "sex", "extent", "differ", "obstruct"),
                           pool = c("primary", "all"),
                           prior_scale = 1,
                           iter = 4000,
                           chains = 4,
                           seed = 20261117,
                           backend = c("brms", "mcmc"),
                           adapt_delta = 0.95,
                           max_treedepth = 12,
                           cores = NULL,
                           ...) {
  family <- match.arg(family, "weibull")
  backend <- match.arg(backend)
  pool <- match.arg(pool)
  if (!is.numeric(prior_scale) || length(prior_scale) != 1L || prior_scale <= 0) {
    rlang::abort("`prior_scale` must be a positive number.",
      class = "trialsurv_input_error"
    )
  }
  if (!is.numeric(iter) || length(iter) != 1L || iter < 100) {
    rlang::abort("`iter` must be at least 100.", class = "trialsurv_input_error")
  }
  if (!is.numeric(chains) || length(chains) != 1L || chains < 1) {
    rlang::abort("`chains` must be a positive integer.", class = "trialsurv_input_error")
  }
  if (!requireNamespace("posterior", quietly = TRUE)) {
    rlang::abort("Package 'posterior' is required for Bayesian summaries.",
      class = "trialsurv_input_error"
    )
  }
  if (is.null(cores)) cores <- chains

  prep <- .prepare_bayes_data(data, subgroup = subgroup, stratum = stratum,
    covariates = covariates, pool = pool
  )
  d <- prep$data

  if (is.null(formula)) {
    form <- prep$brms_formula
    custom_formula <- FALSE
  } else {
    if (!inherits(formula, "formula")) {
      rlang::abort("`formula` must be a formula or NULL.",
        class = "trialsurv_input_error"
      )
    }
    if (identical(backend, "mcmc")) {
      rlang::abort(
        "Custom formulas are only supported with backend = 'brms'.",
        class = "trialsurv_input_error"
      )
    }
    form <- formula
    custom_formula <- TRUE
  }

  start <- Sys.time()
  if (identical(backend, "brms")) {
    if (!requireNamespace("brms", quietly = TRUE)) {
      rlang::abort(
        "Package 'brms' is required for backend = 'brms'; use backend = 'mcmc'.",
        class = "trialsurv_input_error"
      )
    }
    fit_obj <- .fit_bayes_brms(
      form = form, data = d, prior_scale = prior_scale, iter = iter,
      chains = chains, seed = seed, adapt_delta = adapt_delta,
      max_treedepth = max_treedepth, cores = cores, ...
    )
    draws <- posterior::as_draws_array(fit_obj$fit)
    raw_fit <- fit_obj$fit
    n_divergent <- fit_obj$n_divergent
    attempts <- fit_obj$attempts
  } else {
    fit_obj <- .fit_bayes_mcmc(
      prep = prep, iter = iter, chains = chains, seed = seed,
      prior_scale = prior_scale, cores = cores
    )
    draws <- fit_obj$draws
    raw_fit <- fit_obj
    n_divergent <- 0L
    attempts <- fit_obj$attempts
  }
  elapsed <- as.numeric(difftime(Sys.time(), start, units = "secs"))

  hr_draws <- NULL
  if (!custom_formula) {
    hr_draws <- .derive_hr_draws(
      draws = draws, prep = prep, backend = backend
    )
  }

  out <- list(
    fit = raw_fit,
    draws = draws,
    hr_draws = hr_draws,
    data = d,
    prep = prep,
    formula = form,
    family = family,
    subgroup = subgroup,
    pool = pool,
    stratum = if (is.null(stratum)) ".nodal_stratum" else stratum,
    covariates = covariates,
    prior_scale = prior_scale,
    iter = iter,
    chains = chains,
    seed = as.integer(seed),
    backend = backend,
    custom_formula = custom_formula,
    n_divergent = as.integer(n_divergent),
    attempts = attempts,
    fit_seconds = elapsed
  )
  class(out) <- c("trialsurv_bayes", "list")
  out
}

# ---- data preparation -------------------------------------------------------

#' Derive four nodal-risk strata from the number of positive nodes
#' @noRd
.nodal_stratum <- function(nodes) {
  cut(
    nodes,
    breaks = c(-Inf, 1, 3, 7, Inf),
    labels = c("n01", "n23", "n47", "n8p"),
    include.lowest = TRUE,
    right = TRUE
  )
}

#' Prepare data, dummies and design matrices for the Bayesian model
#' @noRd
.prepare_bayes_data <- function(data, subgroup, stratum, covariates,
                                pool = c("primary", "all")) {
  pool <- match.arg(pool)
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.", class = "trialsurv_input_error")
  }
  if (!all(c("time", "status", "endpoint", "nodes") %in% names(data))) {
    rlang::abort(
      "`data` must be output from load_trial_data() with time, status, endpoint and nodes.",
      class = "trialsurv_input_error"
    )
  }
  d <- data[data$endpoint == "death", , drop = FALSE]
  if (nrow(d) == 0L) {
    rlang::abort("No death-endpoint rows in `data`.", class = "trialsurv_input_error")
  }
  d <- as.data.frame(d)
  assert_surv_inputs(d, vars = unique(c(subgroup, covariates,
    if (!is.null(stratum)) stratum else "nodes"
  )))
  d <- d[stats::complete.cases(d[, unique(c(
    "time", "status", subgroup, covariates,
    if (!is.null(stratum)) stratum else "nodes"
  ))]), , drop = FALSE]

  if (!subgroup %in% names(d)) {
    rlang::abort(sprintf("Subgroup column `%s` not found.", subgroup),
      class = "trialsurv_input_error"
    )
  }
  if (!is.factor(d[[subgroup]])) d[[subgroup]] <- factor(d[[subgroup]])
  d[[subgroup]] <- droplevels(d[[subgroup]])
  arm_levels <- levels(d[[subgroup]])
  if (length(arm_levels) < 2L) {
    rlang::abort("The subgroup/treatment factor needs at least two levels.",
      class = "trialsurv_input_error"
    )
  }
  ref_arm <- arm_levels[1]
  nonref <- arm_levels[-1]

  # Readable dummy names for the colon trial arms, generic names otherwise.
  dummy_names <- vapply(nonref, function(a) {
    if (identical(subgroup, "rx") && identical(a, "Lev")) {
      "lev"
    } else if (identical(subgroup, "rx") && identical(a, "Lev+5FU")) {
      "f5"
    } else {
      paste0("arm", match(a, nonref))
    }
  }, character(1))
  dummy_names <- make.unique(dummy_names)
  for (j in seq_along(nonref)) {
    d[[dummy_names[j]]] <- as.integer(d[[subgroup]] == nonref[j])
  }

  if (is.null(stratum)) {
    d$stratum <- .nodal_stratum(d$nodes)
  } else {
    if (!stratum %in% names(d)) {
      rlang::abort(sprintf("Stratum column `%s` not found.", stratum),
        class = "trialsurv_input_error"
      )
    }
    d$stratum <- d[[stratum]]
    if (!is.factor(d$stratum)) d$stratum <- factor(d$stratum)
  }
  d$stratum <- droplevels(factor(d$stratum))
  stratum_levels <- levels(d$stratum)
  if (length(stratum_levels) < 2L) {
    rlang::abort("At least two borrowing strata are required.",
      class = "trialsurv_input_error"
    )
  }

  for (v in covariates) {
    if (is.character(d[[v]])) d[[v]] <- factor(d[[v]])
    if (is.factor(d[[v]])) d[[v]] <- droplevels(d[[v]])
  }

  pooled_flags <- rep(FALSE, length(dummy_names))
  if (identical(pool, "primary")) {
    pooled_flags[length(pooled_flags)] <- TRUE
  } else {
    pooled_flags[] <- TRUE
  }
  arms <- tibble::tibble(
    level = nonref,
    dummy = dummy_names,
    coef = paste0("b_", dummy_names),
    pooled = pooled_flags
  )

  rhs <- c(dummy_names, covariates, "stratum")
  fixed_formula <- stats::reformulate(rhs)
  X <- stats::model.matrix(fixed_formula, data = d)

  ran_dummies <- dummy_names[pooled_flags]
  ran_names <- unlist(lapply(ran_dummies, function(dd) {
    paste(dd, stratum_levels, sep = "_")
  }), use.names = FALSE)
  ran_arm <- rep(seq_along(ran_dummies), each = length(stratum_levels))
  Z <- matrix(0, nrow = nrow(d), ncol = length(ran_names))
  colnames(Z) <- paste0("u_", ran_names)
  for (j in seq_along(ran_dummies)) {
    for (s in seq_along(stratum_levels)) {
      k <- (j - 1L) * length(stratum_levels) + s
      Z[, k] <- as.integer(
        d[[ran_dummies[j]]] == 1L & d$stratum == stratum_levels[s]
      )
    }
  }

  ran_term <- paste0("(0 + ", paste(ran_dummies, collapse = " + "), " || stratum)")
  brms_formula <- stats::as.formula(paste(
    "time | cens(1 - status) ~", paste(c(rhs, ran_term), collapse = " + ")
  ))

  factor_levels <- lapply(d[c(covariates, "stratum")], function(x) {
    if (is.factor(x)) levels(x) else NULL
  })
  factor_levels <- factor_levels[!vapply(factor_levels, is.null, logical(1))]

  list(
    data = d,
    subgroup = subgroup,
    arm_levels = arm_levels,
    ref_arm = ref_arm,
    arms = arms,
    stratum_levels = stratum_levels,
    covariates = covariates,
    fixed_formula = fixed_formula,
    X = X,
    Z = Z,
    ran_names = ran_names,
    ran_arm = ran_arm,
    brms_formula = brms_formula,
    factor_levels = factor_levels,
    dummy_names = dummy_names
  )
}

# ---- brms backend -----------------------------------------------------------

#' Fit the hierarchical Weibull model with brms
#' @noRd
.fit_bayes_brms <- function(form, data, prior_scale, iter, chains, seed,
                            adapt_delta, max_treedepth, cores, ...) {
  stan_backend <- "rstan"
  priors <- c(
    brms::prior_string("normal(0, 2.5)", class = "b"),
    brms::prior_string("normal(4, 2.5)", class = "Intercept"),
    brms::prior_string("gamma(2, 2)", class = "shape"),
    brms::prior_string(sprintf("normal(0, %.10g)", prior_scale), class = "sd")
  )

  attempts <- tibble::tibble(
    attempt = integer(), adapt_delta = numeric(), n_divergent = integer(),
    max_rhat = numeric(), seconds = numeric()
  )
  current_delta <- adapt_delta
  fit <- NULL
  best_fit <- NULL
  best_score <- Inf
  best_n_div <- NA_integer_
  for (attempt in seq_len(2L)) {
    t0 <- Sys.time()
    fit <- brms::brm(
      formula = form, data = data, family = brms::weibull(),
      prior = priors, iter = iter, chains = chains, cores = cores,
      seed = seed, refresh = 0, silent = 2, backend = stan_backend,
      control = list(adapt_delta = current_delta, max_treedepth = max_treedepth),
      ...
    )
    np <- brms::nuts_params(fit)
    n_div <- sum(np$Value[np$Parameter == "divergent__"], na.rm = TRUE)
    rhats <- posterior::rhat(posterior::as_draws_array(fit))
    max_r <- suppressWarnings(max(rhats, na.rm = TRUE))
    attempts <- dplyr::bind_rows(attempts, tibble::tibble(
      attempt = attempt,
      adapt_delta = current_delta,
      n_divergent = as.integer(n_div),
      max_rhat = max_r,
      seconds = as.numeric(difftime(Sys.time(), t0, units = "secs"))
    ))
    # Keep the attempt with the best combination of R-hat and divergences.
    score <- max_r + 0.001 * n_div
    if (is.finite(score) && score < best_score) {
      best_score <- score
      best_fit <- fit
      best_n_div <- as.integer(n_div)
    }
    if (n_div == 0L) break
    # A second, slightly more conservative sampler configuration is worth
    # trying only when the chains otherwise mixed (R-hat <= 1.01). If R-hat
    # is already poor, more iterations rather than a higher target acceptance
    # rate are needed, so the script handles that re-fit.
    if (max_r > 1.01 || attempt >= 2L) break
    current_delta <- min(0.97, current_delta + 0.02)
    message(sprintf(
      "Bayesian fit attempt %d had %d divergent transition(s); refitting with adapt_delta = %.2f.",
      attempt, n_div, current_delta
    ))
  }

  list(fit = best_fit, n_divergent = best_n_div, attempts = attempts)
}

# ---- self-contained MCMC backend -------------------------------------------

#' Adaptive random-walk Metropolis-within-Gibbs sampler for the Weibull model
#' @noRd
.fit_bayes_mcmc <- function(prep, iter, chains, seed, prior_scale, cores) {
  X <- prep$X
  Z <- prep$Z
  time <- prep$data$time
  status <- prep$data$status
  logt <- log(time)
  p <- ncol(X)
  q <- ncol(Z)
  n_arms <- sum(prep$arms$pooled)
  n_strata <- length(prep$stratum_levels)
  arm_of_ran <- prep$ran_arm

  warmup <- max(100L, floor(iter * 0.5))
  nsave <- iter - warmup
  if (nsave < 20L) {
    rlang::abort("`iter` leaves fewer than 20 post-warm-up draws.",
      class = "trialsurv_input_error"
    )
  }

  # Starting values from a Cox fit and an exponential hazard estimate.
  beta0 <- rep(0, p)
  beta0[1] <- log(max(sum(status), 1) / sum(time))
  cox_start <- tryCatch(
    survival::coxph(
      stats::as.formula(paste(
        "survival::Surv(time, status) ~",
        paste(colnames(X)[-1], collapse = " + ")
      )),
      data = prep$data
    ),
    error = function(e) NULL
  )
  if (!is.null(cox_start)) {
    cf <- stats::coef(cox_start)
    idx <- match(names(cf), colnames(X)[-1])
    ok <- !is.na(idx) & is.finite(cf)
    if (any(ok)) beta0[-1][idx[ok]] <- as.numeric(cf)[ok]
  }

  # z is stored arm-major, matching the columns of Z. Within each arm the
  # random effects are centred to sum to zero, so the global treatment effect
  # and the mean subgroup deviation are identified separately.
  centre_z <- function(z) {
    zm <- matrix(z, nrow = n_strata, ncol = n_arms)
    zm <- zm - matrix(colMeans(zm), nrow = n_strata, ncol = n_arms,
      byrow = TRUE
    )
    as.vector(zm)
  }

  log_post <- function(theta, log_sigma, log_k) {
    beta <- theta[seq_len(p)]
    z <- theta[p + seq_len(q)]
    sigma <- exp(log_sigma)
    zc <- centre_z(z)
    eta <- as.vector(X %*% beta + Z %*% (zc * sigma[arm_of_ran]))
    k <- exp(log_k)
    logh <- log_k + (k - 1) * logt + eta
    H <- exp(k * logt + eta)
    ll <- sum(status * logh - H)
    if (!is.finite(ll)) {
      return(-Inf)
    }
    lp <- stats::dnorm(beta[1], -5, 2, log = TRUE) +
      sum(stats::dnorm(beta[-1], 0, 2.5, log = TRUE)) +
      sum(stats::dnorm(z, 0, 1, log = TRUE)) +
      sum(stats::dnorm(sigma, 0, prior_scale, log = TRUE) + log(2) + log_sigma) +
      stats::dnorm(log_k, 0, 0.5, log = TRUE)
    ll + lp
  }

  run_chain <- function(ch) {
    set.seed(seed + ch)
    theta <- c(beta0, rep(0, q))
    log_sigma <- rep(log(0.25), n_arms)
    log_k <- 0

    blocks <- list(
      theta = list(d = p + q, target = 0.25, init = 0.01),
      log_sigma = list(d = n_arms, target = 0.3, init = 0.05),
      log_k = list(d = 1L, target = 0.44, init = 0.05)
    )
    S <- lapply(blocks, function(b) diag(b$init, b$d))
    L <- lapply(S, function(s) t(chol(s)))
    mult <- lapply(blocks, function(b) 1)
    history <- lapply(blocks, function(b) {
      matrix(NA_real_, nrow = warmup, ncol = b$d)
    })
    accept <- lapply(blocks, function(b) logical(warmup))
    out <- matrix(NA_real_, nrow = nsave, ncol = p + q + n_arms + 1L)

    cur <- list(theta = theta, log_sigma = log_sigma, log_k = log_k)
    cur_lp <- log_post(cur$theta, cur$log_sigma, cur$log_k)

    for (it in seq_len(iter)) {
      for (nm in names(blocks)) {
        b <- blocks[[nm]]
        prop <- cur[[nm]] +
          as.vector(mult[[nm]] * L[[nm]] %*% stats::rnorm(b$d))
        prop_list <- cur
        prop_list[[nm]] <- prop
        prop_lp <- log_post(
          prop_list$theta, prop_list$log_sigma, prop_list$log_k
        )
        acc <- is.finite(prop_lp) && log(stats::runif(1)) < prop_lp - cur_lp
        if (acc) {
          cur <- prop_list
          cur_lp <- prop_lp
        }
        if (it <= warmup) {
          history[[nm]][it, ] <- cur[[nm]]
          accept[[nm]][it] <- acc
        }
      }

      if (it > warmup) {
        k <- it - warmup
        sigma <- exp(cur$log_sigma)
        u <- centre_z(cur$theta[p + seq_len(q)]) * sigma[arm_of_ran]
        out[k, ] <- c(
          cur$theta[seq_len(p)], u, sigma, exp(cur$log_k)
        )
      }

      if (it <= warmup && it >= 100 && it %% 25 == 0) {
        for (nm in names(blocks)) {
          b <- blocks[[nm]]
          S_new <- stats::cov(history[[nm]][seq_len(it), , drop = FALSE])
          S_new <- S_new + diag(1e-6, b$d)
          chk <- tryCatch(t(chol(S_new)), error = function(e) NULL)
          if (!is.null(chk)) {
            S[[nm]] <- S_new
            L[[nm]] <- chk
          }
          recent <- accept[[nm]][max(1, it - 24):it]
          mult[[nm]] <- mult[[nm]] * exp(0.5 * (mean(recent) - b$target))
          mult[[nm]] <- min(max(mult[[nm]], 0.05), 10)
        }
      }
    }

    colnames(out) <- c(
      paste0("b_", colnames(X)),
      paste0("u_", prep$ran_names),
      paste0("sd_", prep$arms$dummy[prep$arms$pooled]),
      "shape"
    )
    if (isTRUE(getOption("trialsurv.debug_mcmc"))) {
      message(sprintf(
        "chain %d acceptance: theta %.2f, log_sigma %.2f, log_k %.2f; mult %.2f/%.2f/%.2f",
        ch, mean(accept$theta), mean(accept$log_sigma), mean(accept$log_k),
        mult$theta, mult$log_sigma, mult$log_k
      ))
    }
    out
  }

  if (.Platform$OS.type == "unix" && cores > 1L) {
    chain_list <- parallel::mclapply(seq_len(chains), run_chain,
      mc.cores = min(cores, chains), mc.set.seed = FALSE
    )
  } else {
    chain_list <- lapply(seq_len(chains), run_chain)
  }
  if (any(vapply(chain_list, function(x) !is.matrix(x), logical(1)))) {
    rlang::abort("A fallback MCMC chain failed; see the warnings above.",
      class = "trialsurv_input_error"
    )
  }

  nparam <- ncol(chain_list[[1]])
  arr <- array(NA_real_, dim = c(nsave, chains, nparam))
  dimnames(arr)[[3]] <- colnames(chain_list[[1]])
  for (ch in seq_len(chains)) arr[, ch, ] <- chain_list[[ch]]

  list(
    draws = posterior::as_draws_array(arr),
    n_divergent = 0L,
    attempts = tibble::tibble(
      attempt = 1L, adapt_delta = NA_real_, n_divergent = 0L,
      max_rhat = suppressWarnings(max(posterior::rhat(arr), na.rm = TRUE)),
      seconds = NA_real_
    )
  )
}

# ---- derived hazard ratios --------------------------------------------------

#' Derive subgroup and global hazard-ratio draws
#' @noRd
.derive_hr_draws <- function(draws, prep, backend) {
  arr <- posterior::as_draws_array(draws)
  vars <- dimnames(arr)[[3]]
  it <- dim(arr)[1]
  ch <- dim(arr)[2]
  hr_list <- list()

  for (j in seq_len(nrow(prep$arms))) {
    dummy <- prep$arms$dummy[j]
    bname <- paste0("b_", dummy)
    if (!bname %in% vars) next
    bdraw <- arr[, , bname]
    if (identical(backend, "brms")) {
      # brms Weibull is an AFT parameterisation: the linear predictor is on
      # the log scale, so the log hazard ratio is minus the coefficient.
      log_global <- -bdraw
    } else {
      log_global <- bdraw
    }
    hr_list[[paste0("hr_", dummy, "_global")]] <- exp(log_global)
    if (!isTRUE(prep$arms$pooled[j])) next
    for (s in prep$stratum_levels) {
      rname <- sprintf("r_stratum[%s,%s]", s, dummy)
      if (rname %in% vars) {
        rdraw <- arr[, , rname]
        log_s <- if (identical(backend, "brms")) {
          -(bdraw + rdraw)
        } else {
          bdraw + rdraw
        }
      } else {
        uname <- paste0("u_", dummy, "_", s)
        if (uname %in% vars) {
          log_s <- if (identical(backend, "brms")) {
            -(bdraw + arr[, , uname])
          } else {
            bdraw + arr[, , uname]
          }
        } else {
          log_s <- log_global
        }
      }
      hr_list[[paste0("hr_", dummy, "_", s)]] <- exp(log_s)
    }
  }
  if (length(hr_list) == 0L) {
    return(NULL)
  }
  out <- array(NA_real_, dim = c(it, ch, length(hr_list)),
    dimnames = list(NULL, NULL, names(hr_list))
  )
  for (nm in names(hr_list)) out[, , nm] <- hr_list[[nm]]
  posterior::as_draws_array(out)
}

# ---- summaries --------------------------------------------------------------

#' Tidy posterior summary of a Bayesian survival fit
#'
#' Summarises every population-level and group-level parameter of a
#' `trialsurv_bayes` fit with the posterior mean, median, standard deviation,
#' 95% credible interval, split R-hat and bulk/tail effective sample size. A
#' second block of rows contains the exponentiated hazard ratios for each
#' treatment arm globally and within each borrowing stratum.
#'
#' @param fit A `trialsurv_bayes` fit from [bayes_surv_fit()].
#' @param probs Numeric vector of two probabilities for the credible interval.
#'
#' @return A tibble with columns `parameter`, `type`, `mean`, `median`, `sd`,
#'   `q2.5`, `q97.5`, `rhat`, `ess_bulk`, `ess_tail`, `hazard_ratio`,
#'   `hr_lower` and `hr_upper`. The attributes `max_rhat`, `min_ess_bulk`,
#'   `min_ess_tail` and `n_divergent` record the convergence checks.
#'
#' @examples
#' \dontrun{
#' fit <- bayes_surv_fit(load_trial_data(), iter = 1000, chains = 2)
#' bayes_summary(fit)
#' }
#'
#' @export
bayes_summary <- function(fit, probs = c(0.025, 0.975)) {
  if (!inherits(fit, "trialsurv_bayes")) {
    rlang::abort("`fit` must be a trialsurv_bayes object.",
      class = "trialsurv_input_error"
    )
  }
  if (length(probs) != 2L || any(probs <= 0) || any(probs >= 1)) {
    rlang::abort("`probs` must be two probabilities in (0, 1).",
      class = "trialsurv_input_error"
    )
  }

  raw <- .summarise_draws(fit$draws, probs)
  raw <- raw[!raw$parameter %in% c("lp__", "lprior"), , drop = FALSE]
  raw$type <- "parameter"
  raw$hazard_ratio <- NA_real_
  raw$hr_lower <- NA_real_
  raw$hr_upper <- NA_real_

  if (!is.null(fit$hr_draws)) {
    hr <- .summarise_draws(fit$hr_draws, probs)
    hr$type <- "hazard_ratio"
    hr$hazard_ratio <- hr$mean
    hr$hr_lower <- hr$q2.5
    hr$hr_upper <- hr$q97.5
    out <- dplyr::bind_rows(raw, hr)
  } else {
    out <- raw
  }
  out <- out[, c(
    "parameter", "type", "mean", "median", "sd", "q2.5", "q97.5",
    "rhat", "ess_bulk", "ess_tail", "hazard_ratio", "hr_lower", "hr_upper"
  )]
  attr(out, "max_rhat") <- suppressWarnings(max(out$rhat, na.rm = TRUE))
  attr(out, "min_ess_bulk") <- suppressWarnings(min(out$ess_bulk, na.rm = TRUE))
  attr(out, "min_ess_tail") <- suppressWarnings(min(out$ess_tail, na.rm = TRUE))
  attr(out, "n_divergent") <- fit$n_divergent
  attr(out, "backend") <- fit$backend
  out
}

#' Summarise a draws array with base R and posterior diagnostics
#' @noRd
.summarise_draws <- function(draws, probs = c(0.025, 0.975)) {
  arr <- posterior::as_draws_array(draws)
  vars <- dimnames(arr)[[3]]
  it <- dim(arr)[1]
  ch <- dim(arr)[2]
  mat <- matrix(arr, nrow = it * ch, ncol = length(vars))
  colnames(mat) <- vars
  rh <- vapply(vars, function(v) {
    suppressWarnings(posterior::rhat(arr[, , v, drop = FALSE]))
  }, numeric(1))
  eb <- vapply(vars, function(v) {
    suppressWarnings(posterior::ess_bulk(arr[, , v, drop = FALSE]))
  }, numeric(1))
  et <- vapply(vars, function(v) {
    suppressWarnings(posterior::ess_tail(arr[, , v, drop = FALSE]))
  }, numeric(1))
  qs <- apply(mat, 2, stats::quantile, probs = probs, names = FALSE)
  tibble::tibble(
    parameter = vars,
    mean = colMeans(mat),
    median = apply(mat, 2, stats::median),
    sd = apply(mat, 2, stats::sd),
    q2.5 = as.numeric(qs[1, ]),
    q97.5 = as.numeric(qs[2, ]),
    rhat = as.numeric(rh),
    ess_bulk = as.numeric(eb),
    ess_tail = as.numeric(et)
  )
}

#' Quantify borrowing and shrinkage for each treatment subgroup
#'
#' Compares, for every treatment arm and borrowing stratum, the standalone
#' (no-pooling) Cox log hazard ratio with the partially pooled posterior mean
#' from [bayes_surv_fit()], the global effect, the absolute movement toward the
#' global effect, and the percent reduction in interval width.
#'
#' @param fit A `trialsurv_bayes` fit from [bayes_surv_fit()] using the default
#'   hierarchical formula.
#' @param freq_fit A `trialsurv_cox` object or a [survival::coxph()] fit on the
#'   same endpoint, used for the global frequentist comparison.
#'
#' @return A tibble with one row per arm-by-stratum subgroup and columns for
#'   the standalone and pooled log hazard ratios, their interval widths, the
#'   absolute shrinkage, the percentage movement toward the global effect and
#'   the percentage reduction in interval width.
#'
#' @examples
#' \dontrun{
#' d <- load_trial_data()
#' sp <- trial_split(d)
#' cx <- cox_fit(sp$train)
#' bf <- bayes_surv_fit(sp$train, iter = 1000, chains = 2)
#' shrinkage_report(bf, cx)
#' }
#'
#' @export
shrinkage_report <- function(fit, freq_fit) {
  if (!inherits(fit, "trialsurv_bayes")) {
    rlang::abort("`fit` must be a trialsurv_bayes object.",
      class = "trialsurv_input_error"
    )
  }
  if (is.null(fit$hr_draws)) {
    rlang::abort(
      "shrinkage_report() requires the default hierarchical formula.",
      class = "trialsurv_input_error"
    )
  }
  if (inherits(freq_fit, "trialsurv_cox")) freq_fit <- freq_fit$adjusted
  if (!inherits(freq_fit, "coxph")) {
    rlang::abort("`freq_fit` must be a trialsurv_cox or coxph object.",
      class = "trialsurv_input_error"
    )
  }

  d <- fit$data
  prep <- fit$prep
  hr_arr <- posterior::as_draws_array(fit$hr_draws)
  arms <- prep$arms
  z <- stats::qnorm(0.975)

  global_freq <- stats::coef(freq_fit)
  freq_names <- names(global_freq)

  rows <- list()
  for (j in which(arms$pooled)) {
    dummy <- arms$dummy[j]
    coef_name <- paste0("rx", arms$level[j])
    global_par <- paste0("hr_", dummy, "_global")
    global_draws <- log(as.numeric(hr_arr[, , global_par]))
    global_mean <- mean(global_draws)
    freq_global <- if (coef_name %in% freq_names) {
      as.numeric(global_freq[[coef_name]])
    } else {
      NA_real_
    }

    for (s in prep$stratum_levels) {
      ds <- d[d$stratum == s, , drop = FALSE]
      standalone <- tryCatch({
        fs <- survival::coxph(
          survival::Surv(time, status) ~ subgroup,
          data = transform(ds, subgroup = ds[[prep$subgroup]]),
          ties = "efron"
        )
        cf <- stats::coef(fs)
        se <- sqrt(diag(stats::vcov(fs)))
        nm <- paste0("subgroup", arms$level[j])
        if (nm %in% names(cf)) {
          c(log_hr = as.numeric(cf[[nm]]), se = as.numeric(se[[nm]]))
        } else {
          c(log_hr = NA_real_, se = NA_real_)
        }
      }, error = function(e) c(log_hr = NA_real_, se = NA_real_))

      pooled_par <- paste0("hr_", dummy, "_", s)
      pooled_draws <- log(as.numeric(hr_arr[, , pooled_par]))
      pooled_mean <- mean(pooled_draws)
      pooled_q <- stats::quantile(pooled_draws, c(0.025, 0.975), names = FALSE)
      standalone_width <- 2 * z * standalone[["se"]]
      pooled_width <- diff(pooled_q)
      movement <- if (is.finite(standalone[["log_hr"]]) &&
        is.finite(global_mean) &&
        abs(standalone[["log_hr"]] - global_mean) > 1e-8) {
        100 * (abs(standalone[["log_hr"]] - global_mean) -
          abs(pooled_mean - global_mean)) /
          abs(standalone[["log_hr"]] - global_mean)
      } else {
        NA_real_
      }
      width_reduction <- if (is.finite(standalone_width) &&
        standalone_width > 0) {
        100 * (1 - pooled_width / standalone_width)
      } else {
        NA_real_
      }

      rows[[length(rows) + 1L]] <- tibble::tibble(
        arm = arms$level[j],
        stratum = s,
        n = nrow(ds),
        n_events = sum(ds$status),
        standalone_log_hr = as.numeric(standalone[["log_hr"]]),
        standalone_se = as.numeric(standalone[["se"]]),
        standalone_hr = exp(as.numeric(standalone[["log_hr"]])),
        standalone_lower = exp(as.numeric(standalone[["log_hr"]]) - z *
          as.numeric(standalone[["se"]])),
        standalone_upper = exp(as.numeric(standalone[["log_hr"]]) + z *
          as.numeric(standalone[["se"]])),
        standalone_width = standalone_width,
        pooled_log_hr = pooled_mean,
        pooled_hr = exp(pooled_mean),
        pooled_lower = exp(pooled_q[1]),
        pooled_upper = exp(pooled_q[2]),
        pooled_width = pooled_width,
        global_log_hr = global_mean,
        global_hr = exp(global_mean),
        freq_global_log_hr = freq_global,
        shrinkage_abs = abs(as.numeric(standalone[["log_hr"]]) - pooled_mean),
        pct_toward_global = movement,
        pct_width_reduction = width_reduction,
        pooled_between = is.finite(standalone[["log_hr"]]) &&
          pooled_mean >= min(standalone[["log_hr"]], global_mean) - 1e-8 &&
          pooled_mean <= max(standalone[["log_hr"]], global_mean) + 1e-8
      )
    }
  }
  out <- dplyr::bind_rows(rows)
  attr(out, "prior_scale") <- fit$prior_scale
  out
}

#' Posterior held-out concordance for a Bayesian survival fit
#'
#' Computes Harrell's C-index on new data from posterior draws of the risk
#' score (the negative linear predictor for the brms AFT parameterisation, or
#' the log hazard for the fallback sampler). The posterior mean and a 95%
#' interval over draws are returned.
#'
#' @param fit A `trialsurv_bayes` fit.
#' @param newdata Held-out data with the same columns as the training data.
#' @param n_draws Number of posterior draws used for the C-index interval.
#'
#' @return A one-row tibble with `model`, `c_index`, `conf.low`, `conf.high`,
#'   `n`, `n_events` and `n_draws`.
#'
#' @examples
#' \dontrun{
#' d <- load_trial_data()
#' sp <- trial_split(d)
#' bf <- bayes_surv_fit(sp$train, iter = 1000, chains = 2)
#' bayes_concordance(bf, sp$test)
#' }
#'
#' @export
bayes_concordance <- function(fit, newdata, n_draws = 200L) {
  if (!inherits(fit, "trialsurv_bayes")) {
    rlang::abort("`fit` must be a trialsurv_bayes object.",
      class = "trialsurv_input_error"
    )
  }
  if (!is.data.frame(newdata)) {
    rlang::abort("`newdata` must be a data frame.", class = "trialsurv_input_error")
  }
  nd <- newdata[newdata$endpoint == "death", , drop = FALSE]
  if (nrow(nd) == 0L) {
    rlang::abort("No death-endpoint rows in `newdata`.",
      class = "trialsurv_input_error"
    )
  }
  prep_nd <- .prepare_newdata(fit, nd)

  if (identical(fit$backend, "brms")) {
    lp <- brms::posterior_linpred(
      fit$fit, newdata = prep_nd$data, re_formula = NULL,
      transform = FALSE
    )
    risk <- -as.matrix(lp)
  } else {
    bn <- grep("^b_", dimnames(fit$draws)[[3]])
    beta <- matrix(as.numeric(fit$draws[, , bn]), nrow = dim(fit$draws)[1] *
      dim(fit$draws)[2], ncol = length(bn))
    un <- grep("^u_", dimnames(fit$draws)[[3]])
    u <- matrix(as.numeric(fit$draws[, , un]), nrow = dim(fit$draws)[1] *
      dim(fit$draws)[2], ncol = length(un))
    eta <- beta %*% t(prep_nd$X) + u %*% t(prep_nd$Z)
    risk <- eta
  }

  n_take <- min(as.integer(n_draws), nrow(risk))
  idx <- unique(round(seq(1, nrow(risk), length.out = n_take)))
  cvals <- vapply(idx, function(i) {
    dd <- prep_nd$data
    dd$.risk_score <- risk[i, ]
    tryCatch(
      as.numeric(survival::concordance(
        survival::Surv(time, status) ~ .risk_score,
        data = dd, reverse = TRUE
      )$concordance),
      error = function(e) NA_real_
    )
  }, numeric(1))
  cvals <- cvals[is.finite(cvals)]
  if (length(cvals) == 0L) {
    rlang::abort("Could not compute a posterior concordance.",
      class = "trialsurv_input_error"
    )
  }
  tibble::tibble(
    model = "bayes",
    c_index = mean(cvals),
    conf.low = as.numeric(stats::quantile(cvals, 0.025)),
    conf.high = as.numeric(stats::quantile(cvals, 0.975)),
    n = nrow(prep_nd$data),
    n_events = sum(prep_nd$data$status == 1L),
    n_draws = length(cvals)
  )
}

#' Align new data with the training design matrix
#' @noRd
.prepare_newdata <- function(fit, nd) {
  prep <- fit$prep
  nd <- as.data.frame(nd)
  if (is.null(nd$stratum)) {
    nd$stratum <- .nodal_stratum(nd$nodes)
  }
  nd$stratum <- factor(nd$stratum, levels = prep$stratum_levels)
  if (!is.factor(nd[[prep$subgroup]])) {
    nd[[prep$subgroup]] <- factor(nd[[prep$subgroup]])
  }
  nd[[prep$subgroup]] <- factor(nd[[prep$subgroup]],
    levels = prep$arm_levels
  )
  for (v in names(prep$factor_levels)) {
    if (v %in% names(nd)) {
      nd[[v]] <- factor(nd[[v]], levels = prep$factor_levels[[v]])
    }
  }
  for (j in seq_len(nrow(prep$arms))) {
    nd[[prep$arms$dummy[j]]] <- as.integer(
      nd[[prep$subgroup]] == prep$arms$level[j]
    )
  }
  nd <- nd[stats::complete.cases(nd[, unique(c(
    "time", "status", prep$covariates, "stratum",
    prep$dummy_names
  ))]), , drop = FALSE]

  X <- stats::model.matrix(prep$fixed_formula, data = nd)
  missing_cols <- setdiff(colnames(prep$X), colnames(X))
  if (length(missing_cols) > 0L) {
    add <- matrix(0, nrow = nrow(X), ncol = length(missing_cols),
      dimnames = list(NULL, missing_cols)
    )
    X <- cbind(X, add)
  }
  X <- X[, colnames(prep$X), drop = FALSE]

  Z <- matrix(0, nrow = nrow(nd), ncol = length(prep$ran_names))
  colnames(Z) <- colnames(prep$Z)
  pooled_dummies <- prep$arms$dummy[prep$arms$pooled]
  for (j in seq_along(pooled_dummies)) {
    for (s in seq_along(prep$stratum_levels)) {
      k <- (j - 1L) * length(prep$stratum_levels) + s
      Z[, k] <- as.integer(
        nd[[pooled_dummies[j]]] == 1L &
          nd$stratum == prep$stratum_levels[s]
      )
    }
  }
  list(data = nd, X = X, Z = Z)
}
