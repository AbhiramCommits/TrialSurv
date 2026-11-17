#' Forest plot of posterior hazard ratios
#'
#' Plots the posterior mean and 95% credible interval for each treatment arm
#' globally and within each borrowing stratum, faceted by arm, with a vertical
#' reference line at a hazard ratio of 1. The figure is saved with
#' [fig_save()].
#'
#' @param fit A `trialsurv_bayes` fit using the default hierarchical formula.
#' @param name Output file name (with or without `.png`).
#' @param dir Output directory.
#' @param save Logical; save the figure.
#'
#' @return The plot object, invisibly.
#'
#' @examples
#' \dontrun{
#' fit <- bayes_surv_fit(load_trial_data(), iter = 1000, chains = 2)
#' plot_posterior_hr(fit)
#' }
#'
#' @export
plot_posterior_hr <- function(fit, name = "posterior_hr.png", dir = "figures",
                              save = TRUE) {
  if (!inherits(fit, "trialsurv_bayes")) {
    rlang::abort("`fit` must be a trialsurv_bayes object.",
      class = "trialsurv_input_error"
    )
  }
  if (is.null(fit$hr_draws)) {
    rlang::abort("plot_posterior_hr() requires the default hierarchical formula.",
      class = "trialsurv_input_error"
    )
  }
  bs <- bayes_summary(fit)
  hr <- bs[bs$type == "hazard_ratio", , drop = FALSE]
  if (nrow(hr) == 0L) {
    rlang::abort("No hazard-ratio draws were found in this fit.",
      class = "trialsurv_input_error"
    )
  }

  arm_labels <- stats::setNames(fit$prep$arms$level, fit$prep$arms$dummy)
  stratum_labels <- c(
    n01 = "Nodes 0-1", n23 = "Nodes 2-3", n47 = "Nodes 4-7",
    n8p = "Nodes 8+"
  )
  parsed <- lapply(hr$parameter, function(p) {
    m <- regmatches(p, regexec("^hr_([^_]+)_(.+)$", p))[[1]]
    if (length(m) < 3L) {
      return(c(dummy = NA, stratum = NA))
    }
    c(dummy = m[2], stratum = m[3])
  })
  parsed <- do.call(rbind, parsed)
  hr$dummy <- parsed[, "dummy"]
  hr$stratum <- parsed[, "stratum"]
  hr$arm <- unname(arm_labels[hr$dummy])
  hr$is_global <- hr$stratum == "global"
  hr$stratum_label <- ifelse(hr$is_global, "Pooled (global)",
    unname(stratum_labels[hr$stratum])
  )
  hr$stratum_label[is.na(hr$stratum_label)] <- hr$stratum[is.na(hr$stratum_label)]
  hr$label <- paste0(hr$arm, ": ", hr$stratum_label)
  hr$label <- factor(hr$label, levels = rev(unique(hr$label[order(
    hr$arm, hr$is_global, hr$stratum
  )])))

  p <- ggplot2::ggplot(
    hr,
    ggplot2::aes(x = .data$mean, y = .data$label, colour = .data$is_global)
  ) +
    ggplot2::geom_vline(xintercept = 1, linetype = 2, colour = "grey40") +
    ggplot2::geom_pointrange(
      ggplot2::aes(xmin = .data$q2.5, xmax = .data$q97.5),
      linewidth = 0.6
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_colour_manual(
      values = c(`FALSE` = "#3B6EA5", `TRUE` = "#B22222"),
      labels = c(`FALSE` = "Subgroup", `TRUE` = "Global"),
      name = NULL
    ) +
    ggplot2::facet_wrap(~arm, scales = "free_y", ncol = 1) +
    ggplot2::labs(
      title = "Posterior hazard ratios with partial borrowing",
      subtitle = "Points are posterior means; bars are 95% credible intervals",
      x = "Hazard ratio (log scale)",
      y = NULL,
      caption = "The global estimate borrows across all nodal strata; subgroup estimates are shrunk toward it."
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "top")

  if (isTRUE(save)) {
    fig_save(p, name, dir = dir, width = 9, height = 7)
  }
  invisible(p)
}

#' Trace plot of key posterior quantities
#'
#' Draws trace plots for the global treatment hazard ratios and the Weibull
#' shape parameter, one chain per colour, and saves the figure with
#' [fig_save()].
#'
#' @param fit A `trialsurv_bayes` fit.
#' @param name Output file name (with or without `.png`).
#' @param dir Output directory.
#' @param save Logical; save the figure.
#'
#' @return The plot object, invisibly.
#'
#' @examples
#' \dontrun{
#' fit <- bayes_surv_fit(load_trial_data(), iter = 1000, chains = 2)
#' plot_trace(fit)
#' }
#'
#' @export
plot_trace <- function(fit, name = "trace.png", dir = "figures", save = TRUE) {
  if (!inherits(fit, "trialsurv_bayes")) {
    rlang::abort("`fit` must be a trialsurv_bayes object.",
      class = "trialsurv_input_error"
    )
  }
  draws <- fit$draws
  vars <- dimnames(draws)[[3]]
  it <- dim(draws)[1]
  ch <- dim(draws)[2]

  pieces <- list()
  if (!is.null(fit$hr_draws)) {
    hr <- posterior::as_draws_array(fit$hr_draws)
    hvars <- grep("_global$", dimnames(hr)[[3]], value = TRUE)
    for (nm in hvars) {
      pieces[[length(pieces) + 1L]] <- tibble::tibble(
        parameter = nm,
        iteration = rep(seq_len(it), times = ch),
        chain = rep(seq_len(ch), each = it),
        value = log(as.numeric(hr[, , nm]))
      )
    }
  }
  if ("shape" %in% vars) {
    pieces[[length(pieces) + 1L]] <- tibble::tibble(
      parameter = "shape",
      iteration = rep(seq_len(it), times = ch),
      chain = rep(seq_len(ch), each = it),
      value = as.numeric(draws[, , "shape"])
    )
  }
  if (length(pieces) == 0L) {
    rlang::abort("No traceable parameters were found in this fit.",
      class = "trialsurv_input_error"
    )
  }
  df <- dplyr::bind_rows(pieces)
  arm_labels <- stats::setNames(fit$prep$arms$level, fit$prep$arms$dummy)
  for (j in seq_len(nrow(fit$prep$arms))) {
    old <- paste0("hr_", fit$prep$arms$dummy[j], "_global")
    new <- paste0(arm_labels[[fit$prep$arms$dummy[j]]], " vs ",
      fit$prep$ref_arm, " (log HR)"
    )
    df$parameter[df$parameter == old] <- new
  }

  p <- ggplot2::ggplot(
    df,
    ggplot2::aes(x = .data$iteration, y = .data$value,
      colour = factor(.data$chain))
  ) +
    ggplot2::geom_line(linewidth = 0.25, alpha = 0.75) +
    ggplot2::facet_wrap(~parameter, scales = "free_y", ncol = 1) +
    ggplot2::labs(
      title = "Trace plots for key posterior quantities",
      subtitle = "Well-mixed chains overlap and show no drift",
      x = "Post-warm-up iteration",
      y = "Value",
      colour = "Chain"
    ) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "top")

  if (isTRUE(save)) {
    fig_save(p, name, dir = dir, width = 10, height = 8)
  }
  invisible(p)
}
