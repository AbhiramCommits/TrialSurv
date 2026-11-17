#' Kaplan-Meier plot with risk table
#'
#' Draws Kaplan-Meier curves by group using `ggsurvfit`, optionally with a
#' number-at-risk table underneath, and saves the figure to `figures/`.
#'
#' @param data A tidy trial data frame, usually from [load_trial_data()].
#' @param by Character scalar naming the grouping column, or `NULL`.
#' @param endpoint One of `"death"` or `"recurrence"`.
#' @param name Output file name (with or without `.png`). Defaults to
#'   `"km_by_<by>.png"`.
#' @param dir Output directory.
#' @param save Logical; save the figure with [fig_save()].
#' @param risk_table Logical; include the number-at-risk table.
#'
#' @return The plot object, invisibly.
#'
#' @examples
#' \dontrun{
#' d <- load_trial_data()
#' plot_km(d, by = "rx", name = "km_by_arm.png")
#' }
#'
#' @export
plot_km <- function(data, by = "rx", endpoint = "death", name = NULL,
                    dir = "figures", save = TRUE, risk_table = TRUE) {
  endpoint <- match.arg(endpoint, c("recurrence", "death"))
  if (!requireNamespace("ggsurvfit", quietly = TRUE)) {
    rlang::abort("Package 'ggsurvfit' is required for plot_km().",
      class = "trialsurv_input_error"
    )
  }
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.", class = "trialsurv_input_error")
  }
  d <- data[data$endpoint == endpoint, , drop = FALSE]
  assert_surv_inputs(d, vars = by)

  rhs <- if (is.null(by)) "1" else by
  form <- stats::reformulate(rhs, response = "survival::Surv(time, status)")
  fit <- ggsurvfit::survfit2(form, data = d)

  p <- ggsurvfit::ggsurvfit(fit, linewidth = 1) +
    ggsurvfit::add_confidence_interval(type = "ribbon", alpha = 0.15) +
    ggplot2::scale_y_continuous(
      limits = c(0, 1),
      labels = function(x) sprintf("%.0f%%", 100 * x)
    ) +
    ggplot2::labs(
      title = sprintf(
        "Kaplan-Meier survival by %s (%s endpoint)",
        if (is.null(by)) "cohort" else by, endpoint
      ),
      x = "Time since randomisation (months)",
      y = "Survival probability",
      colour = if (is.null(by)) NULL else by,
      fill = if (is.null(by)) NULL else by
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "top")

  if (isTRUE(risk_table)) {
    p <- p + ggsurvfit::add_risktable(risktable_stats = "{n.risk}")
  }

  if (!is.null(by)) {
    sd <- survival::survdiff(form, data = d)
    df <- length(sd$n) - 1L
    pval <- stats::pchisq(sd$chisq, df = df, lower.tail = FALSE)
    p <- p + ggplot2::labs(
      caption = sprintf("Log-rank p = %s", format.pval(pval, digits = 3, eps = 0.001))
    )
  }

  if (isTRUE(save)) {
    if (is.null(name)) {
      name <- if (is.null(by)) "km_overall.png" else paste0("km_by_", by, ".png")
    }
    fig_save(p, name, dir = dir, width = 9, height = 8)
  }
  invisible(p)
}

#' Schoenfeld residual plot for proportional hazards diagnostics
#'
#' Plots scaled Schoenfeld residuals against time for every coefficient in a
#' Cox model, with a loess trend and a horizontal reference line at zero. A
#' trend away from zero suggests a time-varying effect (a proportional hazards
#' violation). The figure is saved to `figures/`.
#'
#' @param fit A [survival::coxph()] fit or a `trialsurv_cox` object (the
#'   adjusted model is used).
#' @param name Output file name (with or without `.png`).
#' @param dir Output directory.
#' @param save Logical; save the figure with [fig_save()].
#'
#' @return The plot object, invisibly.
#'
#' @examples
#' \dontrun{
#' d <- load_trial_data()
#' cx <- cox_fit(d)
#' plot_schoenfeld(cx)
#' }
#'
#' @export
plot_schoenfeld <- function(fit, name = "schoenfeld.png", dir = "figures",
                            save = TRUE) {
  if (inherits(fit, "trialsurv_cox")) fit <- fit$adjusted
  if (!inherits(fit, "coxph")) {
    rlang::abort("`fit` must be a coxph or trialsurv_cox object.",
      class = "trialsurv_input_error"
    )
  }

  res <- as.matrix(stats::residuals(fit, type = "schoenfeld"))
  tt <- attr(res, "time")
  if (is.null(tt)) tt <- suppressWarnings(as.numeric(rownames(res)))
  tt <- as.numeric(tt)
  if (length(tt) != nrow(res)) {
    rlang::abort("Could not align Schoenfeld residuals with their times.",
      class = "trialsurv_input_error"
    )
  }

  plot_df <- tibble::tibble(
    term = rep(colnames(res), each = nrow(res)),
    time = rep(as.numeric(tt), times = ncol(res)),
    residual = as.vector(res)
  )
  plot_df <- plot_df[is.finite(plot_df$residual), , drop = FALSE]

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$time, y = .data$residual)
  ) +
    ggplot2::geom_hline(yintercept = 0, linetype = 2, colour = "grey40") +
    ggplot2::geom_point(alpha = 0.45, size = 0.9, colour = "#3B6EA5") +
    ggplot2::geom_smooth(
      method = "loess", formula = y ~ x, se = FALSE,
      colour = "#B22222", linewidth = 0.9, span = 0.8
    ) +
    ggplot2::facet_wrap(~term, scales = "free_y", ncol = 2) +
    ggplot2::labs(
      title = "Scaled Schoenfeld residuals over time",
      subtitle = "A flat loess trend supports proportional hazards",
      x = "Time (months)",
      y = "Schoenfeld residual"
    ) +
    ggplot2::theme_bw(base_size = 11)

  if (isTRUE(save)) {
    fig_save(p, name, dir = dir, width = 10, height = 9)
  }
  invisible(p)
}
