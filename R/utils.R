#' Assert that a data frame is ready for survival analysis
#'
#' Internal validation helper used by the modelling functions. It checks that
#' the required columns exist, that follow-up time is strictly positive and
#' finite, and that the event indicator is coded only as 0/1.
#'
#' @param data A data frame.
#' @param time Character scalar naming the follow-up time column.
#' @param status Character scalar naming the event indicator column.
#' @param vars Optional character vector of additional columns that must be
#'   present.
#'
#' @return `TRUE`, invisibly. Errors with class `trialsurv_input_error` when a
#'   check fails.
#' @keywords internal
assert_surv_inputs <- function(data, time = "time", status = "status", vars = NULL) {
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.", class = "trialsurv_input_error")
  }
  needed <- unique(c(time, status, vars))
  absent <- setdiff(needed, names(data))
  if (length(absent) > 0L) {
    rlang::abort(
      sprintf("Missing required column(s): %s.", paste(absent, collapse = ", ")),
      class = "trialsurv_input_error"
    )
  }
  tvec <- data[[time]]
  svec <- data[[status]]
  if (!is.numeric(tvec)) {
    rlang::abort(sprintf("`%s` must be numeric.", time), class = "trialsurv_input_error")
  }
  if (anyNA(tvec) || any(!is.finite(tvec))) {
    rlang::abort(sprintf("`%s` contains missing or non-finite values.", time),
      class = "trialsurv_input_error"
    )
  }
  if (any(tvec <= 0)) {
    rlang::abort(sprintf("`%s` must be strictly positive.", time),
      class = "trialsurv_input_error"
    )
  }
  if (anyNA(svec) || !all(svec %in% c(0, 1))) {
    rlang::abort(sprintf("`%s` must be coded 0/1 with no missing values.", status),
      class = "trialsurv_input_error"
    )
  }
  invisible(TRUE)
}

#' Write an analysis result to the results directory
#'
#' Writes a data frame or a named list to `results/<name>.csv` or
#' `results/<name>.json`, creating the directory if needed. Data frames default
#' to CSV; named lists default to JSON. The object is returned invisibly so
#' that scripts can keep using it.
#'
#' @param x A data frame, tibble, or named list to write.
#' @param name File stem, with or without a `.csv`/`.json` extension.
#' @param dir Output directory, relative to the working directory.
#' @param format One of `"auto"`, `"csv"` or `"json"`.
#' @param digits Numeric digits passed to [jsonlite::write_json()] for JSON
#'   output. `NA` keeps full precision.
#'
#' @return `x`, invisibly. The written path is attached as `attr(x, "path")`.
#' @keywords internal
write_result <- function(x, name, dir = "results",
                         format = c("auto", "csv", "json"), digits = NA) {
  format <- match.arg(format)
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    rlang::abort("`name` must be a non-empty character string.",
      class = "trialsurv_input_error"
    )
  }
  stem <- sub("\\.(csv|json)$", "", name, ignore.case = TRUE)
  supplied_ext <- tolower(sub("^.*\\.(csv|json)$", "\\1", name, ignore.case = TRUE))
  if (grepl("\\.(csv|json)$", name, ignore.case = TRUE) && identical(format, "auto")) {
    format <- supplied_ext
  }
  if (identical(format, "auto")) {
    format <- if (is.data.frame(x)) "csv" else "json"
  }

  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  path <- file.path(dir, paste0(stem, ".", format))

  if (identical(format, "csv")) {
    if (!is.data.frame(x)) {
      if (is.list(x) && all(vapply(x, length, integer(1)) == 1L)) {
        x <- tibble::as_tibble(x)
      } else {
        rlang::abort("CSV output requires a data frame or a list of scalars.",
          class = "trialsurv_input_error"
        )
      }
    }
    utils::write.csv(as.data.frame(x), path, row.names = FALSE, na = "")
  } else {
    if (is.data.frame(x)) {
      x <- as.list(x)
    }
    jsonlite::write_json(x, path, auto_unbox = TRUE, pretty = TRUE, digits = digits)
  }

  message("Wrote ", path)
  attr(x, "path") <- path
  invisible(x)
}

#' Save a ggplot figure to the figures directory
#'
#' A thin wrapper around [ggplot2::ggsave()] that standardises output location,
#' resolution and background colour.
#'
#' @param plot A ggplot object (or a plot object accepted by `ggsave`).
#' @param name File name, with or without a `.png` extension.
#' @param dir Output directory.
#' @param width,height Plot size in inches.
#' @param dpi Resolution in dots per inch.
#'
#' @return The output path, invisibly.
#' @keywords internal
fig_save <- function(plot, name, dir = "figures", width = 8, height = 6, dpi = 300) {
  if (!is.character(name) || length(name) != 1L || !nzchar(name)) {
    rlang::abort("`name` must be a non-empty character string.",
      class = "trialsurv_input_error"
    )
  }
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  }
  if (!grepl("\\.(png|pdf|svg|jpg|jpeg|tif|tiff)$", name, ignore.case = TRUE)) {
    name <- paste0(name, ".png")
  }
  path <- file.path(dir, name)
  ggplot2::ggsave(
    filename = path, plot = plot, width = width, height = height,
    dpi = dpi, bg = "white"
  )
  message("Wrote ", path)
  invisible(path)
}

# Column names used inside tidy-evaluation data masks.
utils::globalVariables(c(
  "time", ".data", "arm", "parameter", "value", "iteration", "subgroup"
))
