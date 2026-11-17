#' Load and tidy the colon cancer trial data
#'
#' Loads the `survival::colon` dataset and reshapes it to one row per patient
#' per endpoint, giving a tidy tibble suitable for Kaplan-Meier and Cox
#' modelling. Two endpoints are returned: `"recurrence"` (time to recurrence,
#' with death without recurrence treated as censoring) and `"death"` (time to
#' death from any cause).
#'
#' @param dataset Character. The dataset to load. Currently only `"colon"`,
#'   the adjuvant chemotherapy trial for colon cancer (Moertel et al., 1990),
#'   is supported.
#' @param na_action Character. `"drop"` (default) removes patients with a
#'   missing value in any analysis variable and records the excluded counts;
#'   `"keep"` retains them and flags them in the `missing` attribute.
#'
#' @return A tibble with columns `id`, `time`, `status`, `endpoint`, `rx`,
#'   `sex`, `age`, `obstruct`, `perfor`, `adhere`, `nodes`, `differ`,
#'   `extent`, `surg` and `node4`. Survival times and `surg` are expressed in
#'   months (days / 30.4375). `rx` is a factor with levels `"Obs"`, `"Lev"`
#'   and `"Lev+5FU"`; `sex` is `"Female"`/`"Male"`; `differ` and `extent` are
#'   ordered factors. The attribute `"cohort"` is a tibble with the number of
#'   patients, number of events and reverse Kaplan-Meier median follow-up for
#'   each endpoint; the attribute `"missing"` records the number of missing
#'   values per variable before any exclusion.
#'
#' @references
#' Moertel CG, Fleming TR, Macdonald JS, et al. (1990). Levamisole and
#' fluorouracil for adjuvant therapy of resected colon carcinoma.
#' *New England Journal of Medicine*, 322(6), 352-358.
#'
#' @examples
#' d <- load_trial_data()
#' attr(d, "cohort")
#' subset(d, endpoint == "death" & rx == "Obs")
#'
#' @export
load_trial_data <- function(dataset = "colon", na_action = c("drop", "keep")) {
  if (!is.character(dataset) || length(dataset) != 1L) {
    rlang::abort("`dataset` must be a single character string.",
      class = "trialsurv_input_error"
    )
  }
  if (!identical(dataset, "colon")) {
    rlang::abort(
      sprintf("Unknown dataset %s; only \"colon\" is available.", encodeString(dataset, quote = "\"")),
      class = "trialsurv_input_error"
    )
  }
  na_action <- match.arg(na_action)

  raw <- survival::colon
  raw <- as.data.frame(raw)

  long <- tibble::tibble(
    id = as.integer(raw$id),
    time = as.numeric(raw$time) / 30.4375,
    status = as.integer(raw$status),
    endpoint = factor(
      ifelse(raw$etype == 1L, "recurrence", "death"),
      levels = c("recurrence", "death")
    ),
    rx = factor(as.character(raw$rx), levels = c("Obs", "Lev", "Lev+5FU")),
    sex = factor(raw$sex, levels = c(0, 1), labels = c("Female", "Male")),
    age = as.numeric(raw$age),
    obstruct = as.integer(raw$obstruct),
    perfor = as.integer(raw$perfor),
    adhere = as.integer(raw$adhere),
    nodes = as.numeric(raw$nodes),
    differ = factor(raw$differ, levels = 1:3, labels = c("Well", "Moderate", "Poor")),
    extent = factor(raw$extent, levels = 1:4, labels = c("Submucosa", "Muscle", "Serosa", "Contiguous")),
    surg = as.numeric(raw$surg) / 30.4375,
    node4 = as.integer(raw$node4)
  )

  analysis_vars <- c(
    "time", "status", "rx", "sex", "age", "obstruct", "perfor",
    "adhere", "nodes", "differ", "extent", "surg", "node4"
  )

  missing_counts <- vapply(long[analysis_vars], function(x) sum(is.na(x)), integer(1))
  n_before <- nrow(long)
  n_patients_before <- length(unique(long$id))

  if (identical(na_action, "drop")) {
    keep <- stats::complete.cases(long[, analysis_vars, drop = FALSE])
    long <- long[keep, , drop = FALSE]
    long <- tibble::as_tibble(long)
  }

  long <- long[order(long$endpoint, long$id), , drop = FALSE]
  rownames(long) <- NULL

  cohort <- purrr::map_dfr(levels(long$endpoint), function(ep) {
    d <- long[long$endpoint == ep, , drop = FALSE]
    rev_km <- survival::survfit(survival::Surv(time, 1 - status) ~ 1, data = d)
    med <- unname(summary(rev_km)$table["median"])
    tibble::tibble(
      endpoint = ep,
      n = nrow(d),
      n_events = sum(d$status == 1L),
      median_followup = as.numeric(med)
    )
  })

  attr(long, "cohort") <- cohort
  missing_tbl <- tibble::tibble(
    variable = names(missing_counts),
    n_missing = as.integer(missing_counts),
    n_rows_before = n_before,
    n_patients_before = n_patients_before,
    n_rows_after = nrow(long),
    n_patients_after = length(unique(long$id))
  )
  attr(missing_tbl, "na_action") <- na_action
  attr(long, "missing") <- missing_tbl
  attr(long, "dataset") <- dataset
  attr(long, "time_unit") <- "months"
  long
}
