#' Reproducible stratified train/test split of trial patients
#'
#' Splits patients, rather than rows, into a training and a test set. The split
#' is stratified by treatment arm and by event status on one endpoint so that
#' both sets contain a comparable mix of arms and events. Because patients are
#' split (not rows), both the recurrence and death endpoints stay together for
#' each patient.
#'
#' @param data A tidy trial data frame, usually from [load_trial_data()].
#' @param prop Proportion of patients assigned to the training set.
#' @param seed Integer seed for the split.
#' @param endpoint Endpoint whose event status defines the strata.
#'
#' @return A list of class `trialsurv_split` with elements `train`, `test`,
#'   `train_ids`, `test_ids`, `seed`, `prop`, `endpoint` and a tidy `summary`
#'   tibble with the number of patients and events in each set.
#'
#' @examples
#' d <- load_trial_data()
#' sp <- trial_split(d, prop = 0.7, seed = 20261117)
#' sp$summary
#'
#' @export
trial_split <- function(data, prop = 0.7, seed = 20261117, endpoint = "death") {
  endpoint <- match.arg(endpoint, c("recurrence", "death"))
  if (!is.data.frame(data)) {
    rlang::abort("`data` must be a data frame.", class = "trialsurv_input_error")
  }
  if (!is.numeric(prop) || length(prop) != 1L || prop <= 0 || prop >= 1) {
    rlang::abort("`prop` must be a single number in (0, 1).",
      class = "trialsurv_input_error"
    )
  }
  if (!is.numeric(seed) || length(seed) != 1L) {
    rlang::abort("`seed` must be a single number.", class = "trialsurv_input_error")
  }
  needed <- c("id", "time", "status", "endpoint", "rx")
  absent <- setdiff(needed, names(data))
  if (length(absent) > 0L) {
    rlang::abort(
      sprintf("Missing required column(s): %s.", paste(absent, collapse = ", ")),
      class = "trialsurv_input_error"
    )
  }

  d <- data[data$endpoint == endpoint, , drop = FALSE]
  if (nrow(d) == 0L) {
    rlang::abort(sprintf("No rows for endpoint %s.", encodeString(endpoint, quote = "\"")),
      class = "trialsurv_input_error"
    )
  }
  if (anyDuplicated(d$id)) {
    rlang::abort("Expected one row per patient for the stratification endpoint.",
      class = "trialsurv_input_error"
    )
  }

  set.seed(seed)
  key <- interaction(d$rx, d$status, drop = TRUE)
  train_rows <- unlist(lapply(split(seq_len(nrow(d)), key), function(i) {
    sample(i, size = floor(prop * length(i)))
  }), use.names = FALSE)
  train_ids <- sort(unique(d$id[train_rows]))
  test_ids <- sort(setdiff(unique(data$id), train_ids))

  train <- data[data$id %in% train_ids, , drop = FALSE]
  test <- data[data$id %in% test_ids, , drop = FALSE]

  make_summary <- function(x, set) {
    dx <- x[x$endpoint == endpoint, , drop = FALSE]
    tibble::tibble(
      set = set,
      n_patients = length(unique(x$id)),
      n_events = sum(dx$status == 1L),
      event_rate = mean(dx$status == 1L)
    )
  }
  summary_tbl <- dplyr::bind_rows(
    make_summary(train, "train"),
    make_summary(test, "test")
  )

  out <- list(
    train = train,
    test = test,
    train_ids = train_ids,
    test_ids = test_ids,
    seed = as.integer(seed),
    prop = prop,
    endpoint = endpoint,
    summary = summary_tbl
  )
  class(out) <- c("trialsurv_split", "list")
  out
}
