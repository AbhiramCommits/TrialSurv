#' Render the TrialSurv go/no-go Quarto report
#'
#' Renders `report/trial_report.qmd` to a self-contained HTML file and checks
#' that the output exists. The report reads the artefacts in `results/` and
#' `figures/`; it never refits a model. All required inputs are validated
#' before Quarto is called so that a missing artefact produces an actionable
#' error instead of a partially rendered report.
#'
#' @param output_dir Directory containing the Quarto template and receiving
#'   the rendered HTML.
#' @param input Path to the `.qmd` file. Defaults to
#'   `file.path(output_dir, "trial_report.qmd")`.
#' @param quiet Passed to [quarto::quarto_render()].
#' @param results_dir Directory containing the result files checked before
#'   rendering.
#' @param figures_dir Directory containing the figures checked before
#'   rendering.
#'
#' @return The path to the rendered HTML file, invisibly.
#'
#' @examples
#' \dontrun{
#' render_report()
#' }
#'
#' @export
render_report <- function(output_dir = "report", input = NULL, quiet = TRUE,
                          results_dir = "results", figures_dir = "figures") {
  if (!requireNamespace("quarto", quietly = TRUE)) {
    rlang::abort(
      "Package 'quarto' is required to render the report.",
      class = "trialsurv_report_error"
    )
  }
  if (is.null(input)) input <- file.path(output_dir, "trial_report.qmd")
  if (!file.exists(input)) {
    rlang::abort(
      sprintf("Report template not found: %s", input),
      class = "trialsurv_report_error"
    )
  }
  .check_report_inputs(results_dir, figures_dir)
  html <- file.path(output_dir, sub("\\.qmd$", ".html", basename(input)))
  quarto::quarto_render(
    input = input,
    output_file = basename(html),
    quiet = quiet
  )
  if (!file.exists(html)) {
    rlang::abort(
      sprintf("Quarto did not produce the expected HTML file: %s", html),
      class = "trialsurv_report_error"
    )
  }
  message("Rendered ", html)
  invisible(normalizePath(html))
}

#' Check that all report inputs exist
#'
#' @param results_dir Directory containing the result CSV/JSON files.
#' @param figures_dir Directory containing the figures referenced by the
#'   report.
#'
#' @return `TRUE`, invisibly. Errors with class `trialsurv_report_error`
#'   listing every missing file.
#' @keywords internal
.check_report_inputs <- function(results_dir = "results",
                                 figures_dir = "figures") {
  required_results <- c(
    "cohort_summary.csv", "km_summary.csv", "cox_hazard_ratios.csv",
    "ph_diagnostics.csv", "concordance_frequentist.csv",
    "bayes_posterior_summary.csv", "shrinkage_report.csv",
    "model_comparison.csv", "bayes_convergence.csv", "power_curve.csv",
    "type1_error.csv", "sample_size_vs_schoenfeld.csv", "simulation_meta.csv"
  )
  required_figures <- c(
    "km_by_arm.png", "schoenfeld.png", "posterior_hr.png", "power_curves.png"
  )
  missing <- c(
    file.path(results_dir, required_results)[
      !file.exists(file.path(results_dir, required_results))
    ],
    file.path(figures_dir, required_figures)[
      !file.exists(file.path(figures_dir, required_figures))
    ]
  )
  if (length(missing) > 0L) {
    rlang::abort(
      paste0(
        "Missing required report input(s):\n  ",
        paste(missing, collapse = "\n  "),
        "\nRun scripts/run_all.R (or the stage scripts) first."
      ),
      class = "trialsurv_report_error"
    )
  }
  invisible(TRUE)
}
