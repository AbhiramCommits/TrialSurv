# Contributing to trialsurv

## Git workflow

1. Branch from `main` for every change:
   `git checkout -b feat/short-description` or `fix/short-description`.
2. Use conventional commit messages:
   `feat:`, `fix:`, `test:`, `docs:`, `chore:`, `refactor:`.
3. Keep commits focused. Analysis changes should include the regenerated
   `results/` and `figures/` artifacts if they change reported numbers.

## Before opening a pull request

Run all of the following from the package root:

```r
devtools::document()
devtools::test()
devtools::check()
lintr::lint_package()
styler::style_pkg()
```

The project targets 0 errors, 0 warnings and 0 notes from `devtools::check()`
in the supported local configuration. The optional Stan smoke test in
`tests/testthat/test-bayes.R` is skipped on CI; run it locally with
`TRIALSURV_TEST_BRMS=true` if you have a working rstan installation.

## Analysis changes

- Keep every reported number traceable to a file in `results/`. Do not edit a
  result CSV by hand.
- Run `Rscript scripts/run_all.R` to regenerate the full pipeline and
  `results/session_info.txt`.
- Re-render the report with `Rscript -e 'trialsurv::render_report()'` and
  commit `report/trial_report.html`.
- Update `NEWS.md` for user-visible changes.

## Code style

- The code is formatted with `styler` and checked with `lintr` (configuration
  in `.lintr`).
- New exported functions need roxygen2 documentation and at least one
  testthat test. Prefer deterministic tests with explicit seeds.
