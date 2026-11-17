# trialsurv 0.1.0

* Initial release.
* Added `load_trial_data()` for tidy, complete-case loading of the
  `survival::colon` adjuvant chemotherapy trial with both the recurrence and
  death endpoints and recorded cohort/missingness metadata.
* Added Kaplan-Meier estimation (`km_fit()`), Cox proportional hazards models
  (`cox_fit()`), proportional hazards diagnostics (`cox_diagnostics()`),
  held-out concordance (`concordance_eval()`) and the corresponding
  `plot_km()` and `plot_schoenfeld()` figures.
* Added a hierarchical Bayesian Weibull survival model with partial borrowing
  across nodal-risk strata (`bayes_surv_fit()`), posterior summaries
  (`bayes_summary()`), shrinkage quantification (`shrinkage_report()`),
  posterior concordance (`bayes_concordance()`) and a frequentist-versus-
  Bayesian comparison (`compare_models()`). The brms/rstan backend is the
  primary implementation; a self-contained Metropolis-within-Gibbs fallback is
  available with `backend = "mcmc"`.
* Added a two-arm trial simulation engine (`simulate_trial()`), Monte Carlo
  power and type I error estimation (`trial_power()`) and a simulated versus
  Schoenfeld sample-size search (`required_n()`).
* Added a reproducible, stratified `trial_split()` used by both the
  frequentist and Bayesian analyses.
* Added `render_report()` and a Quarto GO/NO-GO report that reads only the
  generated `results/` and `figures/` artifacts.
* Added testthat coverage, lintr/styler checks, GitHub Actions CI and a
  workflow-dispatch report workflow.
