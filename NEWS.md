# hetiv 0.1.3

* `plotirf()` now plots up to three sets of impulse responses in the same
  panels. Each set is supplied as an array or as an estimator result (a list
  with `irf` and optionally `se` and `boot`) and is drawn without bands, with
  normal bands from standard errors, or with percentile bands computed from
  bootstrap draws and recentred on the point estimate (`boot_type =
  "recentred"`; `"percentile"` gives raw percentile bands). The default is a
  single 90% band (`ci = 0.90`); up to two confidence levels can be supplied.
  New arguments `names`, `colors`, `boot_type` and `legend`.
  `plot2irf()` is now a thin wrapper around `plotirf()` kept for backward
  compatibility. Fixed: two-dimensional `IRFest` with `IRFse` no longer
  errors.
* `svariv()` returns the bootstrap draws (`boot`) also when `details = FALSE`.
* `svariv()` now computes the reduced-form residual covariance on the full
  VAR sample, as in Mertens and Ravn (2013); only the 2SLS step uses the
  observations with non-missing proxies. With `details = TRUE` it reports,
  for each policy residual, the homoskedastic first-stage F statistic
  (`FStat`) and the robust effective F statistic of Montiel Olea and Pflueger
  (2013) (`FStatRobust`, new argument `cov_type`). Identification failures
  now raise an informative error.
* `hetiv()` and `proxyiv()` with `details = TRUE` now report weak-instrument
  diagnostics on the sample of every local projection: for `E = 1`, a data
  frame `FirstStage` with the homoskedastic first-stage F statistic (squared
  t statistic of the instrument) and the robust effective F statistic of
  Montiel Olea and Pflueger (2013) by horizon; for `E > 1`, a data frame
  `WeakTest` with the generalized Lewis and Mertens (2025) test from
  `gweakivtest()` by horizon and outcome variable. `hetiv()` now also returns
  the OLS residual models in `OLSRes`, and `proxyiv()` computes `Tiv` on
  impact as documented.
* Removed the `Hstep` argument from `hetiv()`, `proxyiv()` and `recols()`.
  All horizons `0, ..., H - 1` are always estimated. Previously, `Hstep > 1`
  combined with `hpredict > 1` selected the wrong horizon for the shock
  prediction.
* Breaking: `hpredict` in `hetiv()`, `proxyiv()` and `recols()` is now a
  horizon on the scale of the row labels of `irf` (`0, ..., H - 1`, default
  `0` = impact) instead of a 1-based row index. Callers who passed a
  non-default value must subtract one. `Psi` equals `irf[hpredict + 1, , ]`.

# hetiv 0.1.2

* New `svariv()`: proxy SVAR (SVAR-IV) following Mertens and Ravn (2013).
  Estimates a VAR(P) by OLS, identifies the first `E` shocks with external
  instruments, computes IRFs by recursion, and optionally returns wild-bootstrap
  percentile bands and standard deviations that work with `plotirf()` and
  `plot2irf()`.
* New `recols()`: recursive (Cholesky) local projections estimated by OLS, with
  the same interface and output as `hetiv()`, as a benchmark identification.
* New `arrangeirf()`: arranges IRF panels from `plotirf()`/`plot2irf()` into one
  figure per shock and optionally saves them as PDF. Adds `gridExtra` to Imports.
* Various smaller improvements and additional options in `hetiv()` and `proxyiv()`

# hetiv 0.1.1

* Changes to default HC settings for local projections

# hetiv 0.1.0

* Initial CRAN-ready release.
* Provides heteroskedasticity-based and proxy-based instrumental variable
  local projection estimators for event-study settings.
* Includes impulse response plotting helpers, weak instrument testing,
  Kalman-filter shock extraction, and simulation utilities.
* Adds an introductory vignette with a worked end-to-end example.
