# Changelog

## hetiv 0.1.4

- **Breaking change:**
  [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md),
  [`proxyiv()`](https://dankaufmann.github.io/hetiv/reference/proxyiv.md)
  and
  [`recols()`](https://dankaufmann.github.io/hetiv/reference/recols.md)
  with `details = TRUE` no longer return the fitted regression objects
  (`IVRes`, `OLSRes`, `ORTHRes`) by default, because they made saved
  results very large (hundreds of megabytes). Use the new argument
  `models = TRUE` to obtain them. Stored models no longer keep the
  environment of the estimator.
- The local projections in
  [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md),
  [`proxyiv()`](https://dankaufmann.github.io/hetiv/reference/proxyiv.md)
  and
  [`recols()`](https://dankaufmann.github.io/hetiv/reference/recols.md)
  are now estimated by matrix algebra instead of
  [`ivreg::ivreg()`](https://zeileis.github.io/ivreg/reference/ivreg.html)/[`lm()`](https://rdrr.io/r/stats/lm.html)
  with `sandwich`, unless `models = TRUE`. Impulse responses and
  standard errors (HC3 and Newey-West) are unchanged; estimation is
  several times faster.
- [`kfpredict()`](https://dankaufmann.github.io/hetiv/reference/kfpredict.md)
  has been renamed to
  [`predict_shocks()`](https://dankaufmann.github.io/hetiv/reference/predict_shocks.md).
  [`kfpredict()`](https://dankaufmann.github.io/hetiv/reference/kfpredict.md)
  is kept as a deprecated wrapper and issues a warning.
- [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md),
  [`proxyiv()`](https://dankaufmann.github.io/hetiv/reference/proxyiv.md)
  and
  [`recols()`](https://dankaufmann.github.io/hetiv/reference/recols.md)
  with `details = TRUE` additionally return the coefficients of the
  one-step (horizon 0) regressions on the information set: `Phi` (N x M
  x P, lags of `O`), `Theta` (N x K, exogenous variables `X`) and
  `const` (intercepts). With `hetiv(interact = TRUE)`, `Phi` and `const`
  refer to control days and `PhiE` and `constE` to event days.
- New function
  [`prediction_corr()`](https://dankaufmann.github.io/hetiv/reference/prediction_corr.md):
  squared correlation between the predicted and the true shocks as a
  measure of the precision of the prediction. It warns if the covariance
  on control days is (nearly) singular. The estimators return it as
  `PredR2` with `details = TRUE`.
- New functions:
  [`plotseries()`](https://dankaufmann.github.io/hetiv/reference/plotseries.md)
  (line plot of several series with a single legend, optional
  standardisation and correlation label; lines connect across missing
  values unless `connect = FALSE`),
  [`gss_pca()`](https://dankaufmann.github.io/hetiv/reference/gss_pca.md)
  (target and path factors of Gurkaynak, Sack and Swanson, 2005, with EM
  principal components for missing values),
  [`lag_opt()`](https://dankaufmann.github.io/hetiv/reference/lag_opt.md)
  (lag length of a VAR by information criterion; requires the suggested
  package `vars`) and
  [`rescale_to_match()`](https://dankaufmann.github.io/hetiv/reference/rescale_to_match.md).
- `plotirf(legend = TRUE)` now adds the legend only once per shock (to
  the first panel) and
  [`arrangeirf()`](https://dankaufmann.github.io/hetiv/reference/arrangeirf.md)
  draws it once below the arranged figure. The legend shows the lines
  only (no separate key for the bands).

## hetiv 0.1.3

- [`plotirf()`](https://dankaufmann.github.io/hetiv/reference/plotirf.md)
  now plots up to three sets of impulse responses in the same panels.
  Each set is supplied as an array or as an estimator result (a list
  with `irf` and optionally `se` and `boot`) and is drawn without bands,
  with normal bands from standard errors, or with percentile bands
  computed from bootstrap draws and recentred on the point estimate
  (`boot_type = "recentred"`; `"percentile"` gives raw percentile
  bands). The default is a single 90% band (`ci = 0.90`); up to two
  confidence levels can be supplied. New arguments `names`, `colors`,
  `boot_type` and `legend`.
  [`plot2irf()`](https://dankaufmann.github.io/hetiv/reference/plot2irf.md)
  is now a thin wrapper around
  [`plotirf()`](https://dankaufmann.github.io/hetiv/reference/plotirf.md)
  kept for backward compatibility. Fixed: two-dimensional `IRFest` with
  `IRFse` no longer errors.
- [`svariv()`](https://dankaufmann.github.io/hetiv/reference/svariv.md)
  returns the bootstrap draws (`boot`) also when `details = FALSE`.
- [`svariv()`](https://dankaufmann.github.io/hetiv/reference/svariv.md)
  now computes the reduced-form residual covariance on the full VAR
  sample, as in Mertens and Ravn (2013); only the 2SLS step uses the
  observations with non-missing proxies. With `details = TRUE` it
  reports, for each policy residual, the homoskedastic first-stage F
  statistic (`FStat`) and the robust effective F statistic of Montiel
  Olea and Pflueger
  2013. (`FStatRobust`, new argument `cov_type`). Identification
        failures now raise an informative error.
- [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md)
  and
  [`proxyiv()`](https://dankaufmann.github.io/hetiv/reference/proxyiv.md)
  with `details = TRUE` now report weak-instrument diagnostics on the
  sample of every local projection: for `E = 1`, a data frame
  `FirstStage` with the homoskedastic first-stage F statistic (squared t
  statistic of the instrument) and the robust effective F statistic of
  Montiel Olea and Pflueger (2013) by horizon; for `E > 1`, a data frame
  `WeakTest` with the generalized Lewis and Mertens (2025) test from
  [`gweakivtest()`](https://dankaufmann.github.io/hetiv/reference/gweakivtest.md)
  by shock dimension (`target`), horizon and outcome variable. Both
  diagnostics always use the HAC (Newey-West) covariance estimator,
  irrespective of `cov_type`. `WeakTest` reports, for each shock
  dimension `e`, the test for the system with the first `e` policy
  variables and instruments, targeting the coefficient on policy
  variable `e`, with the relative bias criterion.
  [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md)
  now also returns the OLS residual models in `OLSRes`, and
  [`proxyiv()`](https://dankaufmann.github.io/hetiv/reference/proxyiv.md)
  computes `Tiv` on impact as documented.
- Removed the `Hstep` argument from
  [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md),
  [`proxyiv()`](https://dankaufmann.github.io/hetiv/reference/proxyiv.md)
  and
  [`recols()`](https://dankaufmann.github.io/hetiv/reference/recols.md).
  All horizons `0, ..., H - 1` are always estimated. Previously,
  `Hstep > 1` combined with `hpredict > 1` selected the wrong horizon
  for the shock prediction.
- New argument `pooled` in
  [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md),
  [`proxyiv()`](https://dankaufmann.github.io/hetiv/reference/proxyiv.md)
  and
  [`recols()`](https://dankaufmann.github.io/hetiv/reference/recols.md):
  with `pooled = TRUE` the shock prediction stacks the forecast errors
  of all horizons `0, ..., hpredict` and uses them jointly (pooled
  prediction), instead of the single horizon `hpredict`.
  [`kfpredict()`](https://dankaufmann.github.io/hetiv/reference/kfpredict.md)
  accepts `Psi` and `et` as arrays with one slice per horizon and stacks
  them, documents the pooled prediction, and warns when the residual
  covariance is singular. Control days now enter the horizon-h residual
  covariance only if no event falls in their window `t + 1, ..., t + h`.
  The columns of `et` and the elements of `OLSRes` are now named by
  variable and horizon (`y1.h0`, `OLS.n1.h0`, …).
- Breaking: `hpredict` in
  [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md),
  [`proxyiv()`](https://dankaufmann.github.io/hetiv/reference/proxyiv.md)
  and
  [`recols()`](https://dankaufmann.github.io/hetiv/reference/recols.md)
  is now a horizon on the scale of the row labels of `irf`
  (`0, ..., H - 1`, default `0` = impact) instead of a 1-based row
  index. Callers who passed a non-default value must subtract one. `Psi`
  equals `irf[hpredict + 1, , ]`.

## hetiv 0.1.2

- New
  [`svariv()`](https://dankaufmann.github.io/hetiv/reference/svariv.md):
  proxy SVAR (SVAR-IV) following Mertens and Ravn (2013). Estimates a
  VAR(P) by OLS, identifies the first `E` shocks with external
  instruments, computes IRFs by recursion, and optionally returns
  wild-bootstrap percentile bands and standard deviations that work with
  [`plotirf()`](https://dankaufmann.github.io/hetiv/reference/plotirf.md)
  and
  [`plot2irf()`](https://dankaufmann.github.io/hetiv/reference/plot2irf.md).
- New
  [`recols()`](https://dankaufmann.github.io/hetiv/reference/recols.md):
  recursive (Cholesky) local projections estimated by OLS, with the same
  interface and output as
  [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md),
  as a benchmark identification.
- New
  [`arrangeirf()`](https://dankaufmann.github.io/hetiv/reference/arrangeirf.md):
  arranges IRF panels from
  [`plotirf()`](https://dankaufmann.github.io/hetiv/reference/plotirf.md)/[`plot2irf()`](https://dankaufmann.github.io/hetiv/reference/plot2irf.md)
  into one figure per shock and optionally saves them as PDF. Adds
  `gridExtra` to Imports.
- Various smaller improvements and additional options in
  [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md)
  and
  [`proxyiv()`](https://dankaufmann.github.io/hetiv/reference/proxyiv.md)

## hetiv 0.1.1

- Changes to default HC settings for local projections

## hetiv 0.1.0

- Initial CRAN-ready release.
- Provides heteroskedasticity-based and proxy-based instrumental
  variable local projection estimators for event-study settings.
- Includes impulse response plotting helpers, weak instrument testing,
  Kalman-filter shock extraction, and simulation utilities.
- Adds an introductory vignette with a worked end-to-end example.
