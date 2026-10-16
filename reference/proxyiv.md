# Estimate impulse responses via proxy-IV local projections

Estimates impulse response functions (IRFs) using external instruments
(proxies) combined with local projections (Jorda, 2005, Stock and
Watson, 2018). The proxy variables serve as instruments for the
endogenous variables. Optionally imposes recursive zero restrictions
across shock dimensions. It also returns weak instrument tests as well
as shock predictions on event days, and supports deterministic controls
following the same conventions as
[`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md).

## Usage

``` r
proxyiv(
  y,
  O,
  Z,
  X = NULL,
  Ind,
  P,
  H,
  E = 1,
  norm = 1,
  cum = FALSE,
  cov_type = "HC3",
  recursive = FALSE,
  hpredict = 0,
  pooled = FALSE,
  details = FALSE
)
```

## Arguments

- y:

  Numeric matrix of outcome variables (T x N). The effect on the first
  variable in each dimension is normalized to `norm` at horizon 0. These
  variables are also used as the endogenous regressors instrumented by
  the columns of `Z`.

- O:

  Numeric matrix of information set variables (T x M). For
  non-stationary variables `y`, it should at least include `y`. Included
  as lags 1 through `P`.

- Z:

  Numeric matrix of external instruments (T x E). Column `e` is used as
  the proxy for shock dimension `e` for the corresponding endogenous
  variable. Following the proxy-SVAR literature (Mertens and Ravn,
  2013), the recommended convention is to set the proxy to zero on
  control days (`Ind == 0`). These days then remain in the IV sample and
  sharpen the estimates of the control coefficients, while the
  identifying moments are unaffected. Alternatively, missing values
  control days are treated as contaminated and excluded from estimation,
  which reproduces the event-day sample of the event-study literature
  (Nakamura and Steinsson, 2018). With zeros, always prefer the robust
  first-stage statistics (`F_rob`, `WeakTest`) over the homoskedastic F,
  which treats the zero observations as informative. Missing values on
  policy days (`Ind == 1`) remain labelled as events for residual-based
  shock prediction, but are necessarily dropped from IV regressions that
  use the missing proxy.

- X:

  Numeric matrix of deterministic variables (T x K), or `NULL`
  (default). May include a constant, time trend, or seasonal dummies.
  Included as-is (no lags).

- Ind:

  Integer vector of length T, event indicator:

  - `0` Control day (no event)

  - `1` Policy day (event)

  - `2` Contaminated control day (excluded from estimation)

- P:

  Integer. Maximum lag order for the information set. Set to `0` for no
  lags (regression on constant only).

- H:

  Integer. Maximum horizon (in periods) up to which IRFs are estimated.

- E:

  Integer. Number of shock dimensions to identify. Default `1`. Needs to
  be equal to the number of columns in `Z`

- norm:

  Numeric scalar. Normalize the impact response of the first variable to
  this value. Defaults to `1` for standard unit-effect normalization.

- cum:

  Logical scalar or vector of length N. If `TRUE` for variable `i`, the
  cumulative IRF is reported. A single value is applied to all
  variables. Default `FALSE`.

- cov_type:

  Covariance estimator for local-projection standard errors: `"HC3"`
  (default) for heteroskedasticity-robust standard errors or `"NW"` for
  Newey-West HAC standard errors. `"HC3"` is the default because Montiel
  Olea et al. (2025) show that heteroskedasticity-robust standard errors
  suffice for local-projection impulse responses under weak conditions,
  even though multi-step forecast errors are typically serially
  correlated. `"NW"` remains available as an optional HAC robustness
  check.

- recursive:

  Logical. If `TRUE`, imposes recursive zero restrictions across shock
  dimensions: for shock `e > 1`, the variables and instruments from
  dimensions `1, ..., e-1` are added as controls. This implies that the
  shock identified for dimension e has no impact effect on the
  endogenous variables from dimensions `1, ..., e-1`. Default `FALSE`.

- hpredict:

  Integer. Horizon of the responses and forecast errors used for shock
  prediction, on the same scale as the row labels of `irf`. Defaults to
  0 (recommended; the impact response and one-step-ahead residual are
  used). Must lie between 0 and `H - 1`.

- pooled:

  Logical. If `FALSE` (default), the shock prediction uses the forecast
  errors and responses at horizon `hpredict` only. If `TRUE`, the
  forecast errors of all horizons `0, ..., hpredict` are stacked and
  used jointly (pooled prediction); see
  [`kfpredict()`](https://dankaufmann.github.io/hetiv/reference/kfpredict.md).
  With `hpredict = 0` both options coincide.

- details:

  Logical. If `TRUE`, saves detailed IV results, weak instrument tests,
  and shock prediction results, which is slightly slower. If `FALSE`,
  returns only impulse responses and standard errors.

## Value

A named list. Always contains:

- `irf`:

  Array (H x N x E) of estimated impulse responses.

- `se`:

  Array (H x N x E) of local-projection standard errors.

- `Method`:

  Character string `"Proxy-IV"`.

With `details = TRUE`, additionally contains:

- `Shocks`:

  Predicted structural shocks on event days from
  [`kfpredict()`](https://dankaufmann.github.io/hetiv/reference/kfpredict.md),
  based on the responses at horizon `hpredict`.

- `IVRes`:

  List of `ivreg` model objects, one per horizon, variable, and shock
  dimension.

- `OLSRes`:

  List of the information-set regressions whose residuals enter the
  shock prediction, named `OLS.n<i>.h<h>` by outcome variable and
  horizon.

- `Obs`:

  Data frame with observation counts: `Tp` (policy days), `Tc` (control
  days), `To` (contaminated days), `Tt` (labelled event or control
  days), and `Tiv` (complete observations available to the IV regression
  of the first shock at the impact horizon).

- `FirstStage`:

  For `E = 1`: data frame with one row per horizon `h` giving the number
  of observations `Tobs`, the homoskedastic first-stage F statistic
  `F_hom` (the squared t statistic of the proxy; Stock and Yogo, 2005)
  and the robust effective F statistic `F_rob` of Montiel Olea and
  Pflueger (2013), always computed with the HAC (Newey-West) covariance
  estimator, irrespective of `cov_type`. The statistics do not depend on
  the outcome variable. `NULL` for `E > 1`.

- `WeakTest`:

  For `E > 1`: data frame with one row per shock dimension `target`
  (`e = 1, ..., E`), horizon `h` and outcome variable `i`, giving the
  results of
  [`gweakivtest()`](https://dankaufmann.github.io/hetiv/reference/gweakivtest.md)
  for the system with the first `e` policy variables as endogenous
  regressors and the first `e` proxies as instruments, targeting the
  coefficient on policy variable `e`: `Tobs`, `gmin` (generalized
  minimum eigenvalue statistic, which depends on the horizon only
  through the sample), `cv` and `cv_simplified` (sharp and simplified
  critical values, which may depend on the outcome and horizon), and the
  Stock-Yogo statistic `sy_stat` with its critical value `sy_cv` (`NaN`
  unless over-identified). The test always uses the HAC (Newey-West)
  covariance estimator, irrespective of `cov_type`, and the relative
  bias criterion (`crit = "rel"`), for comparability with Montiel Olea
  and Pflueger (2013). Rows where the outcome is a policy variable
  (`i <= E`) on impact are omitted. `NULL` for `E = 1`.

- `et`:

  Data frame of forecast errors on event days used for shock prediction,
  one column per variable and horizon (`y1.h0`, ..., `yN.h0`, `y1.h1`,
  ...): horizon `hpredict` only if `pooled = FALSE`, horizons
  `0, ..., hpredict` if `pooled = TRUE`.

- `Sig`:

  Covariance matrix of `et` on event days.

- `SigR`:

  Covariance matrix of the same residuals on control days, or `NA` if
  unavailable. A control day enters horizon h only if no event occurs in
  its window `t + 1, ..., t + h`.

- `Psi`:

  Response matrix used for shock prediction, with rows matching the
  columns of `et`: `irf[hpredict + 1, , ]` (N x E) if `pooled = FALSE`,
  the rows `1, ..., hpredict + 1` of `irf` stacked if `pooled = TRUE`.
  With the default `hpredict = 0` this is the impact response.

- `WeakData`:

  Data frame of endogenous variables and instruments for the
  Lewis-Mertens (2025) weak instrument test.

## Details

For `E > 1`, identification can be order-dependent: the column order of
`y` defines the endogenous shock variables and normalizations, the
column order of `Z` assigns proxies to shock dimensions, and
`recursive = TRUE` imposes zero restrictions in that order. If
`recursive = FALSE` the instruments are assumed to be valid across
dimensions.

## References

Jorda, O. (2005). Estimation and inference of impulse responses by local
projections. *American Economic Review*, 95(1), 161-182.

Lewis, D. J. and Mertens, K. (2025). A robust test for weak instruments
for 2SLS with multiple endogenous regressors. *Review of Economic
Studies*, DOI: 10.1093/restud/rdaf103

Mertens, K. and Ravn, M. O. (2013). The dynamic effects of personal and
corporate income tax changes in the United States. *American Economic
Review*, 103(4), 1212-1247.

Montiel Olea, J. L., M. Plagborg-Moller, E. Qian, and C. K. Wolf (2025).
Local projections or VARs? A primer for macroeconomists. *NBER Working
Paper* No. 33871.

Nakamura, E. and Steinsson, J. (2018). High-frequency identification of
monetary non-neutrality: The information effect. *Quarterly Journal of
Economics*, 133(3), 1283-1330.

Stock, J. H. and Watson, M. W. (2018). Identification and estimation of
dynamic causal effects in macroeconomics using external instruments.
*Economic Journal*, 128(610), 917-948.

## Examples

``` r
set.seed(1)
y <- matrix(rnorm(80), ncol = 2)
Ind <- rep(0L, nrow(y))
Ind[seq(5, nrow(y), by = 5)] <- 1L
Z <- matrix(Ind * y[, 1] + rnorm(nrow(y)), ncol = 1)
res <- proxyiv(
  y = y, O = y, Z = Z, Ind = Ind, P = 1, H = 3,
  details = TRUE
)
dim(res$irf)
#> [1] 3 2 1
```
