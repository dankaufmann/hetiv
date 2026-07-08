# Estimate impulse responses via proxy-IV local projections

Estimates impulse response functions (IRFs) using user-provided external
instruments (proxies) combined with local projections (Jorda, 2005). The
proxy variables serve directly as instruments for the endogenous shock
variables. Optionally imposes recursive zero restrictions across shock
dimensions and supports deterministic controls following the same
conventions as
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
  Hstep = 1,
  cov_type = "HC0",
  recursive = FALSE,
  details = FALSE
)
```

## Arguments

- y:

  Numeric matrix of stationary outcome variables (T x N). The effect on
  the first variable in each dimension is normalized to `norm` at
  horizon 0. These variables are also used as the endogenous regressors
  instrumented by the columns of `Z`.

- O:

  Numeric matrix of information set variables (T x M). May be identical
  to `y`. Included as lags 1 through `P`.

- Z:

  Numeric matrix of external instruments (T x E). Column `e` is used as
  the proxy for shock dimension `e`. Missing values on control days
  (`Ind == 0`) are treated as contaminated and excluded from estimation.
  Missing values on policy days (`Ind == 1`) remain labelled as events
  for residual-based shock prediction, but are necessarily dropped from
  IV regressions that use the missing proxy.

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

  Integer. Number of shock dimensions to identify. Default `1`.

- norm:

  Numeric scalar. Normalize the impact response of the first variable to
  this value. Set to `1` for standard unit-effect normalization.

- cum:

  Logical scalar or vector of length N. If `TRUE` for variable `i`, the
  cumulative IRF is reported. A single value is recycled to all
  variables. Default `FALSE`.

- Hstep:

  Integer. Step size between horizons. The default `1` estimates all
  horizons 0 through H - 1. Values greater than 1 estimate only the
  selected horizons. Default `1`.

- cov_type:

  Covariance estimator for local-projection standard errors: `"HC0"`
  (default) for heteroskedasticity-robust standard errors or `"NW"` for
  Newey-West HAC standard errors. `"HC0"` is the default because Montiel
  Olea et al. (2025) show that heteroskedasticity-robust standard errors
  suffice for local-projection impulse responses under weak conditions,
  even though multi-step forecast errors are typically serially
  correlated. `"NW"` remains available as an optional HAC robustness
  check.

- recursive:

  Logical. If `TRUE`, imposes recursive zero restrictions across shock
  dimensions: for shock `e > 1`, the variables and instruments from
  dimensions `1, ..., e-1` are added as controls. Default `FALSE`.

- details:

  Logical. If `TRUE`, returns detailed results including IV model
  objects, OLS residuals, and covariance matrices. If `FALSE` (default),
  returns only impulse responses and standard errors (faster; use for
  bootstrap).

## Value

A named list. Always contains:

- `irf`:

  Array (H x N x E) of estimated impulse responses.

- `se`:

  Array (H x N x E) of local-projection standard errors.

- `Method`:

  Character string `"Proxy-IV"`.

With `details = TRUE`, additionally contains:

- `IVRes`:

  List of `ivreg` model objects, one per horizon, variable, and shock
  dimension.

- `OLSRes`:

  List of OLS model objects used for residual-based covariance
  estimation, one per outcome variable.

- `Obs`:

  Data frame with observation counts: `Tp` (policy days), `Tc` (control
  days), `To` (contaminated days), `Tt` (labelled event or control
  days), and `Tiv` (complete observations available to the IV regression
  at the impact horizon).

- `et`:

  Data frame of OLS residuals on event days.

- `Sig`:

  Covariance matrix of residuals on event days.

- `SigR`:

  Covariance matrix of residuals on control days, or `NA` if
  unavailable.

- `Psi`:

  Impact matrix (N x E), equal to `irf[1, , ]`. By the package's
  indexing convention `HSeries` starts at 1, so the first LP uses
  `lead(y, 0)` (the contemporaneous value) and is labelled horizon 0;
  `irf[1, , ]` is therefore always the impact response.

- `WeakData`:

  Data frame of endogenous variables and instruments for the
  Lewis-Mertens (2025) weak instrument test.

## Details

For `E > 1`, identification can be order-dependent: the column order of
`y` defines the endogenous shock variables and normalizations, the
column order of `Z` assigns proxies to shock dimensions, and
`recursive = TRUE` imposes restrictions in that order.

## References

Jorda, O. (2005). Estimation and inference of impulse responses by local
projections. *American Economic Review*, 95(1), 161-182.

Lewis, D. J. and Mertens, K. (2025). A robust test for weak instruments
for 2SLS with multiple endogenous regressors. *The Review of Economic
Studies*, DOI: 10.1093/restud/rdaf103

Mertens, K. and Ravn, M. O. (2013). The dynamic effects of personal and
corporate income tax changes in the United States. *American Economic
Review*, 103(4), 1212-1247.

Montiel Olea, J. L., M. Plagborg-Moller, E. Qian, and C. K. Wolf (2025).
Local projections or VARs? A primer for macroeconomists. *NBER Working
Paper* No. 33871.

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
