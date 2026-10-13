# Estimate impulse responses via heteroskedasticity-based IV local projections

Estimates impulse response functions (IRFs) using recursive
heteroskedasticity-IV identification (Rigobon, 2003; Rigobon and Sack,
2004; Lewis, 2022; Burri and Kaufmann, 2026a, 2026b) combined with local
projections (Jorda, 2005). Identification exploits the difference in
variance between policy event days and control days to construct
instruments for the endogenous variables.

## Usage

``` r
hetiv(
  y,
  O,
  X = NULL,
  Ind,
  P,
  H,
  E = 1,
  norm = 1,
  interact = FALSE,
  cum = FALSE,
  cov_type = "HC3",
  hpredict = 0,
  details = FALSE
)
```

## Arguments

- y:

  Numeric matrix of stationary outcome variables (T x N). The effect on
  the first variable in each dimension is normalized to unity at
  horizon 0. These variables are also used to construct
  heteroskedasticity-based instruments and for recursive ordering to
  identify multiple dimensions.

- O:

  Numeric matrix of information set variables (T x M). May be identical
  to `y`. Included as lags 1 through `P`.

- X:

  Numeric matrix of deterministic variables (T x K). For example, time
  trend, seasonal dummies or other deterministic controls. Included as
  is (no lags). A constant is included by default.

- Ind:

  Integer vector of length T, event indicator:

  - `0` Control day (no event)

  - `1` Policy day (event)

  - `2` Contaminated control day (excluded from estimation)

- P:

  Integer. Maximum lag order for the information set. Set to `0` for no
  lags (regression on deterministic terms only).

- H:

  Integer. Maximum horizon (in periods) up to which IRFs are estimated.

- E:

  Integer. Number of shock dimensions to identify via recursive
  ordering.

- norm:

  Numeric scalar. Normalize the impact response of the first variable to
  a specific value. Set to `1` for standard unit-effect normalization.

- interact:

  Logical. If `TRUE`, lagged information set variables are interacted
  with event/non-event dummies.

- cum:

  Logical vector of length N. For each variable in `y`, whether to
  report the cumulative impulse response instead of the level response.
  If only one provided, applied to all impulse responses.

- cov_type:

  Covariance estimator for local-projection standard errors: `"HC3"`
  (default) for heteroskedasticity-robust standard errors or `"NW"` for
  Newey-West HAC standard errors. `"HC3"` is the default because Montiel
  Olea et al. (2025) show that heteroskedasticity-robust standard errors
  suffice for local-projection impulse responses under weak conditions,
  even though multi-step forecast errors are typically serially
  correlated. `"NW"` remains available as an optional HAC robustness
  check.

- hpredict:

  Integer. Horizon of the responses and forecast errors used for shock
  prediction, on the same scale as the row labels of `irf`. Defaults to
  0 (impact response and one-step-ahead residual). Must lie between 0
  and `H - 1`.

- details:

  Logical. If `TRUE`, saves detailed IV results and predicts the shocks,
  which is slightly slower. If `FALSE`, returns only impulse responses
  and standard errors (e.g. for bootstrap).

## Value

A named list. If `details = FALSE`, it contains `irf`, `se` and
`Method`. If `details = TRUE`, it contains:

- `irf`:

  Array (H x N x E) of estimated impulse responses.

- `se`:

  Array (H x N x E) of local-projection standard errors.

- `Shocks`:

  Predicted structural shocks on event days from
  [`kfpredict()`](https://dankaufmann.github.io/hetiv/reference/kfpredict.md),
  based on the responses at horizon `hpredict`.

- `IVRes`:

  List of `ivreg` model objects, one per horizon, variable, and shock
  dimension.

- `OLSRes`:

  List of OLS model objects used for residual-based covariance
  estimation, one per outcome variable.

- `ORTHRes`:

  List of OLS models used to orthogonalize the variable underlying the
  instrument, one per shock dimension.

- `FirstStage`:

  For `E = 1`: data frame with one row per horizon `h` giving the number
  of observations `Tobs`, the homoskedastic first-stage F statistic
  `F_hom` (the squared t statistic of the instrument; Stock and
  Yogo, 2005) and the robust effective F statistic `F_rob` of Montiel
  Olea and Pflueger (2013), computed with the covariance estimator in
  `cov_type`. The statistics do not depend on the outcome variable.
  `NULL` for `E > 1`.

- `WeakTest`:

  For `E > 1`: data frame with one row per horizon `h` and outcome
  variable `i` giving the results of
  [`gweakivtest()`](https://dankaufmann.github.io/hetiv/reference/gweakivtest.md)
  for the system with `E` endogenous policy variables and `E`
  instruments: `gmin` (generalized minimum eigenvalue statistic, which
  depends on the horizon only through the sample), `cv` and
  `cv_simplified` (sharp and simplified critical values, which depend on
  the outcome and horizon), and the Stock-Yogo statistic `sy_stat` with
  its critical value `sy_cv` (`NaN` unless over-identified). Rows where
  the outcome is a policy variable on impact are omitted. `NULL` for
  `E = 1`.

- `Obs`:

  Data frame with observation counts: `Tp` (policy days), `Tc` (control
  days), `To` (contaminated days), `Tt` (total used).

- `Method`:

  Character string `"Heteroscedasticity-IV"`.

- `et`:

  Data frame of OLS residuals on event days (used for covariance
  estimation and shock extraction).

- `Sig`:

  Covariance matrix of residuals on event days (used for shock
  extraction).

- `SigR`:

  Covariance matrix of residuals on control days, or `NA` if unavailable
  (used for shock extraction).

- `Psi`:

  Response matrix (N x E) used for shock prediction: the response at
  horizon `hpredict`, i.e. `irf[hpredict + 1, , ]`. With the default
  `hpredict = 0` this is the impact response.

- `WeakData`:

  Data frame of endogenous variables and instruments for the
  Lewis-Mertens (2025) weak instrument test.

## Details

For `E > 1`, identification is recursive and order-dependent: the column
order of `y` defines both the shock ordering and the normalization
variable for each shock dimension.

## References

Burri, M. and D. Kaufmann (2026a). Measuring monetary policy shocks.
IRENE Working Papers 24-03, IRENE Institute of Economic Research,
University of Neuchatel.

Burri, M. and D. Kaufmann (2026b). Multiple monetary policy shocks from
daily data: A heteroskedasticity IV approach. *Economics Letters*, 268,
113091.
[doi:10.1016/j.econlet.2026.113091](https://doi.org/10.1016/j.econlet.2026.113091)

Jorda, O. (2005). Estimation and inference of impulse responses by local
projections. *American Economic Review*, 95(1), 161-182.

Lewis, D. J. (2022). Robust inference in models identified via
heteroskedasticity. *Review of Economics and Statistics*, 104(3),
510-524.

Lewis, D. J. and Mertens, K. (2025). A robust test for weak instruments
for 2SLS with multiple endogenous regressors. *Review of Economic
Studies*, DOI: 10.1093/restud/rdaf103

Montiel Olea, J. L., M. Plagborg-Moller, E. Qian, and C. K. Wolf (2025).
Local projections or VARs? A primer for macroeconomists. *NBER Working
Paper* No. 33871.

Rigobon, R. (2003). Identification through heteroskedasticity. *Review
of Economics and Statistics*, 85(4), 777-792.

Rigobon, R. and Sack, B. (2004). The impact of monetary policy on asset
prices. *Journal of Monetary Economics*, 51(8), 1553-1575.

## Examples

``` r
set.seed(1)
y <- matrix(rnorm(80), ncol = 2)
Ind <- rep(0L, nrow(y))
Ind[seq(5, nrow(y), by = 5)] <- 1L
res <- hetiv(y = y, O = y, Ind = Ind, P = 1, H = 3, details = TRUE)
#> Warning: Some estimated shock variances are negative (weak heteroskedasticity); using absolute values for scaling.
dim(res$irf)
#> [1] 3 2 1
```
