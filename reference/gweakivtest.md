# Generalised weak instrument test

Tests for weak instruments with multiple endogenous regressors using the
generalised minimum eigenvalue statistic of Lewis and Mertens (2025).
The test is robust to heteroskedasticity and autocorrelation and nests
the classical Stock-Yogo (2005) test as a special case. The function is
a direct port from the Matlab codes by Lewis and Mertens (2025)

## Usage

``` r
gweakivtest(
  y,
  Y,
  X,
  Z,
  cov_type = "EHW",
  alfa = 0.05,
  tau = 0.1,
  points = 1000L,
  target = "beta",
  crit = "abs"
)
```

## Arguments

- y:

  Regressand (T x 1 numeric vector or matrix).

- Y:

  Endogenous regressors (T x N numeric matrix).

- X:

  Exogenous regressors (T x Nx numeric matrix). A constant column is
  added automatically if one is absent or if `X` has zero columns.

- Z:

  Instruments (T x K numeric matrix). Requires K \>= N.

- cov_type:

  HAR covariance estimator: `"EHW"` (Eicker-Huber-White, default) or
  `"NW"` (Newey-West with Lazarus et al. (2018) bandwidth).

- alfa:

  Nominal significance level (default `0.05`).

- tau:

  Bias tolerance: maximum acceptable relative (or absolute, see `crit`)
  bias of the 2SLS estimator (default `0.10`).

- points:

  Number of random starting points for the Stiefel manifold optimisation
  used to compute the sharp critical value when K \> N + 1 (default
  `1000`).

- target:

  Either `"beta"` (default) to test the full coefficient vector, or a
  positive integer `j <= N` to target the single coefficient `beta_j`.

- crit:

  Bias criterion: `"abs"` (absolute bias, default) or `"rel"` (relative
  bias). `"abs"` requires the error covariance matrix; `"rel"` does not.

## Value

A named list with the following elements:

- `nobs`:

  Number of complete observations used.

- `beta_2SLS`:

  2SLS point estimate(s).

- `target`:

  Description of the targeted parameter.

- `criterion`:

  Bias criterion used (`"abs"` or `"rel"`).

- `gmin_generalized`:

  Generalised minimum eigenvalue test statistic (Lewis-Mertens).

- `gmin_generalized_critical_value`:

  Sharp critical value via Stiefel optimisation.

- `gmin_generalized_critical_value_simplified`:

  Conservative simplified critical value (closed-form bound).

- `stock_yogo_test_statistic`:

  Stock-Yogo test statistic under the Nagar approximation.

- `stock_yogo_critical_value_nagar`:

  Stock-Yogo critical value under the Nagar approximation.

## References

Lazarus, E., Lewis, D. J., Stock, J. H. and Watson, M. W. (2018). HAR
Inference: Recommendations for Practice. *Journal of Business & Economic
Statistics*, 36(4), 541–559.

Lewis, D. J. and Mertens, K. (2025). A Robust Test for Weak Instruments
for 2SLS with Multiple Endogenous Regressors. *The Review of Economic
Studies*, DOI: 10.1093/restud/rdaf103.

Stock, J. H. and Yogo, M. (2005). Testing for Weak Instruments in Linear
IV Regression. In D. W. K. Andrews and J. H. Stock (Eds.),
*Identification and Inference for Econometric Models: Essays in Honor of
Thomas Rothenberg*, pp. 80–108. Cambridge University Press.
