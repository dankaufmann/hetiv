# Estimate impulse responses via recursive (Cholesky) local projections

Estimates impulse response functions (IRFs) using recursive zero
restrictions combined with local projections (Jorda, 2005). The shock of
dimension `e` is identified by regressing the outcome at horizon `h` by
OLS on the contemporaneous value of the `e`-th variable in `y`,
controlling for the contemporaneous values of variables `1, ..., e - 1`
in `y`, lags of the information set, and deterministic terms. This is
the local-projection analogue of a Cholesky decomposition and serves as
a benchmark for the heteroskedasticity-based estimator
[`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md).

## Usage

``` r
recols(
  y,
  O,
  X = NULL,
  Ind,
  P,
  H,
  E = 1,
  norm = 1,
  cum = FALSE,
  Hstep = 1,
  cov_type = "HC3",
  hpredict = 1,
  details = FALSE
)
```

## Arguments

- y:

  Numeric matrix of stationary outcome variables (T x N). The effect on
  the `e`-th variable in each dimension `e` is normalized to `norm` at
  horizon 0. The first `E` variables are also used to impose the
  recursive zero restrictions.

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

  At least one policy day is required.

- P:

  Integer. Maximum lag order for the information set. Set to `0` for no
  lags (regression on deterministic terms only).

- H:

  Integer. Maximum horizon (in periods) up to which IRFs are estimated.

- E:

  Integer. Number of shock dimensions to identify via recursive
  ordering.

- norm:

  Numeric scalar. Normalize the impact response of the `e`-th variable
  to shock `e` to a specific value. Set to `1` for standard unit-effect
  normalization.

- cum:

  Logical vector of length N. For each variable in `y`, whether to
  report the cumulative impulse response instead of the level response.
  If only one provided, applied to all impulse responses.

- Hstep:

  Integer. Step size between horizons. The default `1` estimates all
  horizons 0 through H - 1. Values greater than 1 estimate only the
  selected horizons.

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

  Integer. Forecast horizon of the residuals used for shock prediction.
  Defaults to 1 (one-step-ahead residual). Must be one of the estimated
  horizons in `seq(1, H, Hstep)`.

- details:

  Logical. If `TRUE`, saves detailed OLS results and extracts shocks,
  which is slightly slower. If `FALSE`, returns only impulse responses
  and standard errors (e.g. for bootstrap).

## Value

A named list. If `details = FALSE`, it contains `irf`, `se`, and
`Method`. If `details = TRUE`, it contains:

- `irf`:

  Array (H x N x E) of estimated impulse responses.

- `se`:

  Array (H x N x E) of local-projection standard errors.

- `Shocks`:

  Output of
  [`kfpredict()`](https://dankaufmann.github.io/hetiv/reference/kfpredict.md):
  predicted structural shocks.

- `OLSRes`:

  List of `lm` model objects, one per horizon, variable, and shock
  dimension.

- `Obs`:

  Data frame with observation counts: `Tp` (policy days), `Tc` (control
  days), `To` (contaminated days), `Tt` (total used).

- `Method`:

  Character string `"Recursive-OLS"`.

- `et`:

  Matrix of OLS residuals on event days (used for covariance estimation
  and shock extraction).

- `Sig`:

  Covariance matrix of residuals on event days.

- `SigR`:

  Covariance matrix of residuals on control days, or `NA` if
  unavailable.

- `Psi`:

  Impact matrix (N x E) at horizon `hpredict - 1`; with the default
  `hpredict = 1` this equals `irf[1, , ]`.

## Details

For `E > 1`, identification is recursive and order-dependent: the column
order of `y` defines both the shock ordering and the normalization
variable for each shock dimension. By construction, the impact response
of variable `e` to shock `e` equals `norm`, and the impact responses of
variables `1, ..., e - 1` to shock `e` are zero.

Unlike
[`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md), the
distinction between policy days (`Ind == 1`) and control days
(`Ind == 0`) does not affect the IRF estimates: both enter the local
projections, and only contaminated days (`Ind == 2`) are dropped. The
distinction matters only for the shock extraction returned when
`details = TRUE`.

## References

Jorda, O. (2005). Estimation and inference of impulse responses by local
projections. *American Economic Review*, 95(1), 161-182.

Montiel Olea, J. L., M. Plagborg-Moller, E. Qian, and C. K. Wolf (2025).
Local projections or VARs? A primer for macroeconomists. *NBER Working
Paper* No. 33871.

Plagborg-Moller, M. and C. K. Wolf (2021). Local projections and VARs
estimate the same impulse responses. *Econometrica*, 89(2), 955-980.

## See also

[`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md) for
heteroskedasticity-based identification and
[`proxyiv()`](https://dankaufmann.github.io/hetiv/reference/proxyiv.md)
for proxy-based identification.

## Examples

``` r
set.seed(1)
y <- matrix(rnorm(80), ncol = 2)
Ind <- rep(0L, nrow(y))
Ind[seq(5, nrow(y), by = 5)] <- 1L
res <- recols(y = y, O = y, Ind = Ind, P = 1, H = 3, E = 2)
res$irf[1, , ] # impact responses: unit diagonal, zero above it
#>           [,1]         [,2]
#> [1,] 1.0000000 -4.26246e-17
#> [2,] 0.2783758  1.00000e+00
```
