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
  cov_type = "HC3",
  hpredict = 0,
  pooled = FALSE,
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

- pooled:

  Logical. If `FALSE` (default), the shock prediction uses the forecast
  errors and responses at horizon `hpredict` only. If `TRUE`, the
  forecast errors of all horizons `0, ..., hpredict` are stacked and
  used jointly (pooled prediction), which weights every horizon by its
  signal-to-noise ratio; see
  [`kfpredict()`](https://dankaufmann.github.io/hetiv/reference/kfpredict.md).
  With `hpredict = 0` both options coincide.

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

  Matrix of forecast errors on event days used for shock prediction, one
  column per variable and horizon (`y1.h0`, ..., `yN.h0`, `y1.h1`, ...):
  horizon `hpredict` only if `pooled = FALSE`, horizons
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
#> [1,] 1.0000000 5.329117e-17
#> [2,] 0.2783758 1.000000e+00
```
