# Deprecated: extract shocks via Kalman filter

`kfpredict()` has been renamed to
[`predict_shocks()`](https://dankaufmann.github.io/hetiv/reference/predict_shocks.md).
It is kept for backward compatibility, issues a deprecation warning and
will be removed in a future release.

## Usage

``` r
kfpredict(Sig, SigR, Psi, et, tol = sqrt(.Machine$double.eps), scale = TRUE)
```

## Arguments

- Sig:

  Numeric matrix (M x M, with M = N for one horizon). Covariance matrix
  of the residuals in `et` on policy event days (e.g. `Sig` from
  [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md)).

- SigR:

  Numeric matrix (M x M). Covariance matrix of the residuals on control
  (non-event) days. Used to back out the implied shock variances when
  `scale = TRUE`.

- Psi:

  Numeric matrix (M x E) of responses of the residuals in `et` to the E
  structural shocks (e.g. `Psi` from
  [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md)),
  or an array N x E x M of responses at M horizons, which is stacked by
  rows.

- et:

  Numeric matrix or data frame (T x M) of forecast errors on event days
  (e.g. `et` from
  [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md)),
  or an array T x N x M of forecast errors at M horizons, which is
  stacked by columns. Rows correspond to time periods.

- tol:

  Numeric scalar. Tolerance for the generalized inverse
  ([`MASS::ginv()`](https://rdrr.io/pkg/MASS/man/ginv.html)). Clipped
  from below at `sqrt(.Machine$double.eps)`.

- scale:

  Logical. If `TRUE` (default), shocks are rescaled to unit variance
  using the implied shock variances recovered from `Sig` and `SigR`.
  This scaling assumes that structural shock variances are diagonal, so
  that `Sig - SigR = sum_i sigma_i Psi_i Psi_i'`. If `FALSE`, the raw
  Kalman filter projection is returned.

## Value

See
[`predict_shocks()`](https://dankaufmann.github.io/hetiv/reference/predict_shocks.md).
