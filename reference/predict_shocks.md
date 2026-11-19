# Extract monetary policy shocks via Kalman filter

Extracts structural shocks from reduced-form residuals using the Kalman
filter formula. Given the response matrix `Psi` and the covariance
matrices of residuals on event and control days, the function recovers
the latent shocks for each observation. Optionally rescales shocks to
unit variance.

## Usage

``` r
predict_shocks(
  Sig,
  SigR,
  Psi,
  et,
  tol = sqrt(.Machine$double.eps),
  scale = TRUE
)
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

A numeric matrix (T x E) of extracted structural shocks, with the same
number of rows as `et` and one column per shock dimension. Rows with
missing residuals give missing shocks.

## Details

The prediction is the linear projection of the shocks on the residuals,
`eps_t = Sig_eps Psi' Sig^{-1} e_t`, up to the scaling described under
`scale`. It applies to any set of forecast errors that is linear in the
shocks:

- One horizon: `Psi` is the N x E response at horizon h and `et` holds
  the T x N forecast errors of `y_{t+h}` given information dated
  `t - 1`. With h = 0 these are the impact responses and the one-step
  residuals.

- Pooled horizons: the forecast errors of several horizons are stacked,
  `Psi` has one block of N rows per horizon and `et` one block of N
  columns per horizon, in the same order. The pooled prediction weights
  every horizon by its signal-to-noise ratio and accounts for the
  correlation of the forecast errors across horizons. It has a weakly
  smaller mean squared error than the prediction from any single
  horizon, and coincides with the impact prediction when all responses
  are complete on impact.

For convenience, `Psi` may be supplied as an array N x E x M and `et` as
an array T x N x M (one slice per horizon); both are then stacked
internally. `Sig` and `SigR` must be of dimension N\*M in that case.

The pooled prediction requires a nonsingular covariance of the stacked
forecast errors. A sufficient condition is that at least as many
non-policy shocks as variables drive the data (the one-day non-policy
covariance has full rank). If fewer shocks are present, as in a stylised
simulation with fewer regular shocks than variables, `Sig` is singular
for stacked horizons, and the generalized inverse yields a regularized
prediction that is sensitive to estimation error in `Psi`. A warning is
issued when `Sig` is (nearly) singular, including when the number of
complete rows of `et` does not exceed its number of columns.

## References

Burri, M. and Kaufmann, D. (2026a). Measuring monetary policy shocks.
IRENE Working Papers 24-03, IRENE Institute of Economic Research,
University of Neuchatel.

## Examples

``` r
Sig <- diag(2)
SigR <- diag(c(0.5, 0.5))
Psi <- matrix(c(1, 0.5), nrow = 2)
et <- matrix(rnorm(20), ncol = 2)
predict_shocks(Sig = Sig, SigR = SigR, Psi = Psi, et = et)
#>              [,1]
#>  [1,]  0.41915187
#>  [2,] -0.08720511
#>  [3,]  0.13355763
#>  [4,] -0.30746878
#>  [5,] -1.13657435
#>  [6,]  0.01761954
#>  [7,]  0.55182712
#>  [8,] -0.84338832
#>  [9,] -0.92125158
#> [10,]  1.24787258

# Pooled prediction from two horizons (arrays with one slice per horizon)
Psi2 <- array(c(Psi, 0.5 * Psi), dim = c(2, 1, 2))
et2 <- array(rnorm(40), dim = c(10, 2, 2))
et2m <- matrix(et2, nrow = 10)
predict_shocks(Sig = var(et2m), SigR = 0.5 * var(et2m), Psi = Psi2, et = et2)
#>             [,1]
#>  [1,]  1.2184774
#>  [2,]  0.6732302
#>  [3,]  0.3078382
#>  [4,] -0.7063875
#>  [5,] -0.1702423
#>  [6,] -0.4185534
#>  [7,]  1.7734930
#>  [8,]  1.1025325
#>  [9,]  0.7177884
#> [10,]  0.4093612
```
