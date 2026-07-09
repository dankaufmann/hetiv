# Extract monetary policy shocks via Kalman filter

Extracts structural shocks from reduced-form residuals using the Kalman
filter formula. Given the impact matrix `Psi` and the covariance
matrices of residuals on event and control days, the function recovers
the latent shocks for each observation. Optionally rescales shocks to
unit variance.

## Usage

``` r
kfpredict(Sig, SigR, Psi, et, tol = sqrt(.Machine$double.eps), scale = TRUE)
```

## Arguments

- Sig:

  Numeric matrix (N x N). Covariance matrix of reduced-form residuals on
  policy event days (e.g. `et` from
  [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md)).

- SigR:

  Numeric matrix (N x N). Covariance matrix of reduced-form residuals on
  control (non-event) days. Used to back out the implied shock variances
  when `scale = TRUE`.

- Psi:

  Numeric matrix (N x E). Impact matrix, i.e. the contemporaneous
  responses of all N variables to the E structural shocks (e.g. `Psi`
  from
  [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md)).

- et:

  Numeric matrix or data frame (T x N). Reduced-form residuals on event
  days (e.g. `et` from
  [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md)).
  Rows correspond to time periods, columns to variables.

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
number of rows as `et` and one column per shock dimension.

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
kfpredict(Sig = Sig, SigR = SigR, Psi = Psi, et = et)
#>              [,1]
#>  [1,] -0.53121758
#>  [2,]  0.29646687
#>  [3,]  1.11203917
#>  [4,] -0.74216125
#>  [5,]  0.87744535
#>  [6,]  0.38718522
#>  [7,]  0.26866947
#>  [8,] -0.37366525
#>  [9,] -0.15323607
#> [10,]  0.01922568
```
