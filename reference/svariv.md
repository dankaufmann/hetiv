# Estimate impulse responses via a proxy SVAR (SVAR-IV)

Estimates a reduced-form VAR(P) by OLS, identifies the impact matrix of
the first `E` structural shocks with external instruments (proxies)
following Mertens and Ravn (2013) and Stock and Watson (2018), and
computes impulse responses by recursion through the VAR. Optionally
computes confidence bands with a wild bootstrap following Mertens and
Ravn (2013).

## Usage

``` r
svariv(
  y,
  Z,
  X = NULL,
  P,
  H,
  E = 1,
  norm = 1,
  cum = FALSE,
  nboot = 0,
  ci = c(0.9, 0.95),
  cov_type = "HC0",
  details = FALSE
)
```

## Arguments

- y:

  Numeric matrix of VAR variables (T x N) without missing values. The
  first `E` columns are the policy variables on which the responses are
  normalized.

- Z:

  Numeric matrix of external instruments (T x M, `M >= E`). Missing
  values are allowed and drop the corresponding observation from the
  identification step. Set the proxy to zero (not `NA`) on days without
  an event.

- X:

  Numeric matrix of deterministic variables (T x K), or `NULL`
  (default). Included as-is (no lags) in addition to the constant, which
  is always included.

- P:

  Integer. Lag order of the VAR (`P >= 1`).

- H:

  Integer. Number of horizons: IRFs are reported for horizons
  `0, ..., H - 1`.

- E:

  Integer. Number of shocks to identify (`E < N`). Default `1`.

- norm:

  Numeric scalar. Normalize the impact response of variable `e` to shock
  `e` to this value. Default `1`.

- cum:

  Logical scalar or vector of length N. If `TRUE` for variable `i`, the
  cumulative IRF is reported. A single value is recycled to all
  variables. Default `FALSE`.

- nboot:

  Integer. Number of wild-bootstrap replications. `0` (default) skips
  the bootstrap, in which case `se`, `lower`, and `upper` are `NULL`.

- ci:

  Numeric vector of confidence levels in (0, 1) for the percentile
  bands. Default `c(0.90, 0.95)`.

- cov_type:

  Covariance estimator for the robust first-stage F statistic reported
  with `details = TRUE`: `"HC0"` (default, Eicker-Huber-White), `"HC1"`,
  `"HC3"`, or `"NW"` for Newey-West HAC.

- details:

  Logical. If `TRUE`, additionally returns the VAR estimates, residuals
  and first-stage statistics.

## Value

A named list. Always contains:

- `irf`:

  Array (H x N x E) of impulse responses, horizons labelled `0` to
  `H - 1`.

- `se`:

  Array (H x N x E) of bootstrap standard deviations, or `NULL` if
  `nboot = 0`.

- `lower`, `upper`:

  Arrays (H x N x E x length(ci)) of percentile bands, one slice per
  level in `ci` (labelled in the fourth dimension), or `NULL` if
  `nboot = 0`.

- `boot`:

  Array (H x N x E x nboot) of bootstrap draws (`NA` for failed
  replications), or `NULL` if `nboot = 0`.

- `Method`:

  Character string `"SVAR-IV"`.

With `details = TRUE`, additionally contains:

- `Psi`:

  Normalized impact matrix (N x E), equal to `irf[1, , ]` before
  cumulation.

- `Phi`:

  Array (N x N x P) of VAR slope coefficients.

- `C`:

  Vector of length N with the VAR intercepts.

- `D`:

  Matrix (N x K) of coefficients on `X`, or `NULL`.

- `U`:

  Matrix ((T - P) x N) of reduced-form residuals.

- `FStat`, `FStatRobust`, `R2`:

  Vectors of length E with the homoskedastic first-stage F statistic,
  the robust effective F statistic (Montiel Olea and Pflueger, 2013) and
  the R-squared of the regressions of each policy residual on the
  proxies (see Details).

- `Obs`:

  Data frame with `T` (sample size), `Test` (VAR estimation sample), and
  `Tiv` (identification sample).

- `nfail`:

  Number of failed bootstrap replications.

- `ci`:

  The confidence levels used for `lower` and `upper`.

## Details

The first `E` columns of `y` are the policy variables: shock `e` is the
structural shock to variable `e`, and the impact response of variable
`e` to shock `e` is normalized to `norm`. For `E > 1` the impact matrix
of the policy block is identified up to a lower-triangular (Cholesky)
rotation, so shock `e` has no impact effect on policy variables
`1, ..., e-1`; identification is therefore order-dependent, as in
[`proxyiv()`](https://dankaufmann.github.io/hetiv/reference/proxyiv.md)
with `recursive = TRUE`.

All columns of `Z` are used jointly as instruments in the first-stage
regressions of the policy residuals on the proxies. Rows of `Z` with
missing values are excluded from the identification step only: the VAR
and the reduced-form residual covariance `Sigma` are estimated on the
full sample, as in Mertens and Ravn (2013), whereas the 2SLS step that
recovers the relative impact responses uses only the observations with
non-missing proxies (the identification sample).

With `details = TRUE`, first-stage diagnostics are reported for every
policy residual: the homoskedastic F statistic of the proxies (`FStat`,
equal to the squared t statistic when `M = 1`; Stock and Yogo, 2005) and
the robust effective F statistic of Montiel Olea and Pflueger (2013)
(`FStatRobust`), computed with the covariance estimator in `cov_type`.
For `E > 1` these are equation-by-equation statistics; use
[`gweakivtest()`](https://dankaufmann.github.io/hetiv/reference/gweakivtest.md)
for a joint test.

The bootstrap follows Mertens and Ravn (2013): in every replication, the
VAR residuals and the proxies are multiplied by the same Rademacher
draw, artificial data are generated recursively from the first `P`
observations, and the VAR is re-estimated and re-identified. Each draw
is normalized in the same way as the point estimate. Percentile bands
are reported for every level in `ci`; `se` is the bootstrap standard
deviation and `boot` holds the raw draws, so the result can be passed
directly to
[`plotirf()`](https://dankaufmann.github.io/hetiv/reference/plotirf.md),
which draws recentred percentile bands from the draws. Replications in
which identification fails (e.g. a non-positive-definite policy block)
are dropped with a warning.

## References

Mertens, K. and Ravn, M. O. (2013). The dynamic effects of personal and
corporate income tax changes in the United States. *American Economic
Review*, 103(4), 1212-1247.

Montiel Olea, J. L. and Pflueger, C. (2013). A robust test for weak
instruments. *Journal of Business & Economic Statistics*, 31(3),
358-369.

Stock, J. H. and Watson, M. W. (2018). Identification and estimation of
dynamic causal effects in macroeconomics using external instruments.
*Economic Journal*, 128(610), 917-948.

Stock, J. H. and Yogo, M. (2005). Testing for weak instruments in linear
IV regression. In D. W. K. Andrews and J. H. Stock (Eds.),
*Identification and inference for econometric models: Essays in honor of
Thomas Rothenberg* (pp. 80-108). Cambridge University Press.

## See also

[`proxyiv()`](https://dankaufmann.github.io/hetiv/reference/proxyiv.md)
for the local-projection counterpart,
[`computeirf()`](https://dankaufmann.github.io/hetiv/reference/computeirf.md),
[`plotirf()`](https://dankaufmann.github.io/hetiv/reference/plotirf.md).

## Examples

``` r
set.seed(1)
Phi <- array(c(0.5, 0.1, -0.1, 0.4), dim = c(2, 2, 1))
Psi <- matrix(c(1, 0.5, 0, 1), 2, 2)
Nobs <- 400
eps <- matrix(rnorm(2 * Nobs), ncol = 2)
y <- matrix(0, Nobs, 2)
for (t in 2:Nobs) y[t, ] <- Phi[, , 1] %*% y[t - 1, ] + Psi %*% eps[t, ]
Z <- matrix(eps[, 1] + rnorm(Nobs), ncol = 1)
res <- svariv(y = y, Z = Z, P = 1, H = 6, nboot = 50, ci = 0.90)
dim(res$irf)
#> [1] 6 2 1
res$irf[1, , 1] # impact responses, normalized to 1 on the first variable
#> [1] 1.0000000 0.5608632
```
