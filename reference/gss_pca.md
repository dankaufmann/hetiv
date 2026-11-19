# Monetary policy surprise factors from a panel of financial market series

Estimates one or two factors from a panel of event-window changes of
financial market series (rows are events, columns are series), following
Gurkaynak, Sack and Swanson (2005):

- `n_factors = 1`: a single common factor ("overall" surprise), the
  first principal component of the standardised panel.

- `n_factors = 2`: two factors rotated such that the second factor
  ("path") has a zero loading on the first column of `X`. The first
  factor ("target") is then the only factor that moves the first series.

## Usage

``` r
gss_pca(
  X,
  n_factors = 2,
  norm_anchor = seq_len(n_factors),
  tol = 1e-08,
  maxit = 500
)

# S3 method for class 'gss_pca'
print(x, ...)
```

## Arguments

- X:

  Numeric matrix or data frame (T x N) of event-window changes, with
  only the financial market series as columns. For `n_factors = 2`, the
  first column defines the target restriction (the path factor has a
  zero loading on it).

- n_factors:

  Integer. `1` (overall factor) or `2` (target and path factors,
  default).

- norm_anchor:

  Integer vector of length `n_factors` with the column positions in `X`
  against which the factors are normalised. Defaults to `1` or
  `c(1, 2)`. For `n_factors = 2`, `norm_anchor[2]` cannot be `1`.

- tol:

  Numeric. Convergence tolerance of the EM iterations.

- maxit:

  Integer. Maximum number of EM iterations.

- x:

  An object of class `gss_pca`.

- ...:

  Ignored.

## Value

A list of class `gss_pca` with elements:

- `factors`:

  Matrix T x `n_factors` with columns `"overall"` or `"target"` and
  `"path"`, in the units of the anchor series.

- `loadings`:

  Matrix N x `n_factors` of loadings of the standardised series on the
  (rotated and normalised) factors.

- `var_explained`:

  Share of the variance of the standardised panel explained by each
  unrotated principal component.

- `restriction`:

  For `n_factors = 2`: loading of the path factor on the first series
  (zero up to numerical error); `NULL` otherwise.

- `norm_anchor`:

  Named character vector with the anchor series of each factor.

- `norm_coef`:

  Named vector with the regression coefficient to which each factor was
  calibrated.

- `mu`, `sd`:

  Means and standard deviations of the series used for standardisation.

- `X_filled`:

  Standardised panel with missing cells replaced by their EM fit.

- `n_missing`, `iters`, `converged`, `call`:

  Number of missing cells, number of EM iterations, convergence flag,
  and the call.

## Details

Missing values are handled by EM principal components (Stock and Watson,
2002; Swanson, 2021): missing cells are initialised at the mean and
iteratively replaced by their fit from the `n_factors` principal
components until convergence.

The factors are normalised against series in `X` through regression
coefficients, as in Gurkaynak, Sack and Swanson (2005):

- `n_factors = 1`: the factor is rescaled such that the coefficient of a
  regression of `X[, norm_anchor[1]]` on the factor is one (a one-unit
  move of the factor is a one-unit move of its anchor series).

- `n_factors = 2`: the target factor is rescaled in the same way against
  `X[, norm_anchor[1]]`. The path factor is rescaled such that its
  coefficient on `X[, norm_anchor[2]]` equals the coefficient of the
  target factor on the same series.

The anchor regressions use the rows where the anchor series is observed.

## References

Gurkaynak, R. S., Sack, B. and Swanson, E. T. (2005). Do actions speak
louder than words? The response of asset prices to monetary policy
actions and statements. *International Journal of Central Banking*,
1(1), 55-93.

Stock, J. H. and Watson, M. W. (2002). Macroeconomic forecasting using
diffusion indexes. *Journal of Business & Economic Statistics*, 20(2),
147-162.

Swanson, E. T. (2021). Measuring the effects of Federal Reserve forward
guidance and asset purchases on financial markets. *Journal of Monetary
Economics*, 118, 32-53.

## Examples

``` r
# Simulated panel: 5 series driven by a target and a path factor; the first
# series responds to the target factor only
set.seed(1)
f <- matrix(rnorm(400), ncol = 2)
L <- cbind(c(1, 0.8, 0.6, 0.4, 0.2), c(0, 0.3, 0.6, 0.8, 1))
X <- f %*% t(L) + matrix(rnorm(1000, sd = 0.1), ncol = 5)
colnames(X) <- paste0("r", 1:5)
X[3, 4] <- NA
res <- gss_pca(X, n_factors = 2)
res
#> gss_pca: 2 factor(s) [target, path], 200 events x 5 series, 1 missing cells
#> var explained (PC1..PC2): 76.5%, 22.7%
#> path loading on 'r1' (target restriction, should be 0): 2.75e-16
#> normalisation: target -> coef 1.000 on 'r1'; path -> coef 0.794 on 'r2'
#> loadings:
#>    target  path
#> r1  1.061 0.000
#> r2  0.980 0.980
#> r3  0.701 1.932
#> r4  0.411 2.381
#> r5  0.165 2.556
head(res$factors)
#>          target        path
#> [1,] -0.6042134  0.09805403
#> [2,]  0.2541628  0.63802917
#> [3,] -0.8198619  0.54913400
#> [4,]  1.6079106 -0.15266963
#> [5,]  0.2529482 -0.86567612
#> [6,] -0.9417266  0.91834664
```
