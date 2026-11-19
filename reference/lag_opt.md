# Optimal lag length of a VAR

Selects the lag length of a VAR (with a constant and, optionally,
exogenous variables) by an information criterion, using
[`vars::VAR()`](https://rdrr.io/pkg/vars/man/VAR.html). The result can
be used as the number of lags `P` of the information set in
[`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md),
[`proxyiv()`](https://dankaufmann.github.io/hetiv/reference/proxyiv.md),
[`recols()`](https://dankaufmann.github.io/hetiv/reference/recols.md)
and
[`svariv()`](https://dankaufmann.github.io/hetiv/reference/svariv.md).

## Usage

``` r
lag_opt(Y, X = NULL, crit = "AIC", lag.max = 15)
```

## Arguments

- Y:

  Numeric matrix or data frame (T x N) of endogenous variables.

- X:

  Numeric matrix or data frame (T x K) of exogenous variables, or `NULL`
  (default).

- crit:

  Character. Information criterion: `"AIC"` (default, following Montiel
  Olea et al., 2025), `"HQ"`, `"SC"` or `"FPE"`.

- lag.max:

  Integer. Largest lag length considered. Default 15.

## Value

Integer. The lag length selected by the criterion.

## Details

Requires the `vars` package, which is not installed automatically with
`hetiv`.

## References

Montiel Olea, J. L., Plagborg-Moller, M., Qian, E. and Wolf, C. K.
(2025). Local projections or VARs? A primer for macroeconomists. NBER
Macroeconomics Annual, forthcoming.

## Examples

``` r
if (requireNamespace("vars", quietly = TRUE)) {
  set.seed(1)
  Y <- matrix(rnorm(400), ncol = 2)
  for (t in 2:200) Y[t, ] <- 0.5 * Y[t - 1, ] + Y[t, ]
  lag_opt(Y, lag.max = 4)
}
#> [1] 1
```
