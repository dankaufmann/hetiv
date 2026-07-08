# Log-difference transformation

Computes the first log-difference of a numeric vector or time series,
i.e. `log(x_t) - log(x_{t-1})`. Gaps in dates are ignored; the
difference is always taken with respect to the immediately preceding
observation.

## Usage

``` r
logdiff(TS)
```

## Arguments

- TS:

  A numeric vector or time series.

## Value

A numeric vector of the same length as `TS`, with `NA` in the first
position.

## Examples

``` r
logdiff(c(1, 2, 4))
#> [1]        NA 0.6931472 0.6931472
```
