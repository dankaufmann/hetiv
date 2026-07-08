# First-difference transformation

Computes the first difference of a numeric vector or time series, i.e.
`x_t - x_{t-1}`. Gaps in dates are ignored; the difference is always
taken with respect to the immediately preceding observation.

## Usage

``` r
firstdiff(TS)
```

## Arguments

- TS:

  A numeric vector or time series.

## Value

A numeric vector of the same length as `TS`, with `NA` in the first
position.

## Examples

``` r
firstdiff(c(1, 3, 6))
#> [1] NA  2  3
```
