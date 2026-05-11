# Normalize a time series

Standardizes a numeric vector to zero mean and unit standard deviation
(z-score normalization), ignoring missing values.

## Usage

``` r
normalize(x)
```

## Arguments

- x:

  A numeric vector or time series.

## Value

A numeric vector of the same length as `x`, with mean 0 and standard
deviation 1.

## Examples

``` r
normalize(c(1, 2, 3, 4, 5))
#> [1] -1.2649111 -0.6324555  0.0000000  0.6324555  1.2649111
```
