# Rescale a series to match the mean and variance of another series

Applies the linear transformation to `y` that matches the mean and the
standard deviation of `x` on their overlapping (jointly observed)
sample. Useful to compare measures of the same shock that are expressed
in different units, for example a predicted shock and a proxy.

## Usage

``` r
rescale_to_match(x, y)
```

## Arguments

- x:

  Numeric vector. Reference series (mean and variance target).

- y:

  Numeric vector of the same length as `x`. Series to rescale.

## Value

Numeric vector of the same length as `y`: `a + b * y`, with `a` and `b`
chosen such that the result has the mean and standard deviation of `x`
on the overlapping sample. Missing values of `y` remain missing.

## Examples

``` r
set.seed(1)
x <- rnorm(50, mean = 2, sd = 3)
y <- 0.1 * x + rnorm(50, sd = 0.1)
y2 <- rescale_to_match(x, y)
c(mean(y2), sd(y2))
#> [1] 2.301345 2.494182
```
