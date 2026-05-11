# Linearly interpolate gaps in a data frame

Fills missing values (`NA`) in every column of a data frame using linear
interpolation, up to a maximum of `gap` consecutive missing
observations. Gaps longer than `gap` are left as `NA`.

## Usage

``` r
filllinear(DF, gap)
```

## Arguments

- DF:

  A data frame or matrix with numeric columns.

- gap:

  Integer. Maximum number of consecutive `NA`s to fill.

## Value

The data frame `DF` with short gaps linearly interpolated.
