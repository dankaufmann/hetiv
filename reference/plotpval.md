# Plot p-values across horizons

Produces a panel of p-value plots across IRF horizons, with optional
horizontal reference lines at user-supplied significance levels. Returns
a list of `ggplot` objects (one per variable-shock combination) that can
be arranged with e.g.
[`cowplot::plot_grid()`](https://wilkelab.org/cowplot/reference/plot_grid.html).

## Usage

``` r
plotpval(pvals, HTick, Labels, sigLevels = c(0.05, 0.1))
```

## Arguments

- pvals:

  Array (H x N x E) of p-values. Dimnames on the first dimension must be
  the horizon labels (numeric).

- HTick:

  Integer. Step size for x-axis tick marks.

- Labels:

  Character vector of length N. Variable names used as panel titles.

- sigLevels:

  Numeric vector of significance levels at which to draw horizontal
  reference lines. The first value is drawn as a solid line, the second
  as dashed, and any further values as dotted. Defaults to
  `c(0.05, 0.10)`.

## Value

A list of `ggplot` objects, one per variable-shock combination, ordered
by shock (outer loop) then variable (inner loop).

## Examples

``` r
pvals <- array(c(0.2, 0.08, 0.04), dim = c(3, 1, 1))
dimnames(pvals)[[1]] <- 0:2
plotpval(pvals, HTick = 1, Labels = "Output")
#> [[1]]

#> 
```
