# Plot impulse responses from a single estimation approach

Produces a panel of IRF plots for a single set of impulse responses,
with optional confidence bands. Returns a list of `ggplot` objects (one
per variable-shock combination) that can be arranged with e.g.
[`cowplot::plot_grid()`](https://wilkelab.org/cowplot/reference/plot_grid.html).

## Usage

``` r
plotirf(IRFest, IRFse = NULL, HTick, Labels, ci = c(0.9, 0.95))
```

## Arguments

- IRFest:

  Array (H x N x E) of impulse response point estimates. Dimnames on the
  first dimension must be the horizon labels (numeric).

- IRFse:

  Array (H x N x E) of standard errors, or `NULL` to suppress confidence
  bands.

- HTick:

  Integer. Step size for x-axis tick marks.

- Labels:

  Character vector of length N. Variable names used as panel titles.

- ci:

  Numeric vector of confidence levels in (0, 1). One shaded band is
  drawn per level, with wider bands rendered at lower opacity. Defaults
  to `c(0.90, 0.95)`.

## Value

A list of `ggplot` objects, one per variable-shock combination, ordered
by shock (outer loop) then variable (inner loop).

## Examples

``` r
irf <- array(c(1, 0.5, 0.2, 0.1), dim = c(4, 1, 1))
dimnames(irf)[[1]] <- 0:3
se <- array(0.1, dim = dim(irf), dimnames = dimnames(irf))
plotirf(irf, se, HTick = 1, Labels = "Output")
#> [[1]]

#> 
```
