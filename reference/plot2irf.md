# Plot and compare impulse responses from two estimation approaches

Produces a panel of IRF plots comparing two sets of impulse responses
side-by-side. Each panel shows the point estimates and confidence bands
at level `ci` for both approaches. Returns a list of `ggplot` objects
(one per variable-shock combination) that can be arranged with e.g.
[`cowplot::plot_grid()`](https://wilkelab.org/cowplot/reference/plot_grid.html).

## Usage

``` r
plot2irf(IRF1, IRF1se, IRF2, IRF2se, HTick, Labels, ci = 0.9)
```

## Arguments

- IRF1:

  Array (H x N x E) of impulse responses for the first approach.
  Dimnames on the first dimension must be the horizon labels (numeric).

- IRF1se:

  Array (H x N x E) of standard errors for the first approach.

- IRF2:

  Array (H x N x E) of impulse responses for the second approach.

- IRF2se:

  Array (H x N x E) of standard errors for the second approach.

- HTick:

  Integer. Step size for x-axis tick marks.

- Labels:

  Character vector of length N. Variable names used as panel titles.

- ci:

  Numeric scalar in (0, 1). Confidence level for the shaded bands.
  Defaults to `0.90` (90% CI).

## Value

A list of `ggplot` objects, one per variable-shock combination, ordered
by shock (outer loop) then variable (inner loop).
