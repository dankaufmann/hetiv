# Simulate VAR data with heteroskedastic shocks and given parameters

Simulates data from a VAR(P) model with two types of structural shocks:
regular shocks (always present) and event shocks (occurring every `Nevn`
periods). The shock distribution can be normal, Student-t, or
GARCH(1,1).

## Usage

``` r
simulatedata(Phi, SigE, PsiE, PsiR, Nobs, Nbin, N, R, E, Nevn, P, eDist, seed)
```

## Arguments

- Phi:

  Array of VAR coefficient matrices, dimension `N x N x P`.

- SigE:

  Variance of the event shocks. This applies to all E shocks.

- PsiE:

  Impact matrix for event shocks, dimension `N x E`.

- PsiR:

  Impact matrix for regular shocks, dimension `N x R`.

- Nobs:

  Number of observations to retain after burn-in. The returned data have
  exactly `Nobs` rows.

- Nbin:

  Number of burn-in observations. These are simulated to initialise the
  VAR but are discarded before returning.

- N:

  Number of variables in the VAR.

- R:

  Number of regular shocks.

- E:

  Number of event shocks. For GARCH shock distributions
  (`eDist = c(alpha, beta)`), only `E = 1` is supported.

- Nevn:

  Event frequency: an event shock occurs every `Nevn` periods. Set to
  `0` to suppress event shocks entirely (no heteroskedasticity).

- P:

  VAR lag order.

- eDist:

  Shock distribution. Use `0` for standard normal; a positive integer
  for Student-t with that many degrees of freedom; or a numeric vector
  `c(alpha, beta)` for GARCH(1,1) with ARCH effect `alpha` and
  persistence `beta`.

- seed:

  Integer seed passed to
  [`set.seed()`](https://rdrr.io/r/base/Random.html) for
  reproducibility. Use `NA` to skip seeding.

## Value

A named list with components:

- y:

  Simulated VAR data, dimension `Nobs x N` (burn-in discarded).

- IndE:

  Event indicator matrix, dimension `Nobs x 1`.

- eR:

  Simulated regular shocks, dimension `Nobs x R`.

- eE:

  Simulated event shocks, dimension `Nobs x E`.

- e:

  Composite structural shocks, dimension `Nobs x N`.

- Phi:

  VAR coefficient array (returned unchanged).

- PsiE:

  Event shock impact matrix (returned unchanged).

- PsiR:

  Regular shock impact matrix (returned unchanged).

- SigE:

  Event shock variance (returned unchanged).
