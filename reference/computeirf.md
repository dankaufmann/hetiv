# Compute impulse response functions for a VAR(P)

Recursively computes IRFs given an impact matrix and VAR coefficient
array. Optionally cumulates responses for selected variables.

## Usage

``` r
computeirf(Psi, Phi, H, cum)
```

## Arguments

- Psi:

  Impact matrix, dimension `N x R`.

- Phi:

  Array of VAR coefficient matrices, dimension `N x N x P`.

- H:

  Horizon (number of periods) for which to compute IRFs.

- cum:

  Logical scalar or logical vector of length `N`. If `TRUE` for variable
  `i`, the IRF is cumulated via `cumsum`. A single value is applied to
  all variables.

## Value

An array:

- irf:

  Array of impulse responses, dimension `H x N x R`. Row names are
  labelled `0` to `H-1`.
