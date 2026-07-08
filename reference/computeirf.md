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

## Examples

``` r
Psi <- matrix(c(1, 0.5), nrow = 2)
Phi <- array(0, dim = c(2, 2, 1))
Phi[, , 1] <- matrix(c(0.4, 0.1, 0, 0.3), 2, 2)
computeirf(Psi = Psi, Phi = Phi, H = 4, cum = FALSE)
#> , , 1
#> 
#>    [,1]   [,2]
#> 0 1.000 0.5000
#> 1 0.400 0.2500
#> 2 0.160 0.1150
#> 3 0.064 0.0505
#> 
```
