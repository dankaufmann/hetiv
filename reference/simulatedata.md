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

## Examples

``` r
N <- 2
Phi <- array(0, dim = c(N, N, 1))
Phi[, , 1] <- diag(c(0.4, 0.2))
simulatedata(
  Phi = Phi, SigE = 2, PsiE = matrix(c(1, 0.5), N, 1),
  PsiR = diag(N), Nobs = 20, Nbin = 5, N = N, R = N, E = 1,
  Nevn = 5, P = 1, eDist = 0, seed = 1
)
#> $y
#>             [,1]         [,2]
#>  [1,]  0.3284413  1.613786725
#>  [2,]  0.6188056  0.219969618
#>  [3,]  0.9858469  0.431665535
#>  [4,]  0.9701201  0.032528066
#>  [5,] -0.1083364 -1.466051970
#>  [6,]  1.4684466 -0.708204957
#>  [7,]  0.9772219 -0.535930945
#>  [8,] -0.2303518 -0.166499586
#>  [9,] -2.3068406  1.066725455
#> [10,] -0.8489524  0.450947313
#> [11,] -0.3845146 -0.074334134
#> [12,] -0.1699961 -0.268228507
#> [13,]  0.8758378  0.643317674
#> [14,]  1.1715563  0.685326733
#> [15,]  4.1350607  0.984578097
#> [16,]  2.5730017 -0.510579538
#> [17,]  1.8113370  0.262466055
#> [18,]  0.7990998  0.821026135
#> [19,] -1.6697118  0.051859015
#> [20,] -1.8209643  0.005026851
#> 
#> $IndE
#>       [,1]
#>  [1,]    0
#>  [2,]    0
#>  [3,]    0
#>  [4,]    0
#>  [5,]    1
#>  [6,]    0
#>  [7,]    0
#>  [8,]    0
#>  [9,]    0
#> [10,]    1
#> [11,]    0
#> [12,]    0
#> [13,]    0
#> [14,]    0
#> [15,]    1
#> [16,]    0
#> [17,]    0
#> [18,]    0
#> [19,]    0
#> [20,]    1
#> 
#> $eR
#>              [,1]        [,2]
#>  [1,] -0.82046838  1.35867955
#>  [2,]  0.48742905 -0.10278773
#>  [3,]  0.73832471  0.38767161
#>  [4,]  0.57578135 -0.05380504
#>  [5,] -0.30538839 -1.37705956
#>  [6,]  1.51178117 -0.41499456
#>  [7,]  0.38984324 -0.39428995
#>  [8,] -0.62124058 -0.05931340
#>  [9,] -2.21469989  1.10002537
#> [10,]  1.12493092  0.76317575
#> [11,] -0.04493361 -0.16452360
#> [12,] -0.01619026 -0.25336168
#> [13,]  0.94383621  0.69696338
#> [14,]  0.82122120  0.55666320
#> [15,]  0.59390132 -0.68875569
#> [16,]  0.91897737 -0.70749516
#> [17,]  0.78213630  0.36458196
#> [18,]  0.07456498  0.76853292
#> [19,] -1.98935170 -0.11234621
#> [20,]  0.61982575  0.88110773
#> 
#> $eE
#>             [,1]
#>  [1,]  0.0000000
#>  [2,]  0.0000000
#>  [3,]  0.0000000
#>  [4,]  0.0000000
#>  [5,] -0.1909961
#>  [6,]  0.0000000
#>  [7,]  0.0000000
#>  [8,]  0.0000000
#>  [9,]  0.0000000
#> [10,] -1.0511471
#> [11,]  0.0000000
#> [12,]  0.0000000
#> [13,]  0.0000000
#> [14,]  0.0000000
#> [15,]  3.0725369
#> [16,]  0.0000000
#> [17,]  0.0000000
#> [18,]  0.0000000
#> [19,]  0.0000000
#> [20,] -1.7729054
#> 
#> $e
#>              [,1]         [,2]
#>  [1,] -0.82046838  1.358679552
#>  [2,]  0.48742905 -0.102787727
#>  [3,]  0.73832471  0.387671612
#>  [4,]  0.57578135 -0.053805041
#>  [5,] -0.49638444 -1.472557583
#>  [6,]  1.51178117 -0.414994563
#>  [7,]  0.38984324 -0.394289954
#>  [8,] -0.62124058 -0.059313397
#>  [9,] -2.21469989  1.100025372
#> [10,]  0.07378387  0.237602222
#> [11,] -0.04493361 -0.164523596
#> [12,] -0.01619026 -0.253361680
#> [13,]  0.94383621  0.696963375
#> [14,]  0.82122120  0.556663199
#> [15,]  3.66643821  0.847512750
#> [16,]  0.91897737 -0.707495157
#> [17,]  0.78213630  0.364581962
#> [18,]  0.07456498  0.768532925
#> [19,] -1.98935170 -0.112346212
#> [20,] -1.15307961 -0.005344952
#> 
#> $Phi
#> , , 1
#> 
#>      [,1] [,2]
#> [1,]  0.4  0.0
#> [2,]  0.0  0.2
#> 
#> 
#> $PsiE
#>      [,1]
#> [1,]  1.0
#> [2,]  0.5
#> 
#> $PsiR
#>      [,1] [,2]
#> [1,]    1    0
#> [2,]    0    1
#> 
#> $SigE
#> [1] 2
#> 
```
