# hetiv

[![coverage](https://github.com/dankaufmann/hetiv/actions/workflows/coverage.yaml/badge.svg)](https://github.com/dankaufmann/hetiv/actions/workflows/coverage.yaml)

An R package for measuring and identifying structural shocks using heteroskedasticity- and proxy-based instrumental variable (IV) methods with daily financial market data.

**Authors**: Daniel Kaufmann, Marc Burri, Valentin Grob

Comments, bug reports, and feature requests are welcome through the GitHub issue tracker.

## Documentation

Full documentation and a worked example are available at the package website: **https://dankaufmann.github.io/hetiv/**

The introductory vignette — covering data simulation, HET-IV, Proxy-IV, IRF plots, shock extraction, and weak instrument testing — is at: **https://dankaufmann.github.io/hetiv/articles/hetiv-introduction.html**

Replication files for Burri and Kaufmann (2026b), providing a real-world example, are available at: **https://github.com/dankaufmann/bk_2026_eclet_replication/**

`gweakivtest()` is a direct port of the Matlab files by Lewis and Mertens
(2025) available on **https://karelmertens.com/research/**.

## Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("dankaufmann/hetiv")
```

## Functions

| Function | Description |
|---|---|
| `hetiv()` | Estimates impulse response functions via heteroskedasticity-based IV local projections |
| `proxyiv()` | Estimates impulse response functions via proxy-IV local projections |
| `gweakivtest()` | Tests for weak instruments using the generalised minimum eigenvalue statistic of Lewis and Mertens (2025) |
| `kfpredict()` | Extracts structural shocks from reduced-form residuals via Kalman filter |
| `computeirf()` | Recursively computes IRFs from an impact matrix and VAR coefficients |
| `normalize()` | Standardizes a time series to zero mean and unit standard deviation |
| `logdiff()` | Computes the log-difference of a time series |
| `firstdiff()` | Computes the first difference of a time series |
| `filllinear()` | Linearly interpolates short gaps of missing values in a data frame |
| `plotirf()` | Plots IRFs with confidence bands for a single estimation approach |
| `plot2irf()` | Plots and compares IRFs from two estimation approaches |
| `plotpval()` | Plots p-values of IRF coefficients across horizons with significance reference lines |
| `simulatedata()` | Simulates VAR data with heteroskedastic event shocks |

## Usage

### Estimate impulse responses via heteroskedasticity-IV

```r
library(hetiv)

set.seed(1)
N <- 2
E <- 1
Phi <- array(0, dim = c(N, N, 1))
Phi[, , 1] <- diag(c(0.4, 0.2))
sim <- simulatedata(
  Phi = Phi, SigE = 2, PsiE = matrix(c(1, 0.5), N, 1),
  PsiR = diag(N), Nobs = 100, Nbin = 20, N = N, R = N, E = E,
  Nevn = 5, P = 1, eDist = 0, seed = 1
)

y <- sim$y
O <- y
Ind <- as.integer(sim$IndE[, 1])

res <- hetiv(y = y, O = O, Ind = Ind, P = 1, H = 20, E = 1, norm = 1, details = TRUE)

# With a linear time trend
X <- matrix(seq_len(nrow(y)), ncol = 1)
res <- hetiv(y = y, O = O, X = X, Ind = Ind, P = 1, H = 20, E = 1, norm = 1, details = TRUE)
```

### Estimate impulse responses via proxy-IV

```r
# Z: T x E matrix of external instruments (proxies)
Z <- matrix(sim$eE[, 1] + rnorm(nrow(y)), ncol = 1)
Z[Ind == 0, ] <- 0

res <- proxyiv(y = y, O = O, Z = Z, Ind = Ind, P = 1, H = 20, E = 1, norm = 1, details = TRUE)

# With recursive zero restrictions across shock dimensions
res <- proxyiv(y = y, O = O, Z = Z, Ind = Ind, P = 1, H = 20, E = 2,
               norm = 1, recursive = TRUE, details = TRUE)
```

### Plot impulse responses

```r
plots <- plotirf(IRFest = res$irf, IRFse = res$se,
                 HTick = 5, Labels = c("Var 1", "Var 2"))
cowplot::plot_grid(plotlist = plots)
```

### Test for weak instruments (Lewis-Mertens, 2025)

```r
# y: T x 1 outcome variable
# Y: T x N matrix of endogenous regressors
# X: T x Nx matrix of exogenous controls
# Z: T x K matrix of instruments (requires K >= N)

wt <- gweakivtest(y = y, Y = Y, X = X, Z = Z)

# Inspect test statistics and critical values
wt$gmin_generalized                            # generalised min-eigenvalue statistic
wt$gmin_generalized_critical_value             # Lewis-Mertens sharp critical value (Stiefel)
wt$gmin_generalized_critical_value_simplified  # conservative closed-form bound

# With Newey-West HAR standard errors and custom bias tolerance
wt <- gweakivtest(y = y, Y = Y, X = X, Z = Z, cov_type = "NW", tau = 0.10)

# hetiv() and proxyiv() return WeakData when details = TRUE,
# which contains the outcome and instrument columns needed for gweakivtest(). Note that the E endogenous variables have to be ordered first
res     <- hetiv(y = y, O = O, Ind = Ind, P = 1, H = 20, E = 1, details = TRUE)
weakdata <- res$WeakData

# y: outcome variable E + 1 (not used as endogenous regressor)
y_w <- weakdata[, paste0("y", E + 1)]
# Y: first E outcome variables (endogenous regressors)
Y_w <- weakdata[, paste0("y", 1:E)]
# Z: the E instruments
Z_w <- weakdata[, paste0("Z", 1:E), drop = FALSE]
# X: lagged ("o*"), deterministic ("x*"), and indicator ("i*") controls.
# gweakivtest() adds a constant term if one is missing.
ctrl <- startsWith(colnames(weakdata), "o") |
  startsWith(colnames(weakdata), "x") |
  startsWith(colnames(weakdata), "i")
X_w <- if (any(ctrl)) weakdata[, ctrl, drop = FALSE] else matrix(numeric(0), nrow(weakdata), 0)

wt <- gweakivtest(y = y_w, Y = Y_w, X = X_w, Z = Z_w)
```

### Extract structural shocks

```r
# Use estimates from hetiv() or proxyiv() to predict the unobserved underlying
# shocks with the Kalman filter.
shocks <- kfpredict(Sig = res$Sig, SigR = res$SigR,
                    Psi = res$Psi, et = res$et)
```

### Simulate data

```r
sim <- simulatedata(Phi = Phi, SigE = 4, PsiE = PsiE, PsiR = PsiR,
                    Nobs = 500, Nbin = 100, N = 2, R = 2, E = 1,
                    Nevn = 5, P = 1, eDist = 0, seed = 42)
```

## References

Burri, M. and D. Kaufmann (2026a). Measuring monetary policy shocks. IRENE Working Papers 24-03, IRENE Institute of Economic Research, University of Neuchâtel.

Burri, M. and D. Kaufmann (2026b). Multiple monetary policy shocks from daily data: A heteroskedasticity IV approach. IRENE Working Papers 26-06, IRENE Institute of Economic Research, University of Neuchâtel.

Jordà, Ò. (2005). Estimation and inference of impulse responses by local projections. *American Economic Review*, 95(1), 161–182.

Lewis, D. J. (2022). Robust inference in models identified via heteroskedasticity. *Review of Economics and Statistics*, 104(3), 510–524.

Lewis, D. J. and Mertens, K. (2025). A robust test for weak instruments for 2SLS with multiple endogenous regressors. *The Review of Economic Studies*, DOI: 10.1093/restud/rdaf103.

Mertens, K. and Ravn, M. O. (2013). The dynamic effects of personal and corporate income tax changes in the United States. *American Economic Review*, 103(4), 1212–1247.

Rigobon, R. (2003). Identification through heteroskedasticity. *Review of Economics and Statistics*, 85(4), 777–792.

Stock, J. H. and Watson, M. W. (2018). Identification and estimation of dynamic causal effects in macroeconomics using external instruments. *Economic Journal*, 128(610), 917–948.

## License

MIT
