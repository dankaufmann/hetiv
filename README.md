# hetiv

An R package for measuring and identifying structural shocks using heteroskedasticity- and proxy-based instrumental variable (IV) methods with daily financial market data.

**Authors**: Daniel Kaufmann, Marc Burri, Valentin Grob

**Note**: This is work in progress. Installation and use at own risk.

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
| `gweaktest()` | Tests for weak instruments using the generalised minimum eigenvalue statistic of Lewis and Mertens (2025) |
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

# y:   T x N matrix of outcome variables
# O:   T x M matrix of information set variables
# X:   T x K matrix of deterministic variables (optional; e.g. constant, trend, dummies)
# Ind: event indicator (0 = control day, 1 = policy day, 2 = contaminated)

res <- hetiv(y = y, O = O, Ind = Ind, P = 1, H = 20, E = 1, norm = 1, details = TRUE)

# With a linear time trend
X <- matrix(seq_len(nrow(y)), ncol = 1)
res <- hetiv(y = y, O = O, X = X, Ind = Ind, P = 1, H = 20, E = 1, norm = 1, details = TRUE)
```

### Estimate impulse responses via proxy-IV

```r
# Z: T x E matrix of external instruments (proxies)

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
# y: T x 1 regressand
# Y: T x N matrix of endogenous regressors
# X: T x Nx matrix of exogenous controls
# Z: T x K matrix of instruments (requires K >= N)

wt <- gweaktest(y = y, Y = Y, X = X, Z = Z)

# Inspect test statistics and critical values
wt$gmin_generalized                            # generalised min-eigenvalue statistic
wt$gmin_generalized_critical_value             # Lewis-Mertens sharp critical value (Stiefel)
wt$gmin_generalized_critical_value_simplified  # conservative closed-form bound
wt$stock_yogo_test_statistic                   # Stock-Yogo statistic (Nagar approximation)
wt$stock_yogo_critical_value_nagar             # Stock-Yogo critical value

# With Newey-West HAR standard errors and custom bias tolerance
wt <- gweaktest(y = y, Y = Y, X = X, Z = Z, cov_type = "NW", tau = 0.10)

# hetiv() and proxyiv() return WeakData when details = TRUE,
# which contains the outcome and instrument columns needed for gweaktest()
res     <- hetiv(y = y, O = O, Ind = Ind, P = 1, H = 20, E = 1, details = TRUE)
z_cols  <- grep("^Z", names(res$WeakData))
y_cols  <- grep("^y", names(res$WeakData))
wt      <- gweaktest(y = res$WeakData[, y_cols[1], drop = FALSE],
                     Y = res$WeakData[, y_cols[-1], drop = FALSE],
                     X = matrix(numeric(0), nrow(res$WeakData), 0),
                     Z = res$WeakData[, z_cols, drop = FALSE])
```

### Extract structural shocks

```r
shocks <- kfpredict(Sig = res$Sig, SigR = res$SigR,
                    Psi = res$Psi, et = res$et, tol = 1e-10)
```

### Simulate data

```r
sim <- simulatedata(Phi = Phi, SigE = 4, PsiE = PsiE, PsiR = PsiR,
                    Nobs = 500, Nbin = 100, N = 2, R = 2, E = 1,
                    Nevn = 5, P = 1, eDist = 0, seed = 42)
```

## References

Burri, M. and D. Kaufmann (2026). Measuring monetary policy shocks. IRENE Working Papers 24-03, IRENE Institute of Economic Research, University of Neuchâtel. 

Burri, M. and D. Kaufmann (2026). Multiple monetary policy shocks from daily data: A heteroskedasticity IV approach.  IRENE Working Papers 26-06, IRENE Institute of Economic Research, University of Neuchâtel. 

Jordà, Ò. (2005). Estimation and Inference of Impulse Responses by Local Projections. *American Economic Review*, 95(1), 161–182.

Lewis, D. J. (2022). Robust Inference in Models Identified via Heteroskedasticity. *Review of Economics and Statistics*, 104(3), 510–524.

Lewis, D. J. and Mertens, K. (2025). A Robust Test for Weak Instruments for 2SLS with Multiple Endogenous Regressors. *The Review of Economic Studies*, DOI: 10.1093/restud/rdaf103.

Mertens, K. and Ravn, M. O. (2013). The Dynamic Effects of Personal and Corporate Income Tax Changes in the United States. *American Economic Review*, 103(4), 1212–1247.

Rigobon, R. (2003). Identification Through Heteroskedasticity. *Review of Economics and Statistics*, 85(4), 777–792.

Stock, J. H. and Watson, M. W. (2018). Identification and Estimation of Dynamic Causal Effects in Macroeconomics Using External Instruments. *Economic Journal*, 128(610), 917–948.

## License

MIT
