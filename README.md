# hetpxyiv

An R package for measuring and identifying structural shocks using heteroskedasticity- and proxy-based instrumental variable (IV) methods with daily financial market data.

**Authors**: Daniel Kaufmann, Marc Burri, Valentin Grob

**Note**: This is work in progress. Installation and use at own risk.

## Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("dankaufmann/hetpxyiv")
```

## Functions

| Function | Description |
|---|---|
| `hetiv()` | Estimates impulse response functions via heteroskedasticity-based IV local projections |
| `kfpredict()` | Extracts structural shocks from reduced-form residuals via Kalman filter |
| `normalize()` | Standardizes a time series to zero mean and unit standard deviation |
| `plotirf()` | Plots IRFs with confidence bands for a single estimation approach |
| `plot2irf()` | Plots and compares IRFs from two estimation approaches |
| `simulatedata()` | Simulates VAR data with heteroskedastic event shocks |

## Usage

### Estimate impulse responses

```r
library(hetpxyiv)

# y:   T x N matrix of outcome variables
# O:   T x M matrix of information set variables
# X:   T x K matrix of deterministic variables (optional; e.g. constant, trend, dummies)
# Ind: event indicator (0 = control day, 1 = policy day, 2 = contaminated)

# Without deterministic terms
res <- hetiv(y = y, O = O, Ind = Ind, P = 1, H = 20, E = 1, norm = 1, details = TRUE)

# With a linear time trend as deterministic term
X <- matrix(seq_len(nrow(y)), ncol = 1)
res <- hetiv(y = y, O = O, X = X, Ind = Ind, P = 1, H = 20, E = 1, norm = 1, details = TRUE)
```

### Plot impulse responses

```r
plots <- plotirf(IRFest = res$irf, IRFse = res$se,
                 HTick = 5, Labels = c("Var 1", "Var 2"))
cowplot::plot_grid(plotlist = plots)
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

Rigobon, R. (2003). Identification Through Heteroskedasticity. *Review of Economics and Statistics*, 85(4), 777–792.

## License

MIT
