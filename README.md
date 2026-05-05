# hetpxyiv

An R package with tools for measuring and identifying monetary policy shocks using heteroskedasticity- and proxy-based instrumental variable (IV) methods with daily financial market data.

**Authors**: Daniel Kaufmann, Marc Burri, Valentin Grob

**Notes**: This is work in progress. Installation and use at own risk

## Installation

You can install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("dankaufmann/hetpxyiv")
```

## Functions

| Function | Description |
|---|---|
| `normalize()` | Standardizes a time series to zero mean and unit standard deviation |
| `estimLPHet()` | Estimates impulse response functions via heteroscedasticity-based IV local projections |

## Usage

```r
library(hetpxyiv)

# Normalize a time series
x <- c(1.2, 0.8, 1.5, 0.9, 1.1)
normalize(x)
```

## License

MIT
