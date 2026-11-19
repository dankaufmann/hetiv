# Package index

## Estimation

Estimate impulse response functions via IV local projections or a proxy
SVAR

- [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md) :
  Estimate impulse responses via heteroskedasticity-based IV local
  projections
- [`proxyiv()`](https://dankaufmann.github.io/hetiv/reference/proxyiv.md)
  : Estimate impulse responses via proxy-IV local projections
- [`recols()`](https://dankaufmann.github.io/hetiv/reference/recols.md)
  : Estimate impulse responses via recursive (Cholesky) local
  projections
- [`svariv()`](https://dankaufmann.github.io/hetiv/reference/svariv.md)
  : Estimate impulse responses via a proxy SVAR (SVAR-IV)

## Weak instrument test

Test for weak instruments (Lewis-Mertens 2025)

- [`gweakivtest()`](https://dankaufmann.github.io/hetiv/reference/gweakivtest.md)
  : Generalised weak instrument test

## Shock extraction and IRF computation

- [`predict_shocks()`](https://dankaufmann.github.io/hetiv/reference/predict_shocks.md)
  : Extract monetary policy shocks via Kalman filter
- [`prediction_corr()`](https://dankaufmann.github.io/hetiv/reference/prediction_corr.md)
  : Squared correlation of the predicted shocks with the true shocks
- [`computeirf()`](https://dankaufmann.github.io/hetiv/reference/computeirf.md)
  : Compute impulse response functions for a VAR(P)

## Proxies and model specification

Construct surprise factors and choose the lag length

- [`gss_pca()`](https://dankaufmann.github.io/hetiv/reference/gss_pca.md)
  [`print(`*`<gss_pca>`*`)`](https://dankaufmann.github.io/hetiv/reference/gss_pca.md)
  : Monetary policy surprise factors from a panel of financial market
  series
- [`lag_opt()`](https://dankaufmann.github.io/hetiv/reference/lag_opt.md)
  : Optimal lag length of a VAR

## Plotting

- [`plotirf()`](https://dankaufmann.github.io/hetiv/reference/plotirf.md)
  : Plot impulse responses from one or several estimation approaches
- [`plot2irf()`](https://dankaufmann.github.io/hetiv/reference/plot2irf.md)
  : Plot and compare impulse responses from two estimation approaches
- [`plotpval()`](https://dankaufmann.github.io/hetiv/reference/plotpval.md)
  : Plot p-values across horizons
- [`arrangeirf()`](https://dankaufmann.github.io/hetiv/reference/arrangeirf.md)
  : Arrange impulse response plots into one figure per shock
- [`plotseries()`](https://dankaufmann.github.io/hetiv/reference/plotseries.md)
  : Plot one or several series in one panel

## Simulation

- [`simulatedata()`](https://dankaufmann.github.io/hetiv/reference/simulatedata.md)
  : Simulate VAR data with heteroskedastic shocks and given parameters

## Utilities

- [`normalize()`](https://dankaufmann.github.io/hetiv/reference/normalize.md)
  : Normalize a time series
- [`logdiff()`](https://dankaufmann.github.io/hetiv/reference/logdiff.md)
  : Log-difference transformation
- [`firstdiff()`](https://dankaufmann.github.io/hetiv/reference/firstdiff.md)
  : First-difference transformation
- [`filllinear()`](https://dankaufmann.github.io/hetiv/reference/filllinear.md)
  : Linearly interpolate gaps in a data frame
- [`rescale_to_match()`](https://dankaufmann.github.io/hetiv/reference/rescale_to_match.md)
  : Rescale a series to match the mean and variance of another series
