# Estimate impulse responses via heteroskedasticity-based IV local projections

Estimates impulse response functions (IRFs) using recursive
heteroskedasticity-IV identification (Rigobon, 2003; Rigobon and Sack,
2004; Lewis, 2022; Burri and Kaufmann, 2026) combined with local
projections (Jordà, 2005). Identification exploits the difference in
variance between policy event days and control days to construct
instruments for the endogenous variables.

## Usage

``` r
hetiv(
  y,
  O,
  X = NULL,
  Ind,
  P,
  H,
  E = 1,
  norm = 1,
  interact = FALSE,
  cum = FALSE,
  Hstep = 1,
  details = FALSE
)
```

## Arguments

- y:

  Numeric matrix of stationary outcome variables (T x N). The effect on
  the first variable in each dimension is normalized to unity at
  horizon 0. These variables are also used to construct
  heteroskedasticity-based instruments and for recursive ordering to
  identify multiple dimensions.

- O:

  Numeric matrix of information set variables (T x M). May be identical
  to `y`. Included as lags 1 through `P`.

- X:

  Numeric matrix of deterministic variables (T x K). May include a
  constant, time trend, seasonal dummies or other deterministic
  controls. Included as is (no lags).

- Ind:

  Integer vector of length T, event indicator:

  - `0` Control day (no event)

  - `1` Policy day (event)

  - `2` Contaminated control day (excluded from estimation)

- P:

  Integer. Maximum lag order for the information set. Set to `0` for no
  lags (regression on constant only).

- H:

  Integer. Maximum horizon (in periods) up to which IRFs are estimated.

- E:

  Integer. Number of shock dimensions to identify via recursive
  ordering.

- norm:

  Numeric scalar. Normalize the impact response of the first variable to
  a specific value. Set to `1` for standard unit-effect normalization.

- interact:

  Logical. If `TRUE`, lagged information set variables are interacted
  with event/non-event dummies.

- cum:

  Logical vector of length N. For each variable in `y`, whether to
  report the cumulative impulse response instead of the level response.
  If only one provided, applied to all impulse responses.

- Hstep:

  Integer. Step size between horizons. The default `1` estimates all
  horizons 0 through H - 1. Values greater than 1 are intended only for
  fast testing; they are only safe when `Hstep >= H` (a single horizon
  is stored). For complete IRF estimation always use `Hstep = 1`.

- details:

  Logical. If `TRUE`, code saves detailed IV results, which is slower.
  if set to `FALSE`, returns only impulse response and standard error
  (e.g. for bootstrap)

## Value

A named list with the following elements:

- `irf`:

  Array (H x N x E) of estimated impulse responses.

- `se`:

  Array (H x N x E) of HC0 heteroscedasticity-robust standard errors.

- `IVRes`:

  List of `ivreg` model objects, one per horizon, variable, and shock
  dimension.

- `Obs`:

  Data frame with observation counts: `Tp` (policy days), `Tc` (control
  days), `To` (contaminated days), `Tt` (total used).

- `Method`:

  Character string `"Heteroscedasticity-IV"`.

- `et`:

  Data frame of OLS residuals on event days (used for covariance
  estimation and shock extraction).

- `Sig`:

  Covariance matrix of residuals on event days (used for shock
  extraction).

- `SigR`:

  Covariance matrix of residuals on control days, or `NA` if unavailable
  (used for shock extraction).

- `Psi`:

  Impact matrix (N x E), equal to `irf[1, , ]`. By the package's
  indexing convention `HSeries` starts at 1, so the first LP uses
  `lead(y, 0)` (the contemporaneous value) and is labelled horizon 0;
  `irf[1, , ]` is therefore always the impact response.

- `WeakData`:

  Data frame of endogenous variables and instruments for the
  Lewis-Mertens (2025) weak instrument test.

## References

Burri, M. and D. Kaufmann (2026). Measuring monetary policy shocks.
IRENE Working Papers 24-03, IRENE Institute of Economic Research,
University of Neuchâtel.

Burri, M. and D. Kaufmann (2026). Multiple monetary policy shocks from
daily data: A heteroskedasticity IV approach. IRENE Working Papers
26-06, IRENE Institute of Economic Research, University of Neuchâtel.

Jordà, Ò. (2005). Estimation and Inference of Impulse Responses by Local
Projections. *American Economic Review*, 95(1), 161–182.

Lewis, D. J. (2022). Robust Inference in Models Identified via
Heteroskedasticity. *Review of Economics and Statistics*, 104(3),
510–524.

Lewis, D. J. and Mertens, K. (2025). A Robust Test for Weak Instruments
for 2SLS with Multiple Endogenous Regressors. *The Review of Economic
Studies*, DOI: 10.1093/restud/rdaf103

Rigobon, R. (2003). Identification Through Heteroskedasticity. *Review
of Economics and Statistics*, 85(4), 777–792.

Rigobon, R. and Sack, B. (2004). The impact of monetary policy on asset
prices. *Journal of Monetary Economics*, 51(8), 1553–1575.
