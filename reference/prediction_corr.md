# Squared correlation of the predicted shocks with the true shocks

Computes the population squared correlation between the shocks predicted
by
[`predict_shocks()`](https://dankaufmann.github.io/hetiv/reference/predict_shocks.md)
and the true underlying shocks, as a measure of the precision of the
prediction (Burri and Kaufmann, 2026a).

## Usage

``` r
prediction_corr(Psi, SigR, Sig = NULL, squared = TRUE)
```

## Arguments

- Psi:

  Numeric matrix (M x E) of responses of the residuals to the E
  structural shocks (e.g. `Psi` from
  [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md)).

- SigR:

  Numeric matrix (M x M). Covariance matrix of the residuals on control
  (non-event) days.

- Sig:

  Numeric matrix (M x M) or `NULL`. Covariance matrix of the residuals
  on event days, used to rescale `Psi` to one-standard-deviation shocks.
  If `NULL` (default), `Psi` is assumed to be in these units already.

- squared:

  Logical. If `TRUE` (default), the squared correlation is returned. If
  `FALSE`, its square root (the correlation).

## Value

A numeric vector of length E with the squared correlation between each
predicted shock and the true shock (the diagonal of `S (I + S)^{-1}`),
or the correlation if `squared = FALSE`.

## Details

Let the residuals on event days be `e_t = Psi eps_t + v_t`, with
unit-variance shocks `eps_t` and non-policy noise `v_t` with covariance
`SigR`. The linear prediction of the shocks has covariance with the true
shocks `S (I + S)^{-1}`, with the signal-to-noise matrix
`S = Psi' SigR^{-1} Psi`. The diagonal of this matrix is the share of
the shock variance explained by the prediction (one minus the relative
mean squared error), that is, the squared correlation between the
predicted and the true shock.

The correlation itself does not depend on the units of the shocks, but
the formula does: `S` is the signal-to-noise ratio only if `Psi` is the
response to a one-standard-deviation shock. The estimators return `Psi`
normalised to a unit impact response instead, which leaves the variance
of the shocks out of `S`. If `Sig` is supplied, the columns of `Psi` are
therefore rescaled by the implied standard deviations of the shocks,
recovered from `Sig - SigR` in the same way as in
[`predict_shocks()`](https://dankaufmann.github.io/hetiv/reference/predict_shocks.md).
If `Sig` is `NULL`, `Psi` is used as is.

The measure is a population quantity evaluated at the estimated
parameters. It does not account for estimation error in `Psi`, `Sig` and
`SigR`, and it is informative only if the responses in `Psi` are
consistently estimated, that is, if the identifying assumptions of the
estimator hold.

The measure also requires a nonsingular covariance `SigR`. If fewer
non-policy shocks than variables drive the data, `SigR` is (nearly)
singular, the implied signal-to-noise ratio becomes very large and the
measure approaches one whatever the actual precision. A warning is
issued in this case.

## References

Burri, M. and Kaufmann, D. (2026a). Measuring monetary policy shocks.
IRENE Working Papers 24-03, IRENE Institute of Economic Research,
University of Neuchatel.

## See also

[`predict_shocks()`](https://dankaufmann.github.io/hetiv/reference/predict_shocks.md)

## Examples

``` r
# One shock with unit variance, two variables
Psi <- matrix(c(1, 0.5), nrow = 2)
SigR <- diag(c(0.5, 0.5))
prediction_corr(Psi, SigR)
#> [1] 0.7142857

# Psi normalised to unit impact; the shock has variance 4
Sig <- SigR + 4 * Psi %*% t(Psi)
prediction_corr(Psi, SigR, Sig = Sig)
#> [1] 0.9090909
```
