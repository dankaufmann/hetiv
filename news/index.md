# Changelog

## hetiv 0.1.2

- New
  [`recols()`](https://dankaufmann.github.io/hetiv/reference/recols.md):
  recursive (Cholesky) local projections estimated by OLS, with the same
  interface and output as
  [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md),
  as a benchmark identification.
- New
  [`arrangeirf()`](https://dankaufmann.github.io/hetiv/reference/arrangeirf.md):
  arranges IRF panels from
  [`plotirf()`](https://dankaufmann.github.io/hetiv/reference/plotirf.md)/[`plot2irf()`](https://dankaufmann.github.io/hetiv/reference/plot2irf.md)
  into one figure per shock and optionally saves them as PDF. Adds
  `gridExtra` to Imports.
- Various smaller improvements and additional options in
  [`hetiv()`](https://dankaufmann.github.io/hetiv/reference/hetiv.md)
  and
  [`proxyiv()`](https://dankaufmann.github.io/hetiv/reference/proxyiv.md)

## hetiv 0.1.1

- Changes to default HC settings for local projections

## hetiv 0.1.0

- Initial CRAN-ready release.
- Provides heteroskedasticity-based and proxy-based instrumental
  variable local projection estimators for event-study settings.
- Includes impulse response plotting helpers, weak instrument testing,
  Kalman-filter shock extraction, and simulation utilities.
- Adds an introductory vignette with a worked end-to-end example.
