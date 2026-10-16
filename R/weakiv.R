# Internal first-stage diagnostics shared by svariv(), hetiv() and proxyiv().

# One row of the WeakTest data frame returned by hetiv() and proxyiv(): the
# Lewis-Mertens (2025) test for outcome y (horizon h, variable i) in the system
# with endogenous regressors Y (T x e), controls X (or NULL) and instruments Z
# (T x e). cov_type "NW" maps to the HAC version, anything else to EHW.
# target ("beta" or an integer <= e) and crit ("abs" or "rel") are passed to
# gweakivtest(); hetiv() and proxyiv() use cov_type = "NW", crit = "rel" and
# target = e.
.lm_weaktest_row <- function(y, Y, X, Z, h, i, cov_type = "EHW", target = "beta", crit = "abs") {
  if (is.null(X)) X <- matrix(1, nrow = length(y), ncol = 1)
  wt <- gweakivtest(
    y = y, Y = Y, X = X, Z = Z, target = target, crit = crit,
    cov_type = if (cov_type == "NW") "NW" else "EHW"
  )
  data.frame(
    h = h, i = i, target = target, Tobs = wt$nobs,
    gmin = wt$gmin_generalized,
    cv = wt$gmin_generalized_critical_value,
    cv_simplified = wt$gmin_generalized_critical_value_simplified,
    sy_stat = wt$stock_yogo_test_statistic,
    sy_cv = wt$stock_yogo_critical_value_nagar
  )
}

# First-stage F statistics for one endogenous regressor y1 instrumented by Z
# (T x K), with optional exogenous controls X (a constant is always added).
#
# Returns a list with
#   F_hom: homoskedastic F statistic of the K instruments (Wald statistic
#          divided by K); for K = 1 this is the squared t statistic of the
#          instrument (Stock and Yogo, 2005)
#   F_rob: heteroskedasticity- (or HAC-) robust effective F statistic of
#          Montiel Olea and Pflueger (2013),
#          pi' Q_ZZ pi / tr(Q_ZZ V_rob), where Q_ZZ is the second-moment
#          matrix of the instruments after partialling out the controls and
#          V_rob the robust covariance of the first-stage coefficients; for
#          K = 1 this is the squared robust t statistic
#   R2:    R-squared of the first-stage regression
#   Tobs:  number of complete observations used
#   K:     number of instruments
#
# cov_type selects the robust covariance: "HC0", "HC1", "HC3" (sandwich::vcovHC)
# or "NW" (sandwich::NeweyWest, no prewhitening). Rows with missing values
# are dropped. Collinear controls are removed before computing the statistics;
# if the instruments themselves are collinear, NA is returned.
.firststage_f <- function(y1, Z, X = NULL, cov_type = "HC0") {
  cov_type <- .check_choice(cov_type, "cov_type", c("HC0", "HC1", "HC3", "NW"))
  Z <- as.matrix(Z)
  K <- ncol(Z)
  if (!is.null(X)) X <- as.matrix(X)

  # Complete cases across all inputs
  ok <- stats::complete.cases(y1, Z, if (is.null(X)) NULL else X)
  y1 <- y1[ok]
  Z <- Z[ok, , drop = FALSE]
  if (!is.null(X)) X <- X[ok, , drop = FALSE]
  Tobs <- length(y1)
  idx <- 1 + seq_len(K) # positions of the instrument coefficients

  na_result <- list(
    F_hom = NA_real_, F_rob = NA_real_, R2 = NA_real_, Tobs = Tobs, K = K
  )
  if (Tobs <= K + 1 + if (is.null(X)) 0 else ncol(X)) return(na_result)

  fit <- if (is.null(X)) stats::lm(y1 ~ Z) else stats::lm(y1 ~ Z + X)
  cf <- stats::coef(fit)
  if (anyNA(cf[idx])) return(na_result)
  # Drop aliased (collinear) controls and refit, so that vcov() is complete
  if (!is.null(X) && anyNA(cf[-seq_len(K + 1)])) {
    X <- X[, !is.na(cf[-seq_len(K + 1)]), drop = FALSE]
    fit <- if (ncol(X) == 0) stats::lm(y1 ~ Z) else stats::lm(y1 ~ Z + X)
    if (ncol(X) == 0) X <- NULL
    cf <- stats::coef(fit)
  }
  b <- cf[idx]

  # Homoskedastic F (Wald / K)
  V_hom <- stats::vcov(fit)[idx, idx, drop = FALSE]
  F_hom <- drop(crossprod(b, solve(V_hom, b))) / K

  # Robust covariance of the instrument coefficients
  V_rob <- if (cov_type == "NW") {
    sandwich::NeweyWest(fit, prewhite = FALSE, adjust = TRUE)
  } else {
    sandwich::vcovHC(fit, type = cov_type)
  }
  V_rob <- V_rob[idx, idx, drop = FALSE]

  # Q_ZZ of the instruments after partialling out the constant and controls
  Zt <- if (is.null(X)) {
    sweep(Z, 2, colMeans(Z))
  } else {
    stats::resid(stats::lm(Z ~ X))
  }
  Zt <- matrix(Zt, nrow = Tobs, ncol = K)
  Qzz <- crossprod(Zt) / Tobs
  F_rob <- drop(crossprod(b, Qzz %*% b)) / sum(diag(Qzz %*% V_rob))

  list(
    F_hom = unname(F_hom), F_rob = unname(F_rob),
    R2 = summary(fit)$r.squared, Tobs = Tobs, K = K
  )
}
