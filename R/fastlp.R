# Internal fast estimation of one local projection by OLS or 2SLS
#
# hetiv(), proxyiv() and recols() estimate one regression per shock, variable
# and horizon and need only the coefficient on the shock variable and its
# standard error. .lp_fit() computes both by matrix algebra, which is much
# faster than fitting a model object with ivreg()/lm() and passing it to
# sandwich. The results are the same as those of
#   ivreg::ivreg() / lm()  +  sandwich::vcovHC(type = "HC3")  or
#                             sandwich::NeweyWest(prewhite = FALSE, adjust = TRUE).

# Matrix of control variables `vars` (names of columns of `data`) without the
# constant placeholder "1"; NULL if there are none.
.ctrl_matrix <- function(data, vars) {
  vars <- vars[vars != "1"]
  if (length(vars) == 0) return(NULL)
  as.matrix(data[, vars, drop = FALSE])
}

# Coefficient on x and its standard error in the regression of y on x, a
# constant and the controls W.
#
# y, x:  numeric vectors (outcome and shock variable)
# W:     matrix of controls or NULL (a constant is always added)
# z:     instrument for x, or NULL for OLS
# Wz:    controls on the instrument side (defaults to W); they differ from W
#        only under a recursive ordering, where earlier policy variables are
#        instrumented by their own instruments
# cov_type: "HC3" or "NW"
#
# Rows with missing values are dropped. Returns list(coef, se), or NULL if
# the regressors or the instruments are collinear, there are too few
# observations, or the result is not finite; the caller then falls back to
# ivreg()/lm(), which handle these cases (dropping aliased coefficients, with
# their usual warnings).
.lp_fit <- function(y, x, W = NULL, z = NULL, Wz = W, cov_type = "HC3") {
  X <- cbind("(Intercept)" = 1, shockVar = x, W)
  Z <- if (is.null(z)) NULL else cbind("(Intercept)" = 1, Z = z, Wz)

  # Complete cases across outcome, regressors and instruments
  ok <- stats::complete.cases(y, X, if (is.null(Z)) NULL else Z)
  y <- y[ok]
  X <- X[ok, , drop = FALSE]
  n <- length(y)
  k <- ncol(X)
  if (n <= k) return(NULL)

  # Projected regressors: X itself for OLS, the fitted values of the first
  # stage for 2SLS
  if (is.null(Z)) {
    Xh <- X
  } else {
    Z <- Z[ok, , drop = FALSE]
    if (ncol(Z) < k || n <= ncol(Z)) return(NULL)
    qz <- .qr_full_rank(Z)
    if (is.null(qz)) return(NULL)
    Xh <- qz$Q %*% crossprod(qz$Q, X)
    colnames(Xh) <- colnames(X)
  }
  qx <- .qr_full_rank(Xh)
  if (is.null(qx)) return(NULL)

  # Coefficients, and residuals based on the original regressors
  XtXi <- qx$XtXi
  b <- as.vector(XtXi %*% crossprod(Xh, y))
  u <- as.vector(y - X %*% b)

  if (cov_type == "NW") {
    # Newey-West through sandwich on a minimal object, so that the automatic
    # bandwidth and the small-sample adjustment are those of sandwich
    obj <- structure(
      list(estfun = Xh * u, bread = XtXi * n),
      class = "hetiv_lpfit"
    )
    # With an exact fit (e.g. the policy variable regressed on itself on
    # impact) the automatic bandwidth is undefined; fall back to the model path
    V <- tryCatch(
      sandwich::NeweyWest(obj, prewhite = FALSE, adjust = TRUE),
      error = function(e) NULL
    )
    if (is.null(V)) return(NULL)
    se <- sqrt(V[2, 2])
  } else {
    # HC3: squared residuals inflated by (1 - h)^2, with the hat values of
    # the projected regressors; only the element of the shock variable
    h <- rowSums(qx$Q^2)
    a <- as.vector(Xh %*% XtXi[, 2])
    se <- sqrt(sum((a * u / (1 - h))^2))
  }

  if (!is.finite(b[2]) || !is.finite(se)) return(NULL)
  list(coef = unname(b[2]), se = unname(se))
}

# QR decomposition of A (LAPACK, with column pivoting): the orthonormal basis
# Q and the inverse cross-product (A'A)^-1 in the original column order.
# Returns NULL if A does not have full column rank: a column whose part
# orthogonal to the other columns is below 1e-7 of its length counts as
# collinear, the tolerance of lm().
.qr_full_rank <- function(A, tol = 1e-07) {
  storage.mode(A) <- "double"
  k <- ncol(A)
  q <- qr(A, LAPACK = TRUE)
  R <- qr.R(q)
  norms <- sqrt(colSums(A^2))[q$pivot]
  d <- abs(diag(R))
  if (any(norms == 0) || any(d < tol * norms)) return(NULL)
  XtXi <- matrix(0, k, k)
  XtXi[q$pivot, q$pivot] <- chol2inv(R)
  list(Q = qr.Q(q), XtXi = XtXi)
}

# Methods used by sandwich::NeweyWest() on the minimal object of .lp_fit():
# estimating functions (projected regressors times residuals) and bread
# (n times the inverse cross-product), as for ivreg and lm objects

#' @importFrom sandwich estfun
#' @export
estfun.hetiv_lpfit <- function(x, ...) x$estfun

#' @importFrom sandwich bread
#' @export
bread.hetiv_lpfit <- function(x, ...) x$bread
