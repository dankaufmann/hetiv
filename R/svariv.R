#' Estimate impulse responses via a proxy SVAR (SVAR-IV)
#'
#' Estimates a reduced-form VAR(P) by OLS, identifies the impact matrix of the
#' first `E` structural shocks with external instruments (proxies) following
#' Mertens and Ravn (2013) and Stock and Watson (2018), and computes impulse
#' responses by recursion through the VAR. Optionally computes confidence
#' bands with a wild bootstrap following Mertens and Ravn (2013).
#'
#' @details The first `E` columns of `y` are the policy variables: shock `e`
#'   is the structural shock to variable `e`, and the impact response of
#'   variable `e` to shock `e` is normalized to `norm`. For `E > 1` the
#'   impact matrix of the policy block is identified up to a lower-triangular
#'   (Cholesky) rotation, so shock `e` has no impact effect on policy variables
#'   `1, ..., e-1`; identification is therefore order-dependent, as in
#'   [proxyiv()] with `recursive = TRUE`.
#'
#'   All columns of `Z` are used jointly as instruments in the first-stage
#'   regressions of the policy residuals on the proxies. Rows of `Z` with
#'   missing values are excluded from the identification step only: the VAR
#'   and the reduced-form residual covariance `Sigma` are estimated on the
#'   full sample, as in Mertens and Ravn (2013), whereas the 2SLS step that
#'   recovers the relative impact responses uses only the observations with
#'   non-missing proxies (the identification sample).
#'
#'   With `details = TRUE`, first-stage diagnostics are reported for every
#'   policy residual: the homoskedastic F statistic of the proxies (`FStat`,
#'   equal to the squared t statistic when `M = 1`; Stock and Yogo, 2005) and
#'   the robust effective F statistic of Montiel Olea and Pflueger (2013)
#'   (`FStatRobust`), computed with the covariance estimator in `cov_type`.
#'   For `E > 1` these are equation-by-equation statistics; use
#'   [gweakivtest()] for a joint test.
#'
#'   The bootstrap follows Mertens and Ravn (2013): in every replication, the
#'   VAR residuals and the proxies are multiplied by the same Rademacher draw,
#'   artificial data are generated recursively from the first `P`
#'   observations, and the VAR is re-estimated and re-identified. Each draw is
#'   normalized in the same way as the point estimate. Percentile bands are
#'   reported for every level in `ci`; `se` is the bootstrap standard
#'   deviation and `boot` holds the raw draws, so the result can be passed
#'   directly to [plotirf()], which draws recentred percentile bands from the
#'   draws. Replications in which identification fails (e.g. a
#'   non-positive-definite policy block) are dropped with a warning.
#'
#' @param y Numeric matrix of VAR variables (T x N) without missing values.
#'   The first `E` columns are the policy variables on which the responses
#'   are normalized.
#' @param Z Numeric matrix of external instruments (T x M, `M >= E`). Missing
#'   values are allowed and drop the corresponding observation from the
#'   identification step. Set the proxy to zero (not `NA`) on days without an
#'   event.
#' @param X Numeric matrix of deterministic variables (T x K), or `NULL`
#'   (default). Included as-is (no lags) in addition to the constant, which is
#'   always included.
#' @param P Integer. Lag order of the VAR (`P >= 1`).
#' @param H Integer. Number of horizons: IRFs are reported for horizons
#'   `0, ..., H - 1`.
#' @param E Integer. Number of shocks to identify (`E < N`). Default `1`.
#' @param norm Numeric scalar. Normalize the impact response of variable `e` to
#'   shock `e` to this value. Default `1`.
#' @param cum Logical scalar or vector of length N. If `TRUE` for variable `i`,
#'   the cumulative IRF is reported. A single value is recycled to all
#'   variables. Default `FALSE`.
#' @param nboot Integer. Number of wild-bootstrap replications. `0` (default)
#'   skips the bootstrap, in which case `se`, `lower`, and `upper` are `NULL`.
#' @param ci Numeric vector of confidence levels in (0, 1) for the percentile
#'   bands. Default `c(0.90, 0.95)`.
#' @param cov_type Covariance estimator for the robust first-stage F statistic
#'   reported with `details = TRUE`: `"HC0"` (default, Eicker-Huber-White),
#'   `"HC1"`, `"HC3"`, or `"NW"` for Newey-West HAC.
#' @param details Logical. If `TRUE`, additionally returns the VAR estimates,
#'   residuals and first-stage statistics.
#'
#' @return A named list. Always contains:
#'   \describe{
#'     \item{`irf`}{Array (H x N x E) of impulse responses, horizons labelled
#'       `0` to `H - 1`.}
#'     \item{`se`}{Array (H x N x E) of bootstrap standard deviations, or
#'       `NULL` if `nboot = 0`.}
#'     \item{`lower`, `upper`}{Arrays (H x N x E x length(ci)) of percentile
#'       bands, one slice per level in `ci` (labelled in the fourth
#'       dimension), or `NULL` if `nboot = 0`.}
#'     \item{`boot`}{Array (H x N x E x nboot) of bootstrap draws (`NA` for
#'       failed replications), or `NULL` if `nboot = 0`.}
#'     \item{`Method`}{Character string `"SVAR-IV"`.}
#'   }
#'   With `details = TRUE`, additionally contains:
#'   \describe{
#'     \item{`Psi`}{Normalized impact matrix (N x E), equal to `irf[1, , ]`
#'       before cumulation.}
#'     \item{`Phi`}{Array (N x N x P) of VAR slope coefficients.}
#'     \item{`C`}{Vector of length N with the VAR intercepts.}
#'     \item{`D`}{Matrix (N x K) of coefficients on `X`, or `NULL`.}
#'     \item{`U`}{Matrix ((T - P) x N) of reduced-form residuals.}
#'     \item{`FStat`, `FStatRobust`, `R2`}{Vectors of length E with the
#'       homoskedastic first-stage F statistic, the robust effective F
#'       statistic (Montiel Olea and Pflueger, 2013) and the R-squared of the
#'       regressions of each policy residual on the proxies (see Details).}
#'     \item{`Obs`}{Data frame with `T` (sample size), `Test` (VAR estimation
#'       sample), and `Tiv` (identification sample).}
#'     \item{`nfail`}{Number of failed bootstrap replications.}
#'     \item{`ci`}{The confidence levels used for `lower` and `upper`.}
#'   }
#'
#' @references
#' Mertens, K. and Ravn, M. O. (2013). The dynamic effects of personal and
#' corporate income tax changes in the United States. *American Economic
#' Review*, 103(4), 1212-1247.
#'
#' Montiel Olea, J. L. and Pflueger, C. (2013). A robust test for weak
#' instruments. *Journal of Business & Economic Statistics*, 31(3), 358-369.
#'
#' Stock, J. H. and Watson, M. W. (2018). Identification and estimation of
#' dynamic causal effects in macroeconomics using external instruments.
#' *Economic Journal*, 128(610), 917-948.
#'
#' Stock, J. H. and Yogo, M. (2005). Testing for weak instruments in linear IV
#' regression. In D. W. K. Andrews and J. H. Stock (Eds.), *Identification
#' and inference for econometric models: Essays in honor of Thomas
#' Rothenberg* (pp. 80-108). Cambridge University Press.
#'
#' @seealso [proxyiv()] for the local-projection counterpart, [computeirf()],
#'   [plotirf()].
#'
#' @importFrom stats complete.cases lm.fit runif sd quantile
#'
#' @examples
#' set.seed(1)
#' Phi <- array(c(0.5, 0.1, -0.1, 0.4), dim = c(2, 2, 1))
#' Psi <- matrix(c(1, 0.5, 0, 1), 2, 2)
#' Nobs <- 400
#' eps <- matrix(rnorm(2 * Nobs), ncol = 2)
#' y <- matrix(0, Nobs, 2)
#' for (t in 2:Nobs) y[t, ] <- Phi[, , 1] %*% y[t - 1, ] + Psi %*% eps[t, ]
#' Z <- matrix(eps[, 1] + rnorm(Nobs), ncol = 1)
#' res <- svariv(y = y, Z = Z, P = 1, H = 6, nboot = 50, ci = 0.90)
#' dim(res$irf)
#' res$irf[1, , 1] # impact responses, normalized to 1 on the first variable
#'
#' @export
svariv <- function(y, Z, X = NULL, P, H, E = 1, norm = 1, cum = FALSE,
                   nboot = 0, ci = c(0.90, 0.95), cov_type = "HC0",
                   details = FALSE) {
  y <- .as_numeric_matrix(y, "y", allow_na = FALSE)
  Z <- .as_numeric_matrix(Z, "Z", nrow = nrow(y))
  if (!is.null(X)) {
    X <- .as_numeric_matrix(X, "X", nrow = nrow(y), allow_na = FALSE)
  }
  P <- .check_integerish_scalar(P, "P", min = 1)
  H <- .check_integerish_scalar(H, "H", min = 1)
  E <- .check_integerish_scalar(E, "E", min = 1)
  norm <- .check_numeric_scalar(norm, "norm")
  cum <- .check_cum(cum, ncol(y))
  nboot <- .check_integerish_scalar(nboot, "nboot", min = 0)
  cov_type <- .check_choice(cov_type, "cov_type", c("HC0", "HC1", "HC3", "NW"))
  details <- .check_logical_scalar(details, "details")
  if (
    !is.numeric(ci) || length(ci) < 1 || anyNA(ci) ||
      any(ci <= 0 | ci >= 1)
  ) {
    stop("ci must contain confidence levels between 0 and 1.", call. = FALSE)
  }

  Nobs <- nrow(y)
  N <- ncol(y)
  M <- ncol(Z)
  K <- if (!is.null(X)) ncol(X) else 0

  if (E >= N) {
    stop("E must be smaller than the number of columns in y.", call. = FALSE)
  }
  if (M < E) {
    stop("Z must have at least E columns (one instrument per shock).",
      call. = FALSE
    )
  }
  if (Nobs <= P + N * P + K + 1) {
    stop("Not enough observations to estimate the VAR.", call. = FALSE)
  }

  # 1) Reduced-form VAR(P) by OLS (constant always included)
  fit <- .svariv_var(y, X, P)

  # 2) Identification: Sigma on the full VAR sample, 2SLS on the observations
  #    with a non-missing proxy
  Zs <- Z[(P + 1):Nobs, , drop = FALSE]
  ident <- stats::complete.cases(Zs)
  Tiv <- sum(ident)
  if (Tiv <= M + 1) {
    stop("Not enough non-missing proxy observations for identification.",
      call. = FALSE
    )
  }
  id <- .svariv_identify(fit$U, Zs, ident, E, fstat = details, cov_type = cov_type)

  # 3) Impulse responses by recursion, normalized on the policy variables
  irf <- .svariv_irf(id$Psi, fit$Phi, H, E, norm, cum)
  Psi <- id$Psi
  for (e in seq_len(E)) {
    Psi[, e] <- Psi[, e] / Psi[e, e] * norm
  }

  # 4) Wild bootstrap (Rademacher weights on residuals and proxies jointly)
  se <- lower <- upper <- boot <- NULL
  nfail <- 0L
  if (nboot > 0) {
    boot <- array(NA_real_, dim = c(H, N, E, nboot))
    Teff <- Nobs - P
    for (j in seq_len(nboot)) {
      rr <- 1 - 2 * (stats::runif(Teff) > 0.5)
      bootU <- fit$U * rr
      bootZ <- Zs * rr
      yb <- .svariv_simulate(y, X, fit, bootU, P)
      fb <- .svariv_var(yb, X, P)
      draw <- tryCatch(
        {
          idb <- .svariv_identify(fb$U, bootZ, ident, E)
          .svariv_irf(idb$Psi, fb$Phi, H, E, norm, cum)
        },
        error = function(e) NULL
      )
      if (is.null(draw)) {
        nfail <- nfail + 1L
      } else {
        boot[, , , j] <- draw
      }
    }
    if (nfail > 0) {
      warning(nfail, " of ", nboot, " bootstrap replications failed and were dropped.",
        call. = FALSE
      )
    }
    if (nfail == nboot) {
      stop("All bootstrap replications failed.", call. = FALSE)
    }

    se <- apply(boot, c(1, 2, 3), stats::sd, na.rm = TRUE)
    lower <- array(NA_real_, dim = c(H, N, E, length(ci)))
    upper <- array(NA_real_, dim = c(H, N, E, length(ci)))
    for (k in seq_along(ci)) {
      lower[, , , k] <- apply(boot, c(1, 2, 3), stats::quantile,
        probs = (1 - ci[k]) / 2, na.rm = TRUE, names = FALSE
      )
      upper[, , , k] <- apply(boot, c(1, 2, 3), stats::quantile,
        probs = (1 + ci[k]) / 2, na.rm = TRUE, names = FALSE
      )
    }
    dimnames(se) <- dimnames(irf)
    dimnames(lower) <- dimnames(upper) <- c(dimnames(irf), list(as.character(ci)))
    dimnames(boot) <- c(dimnames(irf), list(NULL))
  }

  Method <- "SVAR-IV"

  if (details == TRUE) {
    Obs <- data.frame(T = Nobs, Test = Nobs - P, Tiv = Tiv)
    return(list(
      irf = irf, se = se, lower = lower, upper = upper, boot = boot,
      Method = Method,
      Psi = Psi, Phi = fit$Phi, C = fit$C, D = fit$D, U = fit$U,
      FStat = id$FStat, FStatRobust = id$FStatRobust, R2 = id$R2, Obs = Obs,
      nfail = nfail, ci = ci
    ))
  } else {
    # boot is returned also without details so that plotirf() can draw
    # percentile bands from the draws
    return(list(
      irf = irf, se = se, lower = lower, upper = upper, boot = boot,
      Method = Method
    ))
  }
}

# Reduced-form VAR(P) with constant and optional deterministics by OLS.
# Returns intercepts C (N), slopes Phi (N x N x P, Phi[, , p] is the matrix on
# lag p), coefficients on X (N x K or NULL), and residuals U ((T - P) x N).
.svariv_var <- function(y, X, P) {
  Nobs <- nrow(y)
  N <- ncol(y)
  Ylag <- do.call(cbind, lapply(seq_len(P), function(p) {
    y[(P + 1 - p):(Nobs - p), , drop = FALSE]
  }))
  W <- cbind(1, Ylag)
  if (!is.null(X)) {
    W <- cbind(W, X[(P + 1):Nobs, , drop = FALSE])
  }
  Ydep <- y[(P + 1):Nobs, , drop = FALSE]
  coef <- stats::lm.fit(W, Ydep)$coefficients
  coef <- matrix(coef, nrow = ncol(W), ncol = N)
  if (anyNA(coef)) {
    stop("VAR regressors are collinear; check y, X, and P.", call. = FALSE)
  }
  U <- Ydep - W %*% coef

  Phi <- array(0, dim = c(N, N, P))
  for (p in seq_len(P)) {
    Phi[, , p] <- t(coef[1 + ((p - 1) * N + 1):(p * N), , drop = FALSE])
  }
  D <- if (!is.null(X)) t(coef[(2 + N * P):nrow(coef), , drop = FALSE]) else NULL

  list(C = coef[1, ], Phi = Phi, D = D, U = U)
}

# Identify the impact matrix of the first E shocks with external instruments
# (Mertens and Ravn, 2013). U: reduced-form residuals of the full VAR sample
# (T x N), the first E columns being the policy variables; Z: proxies (T x M,
# may contain NA); ident: logical vector marking the identification sample
# (rows with non-missing proxies). The residual covariance Sigma is computed
# on the full sample, the 2SLS step on the identification sample only.
# Returns Psi (N x E, policy block lower triangular with positive diagonal)
# and, if fstat = TRUE, the first-stage statistics of each policy residual
# (see .firststage_f()).
.svariv_identify <- function(U, Z, ident, E, fstat = FALSE, cov_type = "HC0") {
  N <- ncol(U)
  k <- E

  # Reduced-form covariance on the full VAR sample
  Sigma <- crossprod(U) / nrow(U)

  # Identification sample
  Uid <- U[ident, , drop = FALSE]
  Zid <- Z[ident, , drop = FALSE]
  Tn <- nrow(Uid)

  # First stage: policy residuals on constant and proxies
  fs <- stats::lm.fit(cbind(1, Zid), Uid[, seq_len(k), drop = FALSE])
  uhat <- matrix(fs$fitted.values, nrow = Tn, ncol = k)

  # Second stage: non-policy residuals on fitted policy residuals (2SLS)
  ss <- stats::lm.fit(cbind(1, uhat), Uid[, (k + 1):N, drop = FALSE])
  ss_coef <- matrix(ss$coefficients, nrow = k + 1, ncol = N - k)
  b21ib11 <- t(ss_coef[2:(k + 1), , drop = FALSE]) # (N - k) x k

  Sig11 <- Sigma[seq_len(k), seq_len(k), drop = FALSE]
  Sig21 <- Sigma[(k + 1):N, seq_len(k), drop = FALSE]
  Sig22 <- Sigma[(k + 1):N, (k + 1):N, drop = FALSE]

  ZZp <- b21ib11 %*% Sig11 %*% t(b21ib11) -
    (Sig21 %*% t(b21ib11) + b21ib11 %*% t(Sig21)) + Sig22
  b12b12p <- t(Sig21 - b21ib11 %*% Sig11) %*%
    solve(ZZp, Sig21 - b21ib11 %*% Sig11)
  b11b11p <- Sig11 - b12b12p
  L1 <- tryCatch(
    t(chol(b11b11p)),
    error = function(e) {
      stop("Identification failed: the implied covariance matrix of the ",
        "policy shocks is not positive definite. Check the instruments ",
        "and the ordering of the variables in y.", call. = FALSE)
    }
  )

  Psi <- rbind(L1, b21ib11 %*% L1)
  dimnames(Psi) <- NULL

  # First-stage diagnostics, equation by equation (skipped in the bootstrap)
  FStat <- FStatRobust <- R2 <- NULL
  if (fstat) {
    fsx <- lapply(seq_len(k), function(e) {
      .firststage_f(Uid[, e], Zid, cov_type = cov_type)
    })
    FStat <- vapply(fsx, `[[`, numeric(1), "F_hom")
    FStatRobust <- vapply(fsx, `[[`, numeric(1), "F_rob")
    R2 <- vapply(fsx, `[[`, numeric(1), "R2")
  }

  list(Psi = Psi, FStat = FStat, FStatRobust = FStatRobust, R2 = R2)
}

# IRFs by recursion through the VAR, normalized so that the impact response of
# policy variable e to shock e equals norm.
.svariv_irf <- function(Psi, Phi, H, E, norm, cum) {
  irf <- computeirf(Psi = Psi, Phi = Phi, H = H, cum = cum)
  for (e in seq_len(E)) {
    irf[, , e] <- irf[, , e] / irf[1, e, e] * norm
  }
  irf
}

# Generate artificial data from the estimated VAR with given residuals,
# starting from the first P observations of y.
.svariv_simulate <- function(y, X, fit, U, P) {
  Nobs <- nrow(y)
  N <- ncol(y)
  B <- matrix(fit$Phi, nrow = N, ncol = N * P) # [Phi_1, ..., Phi_P]
  yb <- y
  for (t in (P + 1):Nobs) {
    lags <- c(t(yb[(t - 1):(t - P), , drop = FALSE])) # lag 1 vars, lag 2 vars, ...
    det <- fit$C
    if (!is.null(X)) {
      det <- det + fit$D %*% X[t, ]
    }
    yb[t, ] <- det + B %*% lags + U[t - P, ]
  }
  yb
}
