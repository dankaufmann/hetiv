# Internal helpers -------------------------------------------------------

# Symmetric matrix power A^p via eigendecomposition (A must be symmetric PSD)
.mat_pow <- function(A, p) {
  e <- eigen(A, symmetric = TRUE)
  e$vectors %*% diag(e$values^p, nrow = length(e$values)) %*% t(e$vectors)
}

# Commutation matrix K_{m,n}: transforms vec(A_{m x n}) to vec(A'_{n x m})
.spKgen <- function(m, n) {
  idx <- as.vector(t(matrix(seq_len(m * n), m, n)))
  diag(m * n)[idx, , drop = FALSE]
}

# Nearest symmetric positive definite matrix (Higham 1988)
.nearestSPD <- function(A) {
  B    <- (A + t(A)) / 2
  svdB <- svd(B)
  H    <- svdB$v %*% diag(svdB$d, nrow = length(svdB$d)) %*% t(svdB$v)
  Ahat <- (B + H) / 2
  Ahat <- (Ahat + t(Ahat)) / 2
  k_iter <- 0L
  repeat {
    R_try <- tryCatch(chol(Ahat), error = function(e) NULL)
    if (!is.null(R_try)) break
    k_iter <- k_iter + 1L
    mineig <- min(eigen(Ahat, symmetric = TRUE, only.values = TRUE)$values)
    Ahat   <- Ahat + (-mineig * k_iter^2 + .Machine$double.eps * abs(mineig)) * diag(nrow(A))
  }
  Ahat
}

# Thin QR decomposition with sign convention (diagonal of R >= 0)
.myQR <- function(XX, k) {
  qr_res <- qr(XX)
  Q <- qr.Q(qr_res)[, seq_len(k), drop = FALSE]
  R <- qr.R(qr_res)[seq_len(k), seq_len(k), drop = FALSE]
  d <- sign(diag(R))
  d[d == 0] <- 1
  Q %*% diag(d, k, k)
}

# Matrix (Frobenius) inner product
.iprod <- function(x, y) Re(sum(Conj(x) * y))

# Objective + gradient for Stiefel optimisation; x is a K x N matrix on the manifold
.objL0 <- function(x, M1, M2PsiM2, X1, N, K) {
  L0    <- t(x)
  vecL0 <- as.vector(L0)
  QLL   <- kronecker(kronecker(diag(N), L0), L0)

  Mobj <- M1 %*% QLL %*% M2PsiM2 %*% t(QLL) %*% t(M1) / K
  Mobj <- 0.5 * (Mobj + t(Mobj))
  Mobj <- .nearestSPD(Mobj)

  ev   <- eigen(Mobj, symmetric = TRUE)$vectors[, 1L, drop = TRUE]
  fval <- -as.numeric(t(ev) %*% Mobj %*% ev)

  left1   <- matrix(as.vector(t(ev) %*% M1 %*% QLL %*% M2PsiM2), 1L)
  left2   <- matrix(as.vector(t(ev) %*% M1), 1L)
  kron_L  <- kronecker(left1, left2)
  kron_R  <- kronecker(diag(N * K), matrix(vecL0, ncol = 1L))
  grad_v  <- 2 * (kron_L %*% X1 %*% kron_R)
  gradient <- t(-matrix(grad_v, N, K))

  list(fval = fval, gradient = gradient)
}

# Curvilinear search on the Stiefel manifold (Wen & Yin 2013, OptStiefelGBB)
.OptStiefelGBB <- function(X, fun, opts = list()) {
  n <- nrow(X); k <- ncol(X)

  xtol  <- if (!is.null(opts$xtol))   opts$xtol   else 1e-6
  gtol  <- if (!is.null(opts$gtol))   opts$gtol   else 1e-6
  ftol  <- if (!is.null(opts$ftol))   opts$ftol   else 1e-12
  tau   <- if (!is.null(opts$tau))    opts$tau    else 1e-3
  rhols <- if (!is.null(opts$rhols))  opts$rhols  else 1e-4
  eta   <- if (!is.null(opts$eta))    opts$eta    else 0.1
  gamma <- if (!is.null(opts$gamma))  opts$gamma  else 0.85
  nt    <- if (!is.null(opts$nt))     opts$nt     else 5L
  mxitr <- if (!is.null(opts$mxitr)) opts$mxitr  else 1000L
  tiny  <- if (!is.null(opts$tiny))   opts$tiny   else 1e-13

  res  <- fun(X)
  F    <- res$fval
  G    <- res$gradient
  GX   <- t(G) %*% X
  dtX  <- G - X %*% GX
  nrmG <- norm(dtX, "F")

  Q    <- 1; Cval <- F
  crit_mat <- matrix(1, mxitr, 3)
  itr      <- mxitr

  for (itr in seq_len(mxitr)) {
    XP <- X; FP <- F; dtXP <- dtX

    deriv <- rhols * nrmG^2
    nls   <- 1L
    repeat {
      X_try <- .myQR(XP - tau * dtX, k)
      if (norm(t(X_try) %*% X_try - diag(k), "F") > tiny)
        X_try <- .myQR(X_try, k)
      res <- fun(X_try); F <- res$fval; G <- res$gradient
      if (F <= Cval - tau * deriv || nls >= 5L) break
      tau <- eta * tau; nls <- nls + 1L
    }
    X <- X_try

    GX   <- t(G) %*% X
    dtX  <- G - X %*% GX
    nrmG <- norm(dtX, "F")
    S    <- X - XP
    XDiff <- norm(S, "F") / sqrt(n)
    FDiff <- abs(FP - F) / (abs(FP) + 1)

    Y  <- dtX - dtXP
    SY <- abs(.iprod(S, Y))
    if (itr %% 2L == 0L) {
      tau <- norm(S, "F")^2 / SY
    } else {
      tau <- SY / norm(Y, "F")^2
    }
    tau <- max(min(tau, 1e20), 1e-20)

    crit_mat[itr, ] <- c(nrmG, XDiff, FDiff)
    start_r <- max(1L, itr - nt + 1L)
    mcrit   <- colMeans(crit_mat[start_r:itr, , drop = FALSE])

    if ((XDiff < xtol && FDiff < ftol) || nrmG < gtol ||
        all(mcrit[2:3] < 10 * c(xtol, ftol))) break
    Qp <- Q; Q <- gamma * Qp + 1; Cval <- (gamma * Qp * Cval + F) / Q
  }

  if (norm(t(X) %*% X - diag(k), "F") > 1e-13) {
    X   <- .myQR(X, k)
    res <- fun(X); F <- res$fval
  }

  list(X = X, out = list(fval = F, itr = itr))
}


# Critical values (internal) ---------------------------------------------

# Critical values for gweaktest(); not exported — called internally.
# Reference: Lewis & Mertens (2025).
gweaktest_critical_values <- function(W, K,
                                      Sig     = NULL,
                                      alfa    = 0.05,
                                      tau     = 0.10,
                                      points  = 1000L,
                                      target  = "beta",
                                      crit    = "abs") {

  N <- as.integer(round(nrow(W) / K - 1))

  if (K < N) stop("Not identified: require K >= N instruments.")

  vecIK  <- as.vector(diag(K))
  vecIN  <- as.vector(diag(N))

  RNK  <- kronecker(diag(N),     matrix(vecIK, K^2, 1))
  RNN  <- kronecker(diag(N),     matrix(vecIN, N^2, 1))
  RNpK <- kronecker(diag(N + 1), matrix(vecIK, K^2, 1))

  M1 <- t(RNN) %*% (diag(N^3) + kronecker(.spKgen(N, N), diag(N)))
  M2 <- RNK %*% t(RNK) / (1 + N) - diag(N * K^2)

  W2  <- W[(K + 1):nrow(W), (K + 1):ncol(W), drop = FALSE]
  W12 <- W[1:K,              (K + 1):ncol(W), drop = FALSE]

  Phi          <- t(RNK) %*% kronecker(W2, diag(K)) %*% RNK
  Phi_sqrt_inv <- .mat_pow(Phi / K, -0.5)
  KPhiInv      <- kronecker(Phi_sqrt_inv, diag(K))
  Sigma        <- KPhiInv %*% W2 %*% t(KPhiInv)

  inner  <- KPhiInv %*% t(rbind(W12, W2))
  Psibar <- kronecker(inner, diag(K)) %*% RNpK

  if (crit == "rel") {
    RWR <- t(RNpK) %*% kronecker(W, diag(K)) %*% RNpK
    Psi <- Psibar %*% .mat_pow(RWR, -0.5)
  } else if (crit == "abs") {
    if (is.null(Sig))
      stop("Error covariance matrix Sig required for crit = 'abs'.")
    Sig22 <- Sig[2:(N + 1), 2:(N + 1), drop = FALSE]
    scl   <- norm(.mat_pow(Phi, -0.5) %*% .mat_pow(Sig22, 0.5), "2")
    Psi   <- Psibar %*% .mat_pow(Sig, -0.5) * scl
  } else {
    stop("crit must be 'abs' or 'rel'.")
  }

  M2PsiM2 <- M2 %*% tcrossprod(Psi) %*% t(M2)

  X1_d <- diag(N^2 * K^2) + .spKgen(N * K, N * K)
  X1_c <- kronecker(kronecker(diag(K), .spKgen(K, N)), diag(N))  %*% X1_d
  X1_b <- kronecker(matrix(vecIN, N^2, 1), diag(K^2 * N^2))     %*% X1_c
  X1   <- kronecker(kronecker(diag(N), .spKgen(K^2, N)), diag(N^2)) %*% X1_b

  Bmax       <- numeric(3)
  M2Psi_norm <- norm(M2 %*% Psi, "2")
  Psi_norm   <- norm(Psi, "2")

  Bmax[2] <- if (K > N + 1L) {
    val <- min(sqrt(2 * (N + 1) / K) * M2Psi_norm, Psi_norm)
    if (N == 1L && crit == "rel") min(val, 1) else val
  } else {
    val <- max(sqrt(2 * (N + 1) / K) * M2Psi_norm, Psi_norm)
    if (N == 1L && crit == "rel") min(val, 1) else val
  }

  if (K > N + 1L) {
    opts_s   <- list(mxitr = 1000L, xtol = 1e-5, gtol = 1e-5, ftol = 1e-7)
    obj_fun  <- function(x) .objL0(x, M1, M2PsiM2, X1, N, K)
    Bmax_vec <- numeric(points)
    for (iter in seq_len(points)) {
      L0_init        <- qr.Q(qr(matrix(rnorm(K * K), K, K)))[, seq_len(N), drop = FALSE]
      res_s          <- .OptStiefelGBB(L0_init, obj_fun, opts_s)
      Bmax_vec[iter] <- sqrt(max(-res_s$out$fval, 0))
    }
    Bmax[1] <- max(Bmax_vec)
    Bmax[3] <- (K - (1 + N)) / K
  } else {
    Bmax[1] <- Bmax[2]
    Bmax[3] <- NaN
  }

  if (K == N) {
   #message("Model is just-identified: test is for median bias.")
    if (N == 1L) tau <- tau / 0.455
  }

  if (is.numeric(target) && crit == "abs") {
    j_idx    <- as.integer(target)
    iPhi     <- .mat_pow(Phi, -0.5)
    row_norm <- sqrt(sum(iPhi[j_idx, ]^2))
    Sig22    <- Sig[2:(N + 1), 2:(N + 1), drop = FALSE]
    full_nrm <- norm(iPhi %*% .mat_pow(Sig22, 0.5), "2")
    tau      <- tau / (sqrt(Sig[j_idx + 1, j_idx + 1]) * row_norm) * full_nrm
  }

  Sig_pows <- list(Sigma, Sigma %*% Sigma, Sigma %*% Sigma %*% Sigma)
  Sig_norm <- norm(Sigma, "2")

  cv <- numeric(3)

  for (j in 1:3) {
    lmin <- Bmax[j] / tau

    if (j < 3) {
      cums <- numeric(3)
      for (nn in 1:3) {
        inner_n  <- t(RNK) %*% kronecker(Sig_pows[[nn]], diag(K)) %*% RNK
        cums[nn] <- 2^(nn - 1) * factorial(nn - 1) *
          (norm(inner_n, "2") + nn * K * lmin * Sig_norm^(nn - 1))
      }

      ome <- cums[2] / cums[3]
      nu  <- 8 * cums[2] * ome^2
      cc  <- qchisq(1 - alfa, nu)

      fun_phiz <- function(z) {
        arg <- 1 + (z - cums[1]) / (2 * cums[2] * ome)
        ome * arg^(nu / 2 - 1) * exp(-nu / 2 * arg) *
          nu^(nu / 2 - 1) / (2^(nu / 2 - 2)) / gamma(nu / 2)
      }
      G1fun <- function(q) {
        -0.5 * (q - 2 * nu * (nu - 2) / q + nu) +
          1.5 * nu * (log(q / 2) - digamma(nu / 2))
      }
      G2fun <- function(q) {
        0.5 * (q - nu * (nu - 2) / q) -
          nu * (log(q / 2) - digamma(nu / 2))
      }
      D1fun <- function(q) {
        u <- (q - cums[1]) / (2 * cums[2] * ome)
        (1 + (q - cums[1]) * 2 * ome) / (2 * cums[2] * ome) *
          (1 + u)^(-1) * fun_phiz(q)
      }
      D2fun <- function(q) {
        fun_phiz(q) / cums[2] * G1fun(nu + (q - cums[1]) * 4 * ome)
      }
      D3fun <- function(q) {
        G2fun(nu + (q - cums[1]) * 4 * ome) / cums[3] * fun_phiz(q)
      }

      lb       <- (cc - nu) / (4 * ome) + cums[1]
      kt_cond1 <- suppressWarnings(tryCatch(integrate(D1fun, lb, Inf)$value, error = function(e) 0))
      kt_cond2 <- suppressWarnings(tryCatch(integrate(D2fun, lb, Inf)$value, error = function(e) 0))
      kt_cond3 <- suppressWarnings(tryCatch(integrate(D3fun, lb, Inf)$value, error = function(e) 0))
      kt_cond  <- (kt_cond1 >= 0) && (kt_cond2 >= 0) && (kt_cond3 >= 0)

      if (!kt_cond) {
        cums_old <- cums
        if (N > 1L) {
          fn_opt <- function(x) {
            nu_x <- 8 * x[2] * (x[2] / x[3])^2
            -((qchisq(1 - alfa, nu_x) - nu_x) / 4 / (x[2] / x[3]) + x[1])
          }
          res_opt <- suppressWarnings(
            optim(cums, fn_opt, method = "L-BFGS-B",
                  lower = rep(0.01, 3), upper = cums_old)
          )
          cums <- res_opt$par
          ome  <- cums[2] / cums[3]
          nu   <- 8 * cums[2] * ome^2
        } else {
          fn_opt1 <- function(x) {
            nu_x <- 8 * x[1] * (x[1] / x[2])^2
            -((qchisq(1 - alfa, nu_x) - nu_x) / 4 / (x[1] / x[2]) + cums[1])
          }
          res_opt <- suppressWarnings(
            optim(cums[2:3], fn_opt1, method = "L-BFGS-B",
                  lower = rep(0.01, 2), upper = cums_old[2:3])
          )
          knew      <- res_opt$par
          ome       <- knew[1] / knew[2]
          nu        <- 8 * knew[1] * ome^2
          cums[2:3] <- knew
        }
        cc <- qchisq(1 - alfa, nu)
      }

      cv[j] <- ((cc - nu) / 4 / ome + cums[1]) / K

    } else {
      cv[j] <- if (!is.nan(lmin) && is.finite(lmin)) {
        qchisq(1 - alfa, K, ncp = K * lmin) / K
      } else {
        NaN
      }
    }
  }

  list(
    wiv_cv             = cv[1],
    wiv_cv_simplified  = cv[2],
    wiv_cv_sy          = cv[3]
  )
}


# Exported function ------------------------------------------------------

#' Generalised weak instrument test
#'
#' Tests for weak instruments with multiple endogenous regressors using the
#' generalised minimum eigenvalue statistic of Lewis and Mertens (2025). The
#' test is robust to heteroscedasticity and serial correlation and nests the
#' classical Stock-Yogo (2005) test as a special case.
#'
#' @param y Regressand (T x 1 numeric vector or matrix).
#' @param Y Endogenous regressors (T x N numeric matrix).
#' @param X Exogenous regressors (T x Nx numeric matrix). A constant column is
#'   added automatically if one is absent or if `X` has zero columns.
#' @param Z Instruments (T x K numeric matrix). Requires K >= N.
#' @param cov_type HAR covariance estimator: `"EHW"` (Eicker-Huber-White,
#'   default) or `"NW"` (Newey-West with Lazarus et al. (2018) bandwidth).
#' @param alfa Nominal significance level (default `0.05`).
#' @param tau Bias tolerance: maximum acceptable relative (or absolute, see
#'   `crit`) bias of the 2SLS estimator (default `0.10`).
#' @param points Number of random starting points for the Stiefel manifold
#'   optimisation used to compute the sharp critical value when K > N + 1
#'   (default `1000`).
#' @param target Either `"beta"` (default) to test the full coefficient vector,
#'   or a positive integer `j <= N` to target the single coefficient
#'   `beta_j`.
#' @param crit Bias criterion: `"abs"` (absolute bias, default) or `"rel"`
#'   (relative bias). `"abs"` requires the error covariance matrix; `"rel"`
#'   does not.
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{`nobs`}{Number of complete observations used.}
#'     \item{`beta_2SLS`}{2SLS point estimate(s).}
#'     \item{`target`}{Description of the targeted parameter.}
#'     \item{`criterion`}{Bias criterion used (`"abs"` or `"rel"`).}
#'     \item{`gmin_generalized`}{Generalised minimum eigenvalue test statistic
#'       (Lewis-Mertens).}
#'     \item{`gmin_generalized_critical_value`}{Sharp critical value via
#'       Stiefel optimisation.}
#'     \item{`gmin_generalized_critical_value_simplified`}{Conservative
#'       simplified critical value (closed-form bound).}
#'     \item{`stock_yogo_test_statistic`}{Stock-Yogo test statistic under the
#'       Nagar approximation.}
#'     \item{`stock_yogo_critical_value_nagar`}{Stock-Yogo critical value under
#'       the Nagar approximation.}
#'   }
#'
#' @references
#' Lazarus, E., Lewis, D. J., Stock, J. H. and Watson, M. W. (2018).
#' HAR Inference: Recommendations for Practice. *Journal of Business &
#' Economic Statistics*, 36(4), 541–559.
#'
#' Lewis, D. J. and Mertens, K. (2025). A Robust Test for Weak Instruments for
#' 2SLS with Multiple Endogenous Regressors. *The Review of Economic Studies*,
#' DOI: 10.1093/restud/rdaf103.
#'
#' Stock, J. H. and Yogo, M. (2005). Testing for Weak Instruments in Linear IV
#' Regression. In D. W. K. Andrews and J. H. Stock (Eds.), *Identification and
#' Inference for Econometric Models: Essays in Honor of Thomas Rothenberg*,
#' pp. 80–108. Cambridge University Press.
#'
#' @export
gweaktest <- function(y, Y, X, Z,
                      cov_type = "EHW",
                      alfa     = 0.05,
                      tau      = 0.10,
                      points   = 1000L,
                      target   = "beta",
                      crit     = "abs") {

  # Force column orientation
  y <- as.matrix(y);  if (ncol(y) > nrow(y)) y <- t(y)
  Y <- as.matrix(Y);  if (ncol(Y) > nrow(Y)) Y <- t(Y)
  Z <- as.matrix(Z);  if (ncol(Z) > nrow(Z)) Z <- t(Z)
  X <- as.matrix(X);  if (ncol(X) > nrow(X)) X <- t(X)

  # Drop missing observations
  obs_data <- cbind(y, Y, Z, X)
  sel      <- complete.cases(obs_data)
  y <- y[sel, , drop = FALSE]
  Y <- Y[sel, , drop = FALSE]
  Z <- Z[sel, , drop = FALSE]
  X <- X[sel, , drop = FALSE]

  Tobs <- nrow(y)

  # Add constant to X if absent; remove zero-variance columns
  if (ncol(X) > 0) {
    col_var <- apply(X, 2, var)
    X       <- cbind(X[, col_var != 0, drop = FALSE], 1)
  } else {
    X <- matrix(1, Tobs, 1)
  }

  Nx <- ncol(X)
  N  <- ncol(Y)
  K  <- ncol(Z)

  # Partial out exogenous regressors
  Zo  <- Z - X %*% qr.solve(X, Z)
  Zo  <- sweep(Zo, 2, colMeans(Zo), "-")
  Zo  <- Zo %*% .mat_pow(crossprod(Zo) / T, -0.5)   # normalise: Zo'Zo/T = I_K

  Yo  <- Y - X %*% qr.solve(X, Y)
  yo  <- y - X %*% qr.solve(X, y)

  PYo <- Zo %*% qr.solve(Zo, Yo)
  Pyo <- Zo %*% qr.solve(Zo, yo)

  betahat <- qr.solve(PYo, Pyo)   # 2SLS estimate

  # Residuals
  v1 <- yo  - Pyo
  v2 <- Yo  - PYo
  e  <- cbind(v1, v2)

  # Score matrix ZV: T x K(N+1)
  ZV <- kronecker(matrix(1, 1, N + 1), Zo) *
        e[, rep(seq_len(N + 1), each = K), drop = FALSE]

  # HAR covariance matrices
  if (cov_type == "NW") {
    L_bw    <- ceiling(1.3 * Tobs^0.5)
    W_mat   <- crossprod(ZV) / Tobs
    Sig_mat <- crossprod(e)  / Tobs
    for (jj in seq_len(L_bw)) {
      w_l  <- 1 - jj / L_bw
      acv  <- crossprod(ZV[(jj + 1):Tobs, ], ZV[1:(Tobs - jj), ]) / Tobs +
               crossprod(ZV[1:(Tobs - jj), ], ZV[(jj + 1):Tobs, ]) / Tobs
      W_mat   <- W_mat   + w_l * acv
      acve <- crossprod(e[(jj + 1):Tobs, ], e[1:(Tobs - jj), ]) / Tobs +
               crossprod(e[1:(Tobs - jj), ], e[(jj + 1):Tobs, ]) / Tobs
      Sig_mat <- Sig_mat + w_l * acve
    }
  } else {
    W_mat   <- crossprod(ZV) / Tobs
    Sig_mat <- crossprod(e)  / Tobs
  }

  df_corr <- Tobs / (Tobs - K - Nx)
  W_mat   <- W_mat   * df_corr
  Sig_mat <- Sig_mat * df_corr

  # Generalised minimum eigenvalue statistic
  vecIK  <- as.vector(diag(K))
  RNK    <- kronecker(diag(N), matrix(vecIK, K^2, 1))
  W2     <- W_mat[(K + 1):nrow(W_mat), (K + 1):ncol(W_mat), drop = FALSE]
  Phi    <- t(RNK) %*% kronecker(W2, diag(K)) %*% RNK

  Phi_half_inv <- .mat_pow(Phi, -0.5)
  A_stat       <- Phi_half_inv %*% crossprod(PYo) %*% Phi_half_inv
  gmin_gen     <- min(eigen(A_stat, symmetric = TRUE)$values)

  # Critical values
  Sig_arg <- if (crit == "rel") NULL else Sig_mat

  cv_list <- gweaktest_critical_values(
    W      = W_mat,
    K      = K,
    Sig    = Sig_arg,
    alfa   = alfa,
    tau    = tau,
    points = points,
    target = target,
    crit   = crit
  )

  # Stock-Yogo statistic (Nagar approximation)
  Svv      <- crossprod(v2) / (Tobs - K - Nx)
  ZtZ_inv  <- solve(crossprod(Zo))
  inner_sy <- .mat_pow(Svv, -0.5) %*%
               (t(Yo) %*% Zo %*% ZtZ_inv %*% t(Zo) %*% Yo) %*%
               .mat_pow(Svv, -0.5) / K
  gmin_sy  <- min(eigen(inner_sy, symmetric = TRUE)$values)

  list(
    nobs                                        = Tobs,
    beta_2SLS                                   = as.vector(betahat),
    target                                      = if (identical(target, "beta")) "beta vector"
                                                  else paste0("beta_", target),
    criterion                                   = crit,
    gmin_generalized                            = gmin_gen,
    gmin_generalized_critical_value             = cv_list$wiv_cv,
    gmin_generalized_critical_value_simplified  = cv_list$wiv_cv_simplified,
    stock_yogo_test_statistic                   = gmin_sy,
    stock_yogo_critical_value_nagar             = cv_list$wiv_cv_sy
  )
}
