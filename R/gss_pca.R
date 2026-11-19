#' Monetary policy surprise factors from a panel of financial market series
#'
#' Estimates one or two factors from a panel of event-window changes of
#' financial market series (rows are events, columns are series), following
#' Gurkaynak, Sack and Swanson (2005):
#' \itemize{
#'   \item `n_factors = 1`: a single common factor ("overall" surprise), the
#'     first principal component of the standardised panel.
#'   \item `n_factors = 2`: two factors rotated such that the second factor
#'     ("path") has a zero loading on the first column of `X`. The first
#'     factor ("target") is then the only factor that moves the first series.
#' }
#'
#' @details Missing values are handled by EM principal components (Stock and
#'   Watson, 2002; Swanson, 2021): missing cells are initialised at the mean
#'   and iteratively replaced by their fit from the `n_factors` principal
#'   components until convergence.
#'
#'   The factors are normalised against series in `X` through regression
#'   coefficients, as in Gurkaynak, Sack and Swanson (2005):
#'   \itemize{
#'     \item `n_factors = 1`: the factor is rescaled such that the coefficient
#'       of a regression of `X[, norm_anchor[1]]` on the factor is one (a
#'       one-unit move of the factor is a one-unit move of its anchor series).
#'     \item `n_factors = 2`: the target factor is rescaled in the same way
#'       against `X[, norm_anchor[1]]`. The path factor is rescaled such that
#'       its coefficient on `X[, norm_anchor[2]]` equals the coefficient of
#'       the target factor on the same series.
#'   }
#'   The anchor regressions use the rows where the anchor series is observed.
#'
#' @param X Numeric matrix or data frame (T x N) of event-window changes, with
#'   only the financial market series as columns. For `n_factors = 2`, the
#'   first column defines the target restriction (the path factor has a zero
#'   loading on it).
#' @param n_factors Integer. `1` (overall factor) or `2` (target and path
#'   factors, default).
#' @param norm_anchor Integer vector of length `n_factors` with the column
#'   positions in `X` against which the factors are normalised. Defaults to
#'   `1` or `c(1, 2)`. For `n_factors = 2`, `norm_anchor[2]` cannot be `1`.
#' @param tol Numeric. Convergence tolerance of the EM iterations.
#' @param maxit Integer. Maximum number of EM iterations.
#'
#' @return A list of class `gss_pca` with elements:
#'   \describe{
#'     \item{`factors`}{Matrix T x `n_factors` with columns `"overall"` or
#'       `"target"` and `"path"`, in the units of the anchor series.}
#'     \item{`loadings`}{Matrix N x `n_factors` of loadings of the
#'       standardised series on the (rotated and normalised) factors.}
#'     \item{`var_explained`}{Share of the variance of the standardised panel
#'       explained by each unrotated principal component.}
#'     \item{`restriction`}{For `n_factors = 2`: loading of the path factor on
#'       the first series (zero up to numerical error); `NULL` otherwise.}
#'     \item{`norm_anchor`}{Named character vector with the anchor series of
#'       each factor.}
#'     \item{`norm_coef`}{Named vector with the regression coefficient to
#'       which each factor was calibrated.}
#'     \item{`mu`, `sd`}{Means and standard deviations of the series used for
#'       standardisation.}
#'     \item{`X_filled`}{Standardised panel with missing cells replaced by
#'       their EM fit.}
#'     \item{`n_missing`, `iters`, `converged`, `call`}{Number of missing
#'       cells, number of EM iterations, convergence flag, and the call.}
#'   }
#'
#' @references
#' Gurkaynak, R. S., Sack, B. and Swanson, E. T. (2005). Do actions speak
#' louder than words? The response of asset prices to monetary policy actions
#' and statements. *International Journal of Central Banking*, 1(1), 55-93.
#'
#' Stock, J. H. and Watson, M. W. (2002). Macroeconomic forecasting using
#' diffusion indexes. *Journal of Business & Economic Statistics*, 20(2),
#' 147-162.
#'
#' Swanson, E. T. (2021). Measuring the effects of Federal Reserve forward
#' guidance and asset purchases on financial markets. *Journal of Monetary
#' Economics*, 118, 32-53.
#'
#' @importFrom stats coef setNames
#'
#' @examples
#' # Simulated panel: 5 series driven by a target and a path factor; the first
#' # series responds to the target factor only
#' set.seed(1)
#' f <- matrix(rnorm(400), ncol = 2)
#' L <- cbind(c(1, 0.8, 0.6, 0.4, 0.2), c(0, 0.3, 0.6, 0.8, 1))
#' X <- f %*% t(L) + matrix(rnorm(1000, sd = 0.1), ncol = 5)
#' colnames(X) <- paste0("r", 1:5)
#' X[3, 4] <- NA
#' res <- gss_pca(X, n_factors = 2)
#' res
#' head(res$factors)
#'
#' @export
gss_pca <- function(X, n_factors = 2, norm_anchor = seq_len(n_factors),
                    tol = 1e-8, maxit = 500) {

  cl <- match.call()

  # --- input checks ----------------------------------------------------------
  if (!n_factors %in% 1:2)
    stop("n_factors must be 1 (overall) or 2 (target/path).")
  X <- as.matrix(X)
  norm_anchor <- as.integer(norm_anchor)
  if (length(norm_anchor) != n_factors ||
      any(norm_anchor < 1) || any(norm_anchor > ncol(X)))
    stop("norm_anchor must give one column position of X per factor ",
         "(length ", n_factors, ", values in 1..", ncol(X), ").")
  if (n_factors == 2 && norm_anchor[2] == 1)
    stop("norm_anchor[2] cannot be 1: the path factor has zero loading on ",
         "column 1 by construction, so a coefficient on it is undefined.")
  if (!is.numeric(X))
    stop("X must be numeric: pass only the financial market series ",
         "(drop Date/meta columns before calling gss_pca).")
  if (ncol(X) <= n_factors)
    stop("Need more series (columns) than factors.")
  if (any(colSums(!is.na(X)) < 2))
    stop("Column(s) with fewer than 2 observations: ",
         paste(colnames(X)[colSums(!is.na(X)) < 2], collapse = ", "))
  empty_rows <- rowSums(!is.na(X)) == 0
  if (any(empty_rows))
    warning(sum(empty_rows), " row(s) have no observed series at all; ",
            "their factor values are pure EM interpolation (~0).")
  k <- n_factors
  nT <- nrow(X)
  N <- ncol(X)
  sn <- if (is.null(colnames(X))) paste0("x", seq_len(N)) else colnames(X)

  # --- EM principal components  --
  mu <- apply(X, 2, mean, na.rm = TRUE)
  sdev <- apply(X, 2, stats::sd, na.rm = TRUE)
  if (any(sdev == 0)) stop("Constant column(s): ", paste(sn[sdev == 0], collapse = ", "))
  Z    <- sweep(sweep(X, 2, mu, "-"), 2, sdev, "/")
  miss <- is.na(Z)
  Zf   <- Z; Zf[miss] <- 0                    # initialise missing at the mean (0)

  prev <- Inf; converged <- FALSE; it <- 0
  for (it in seq_len(maxit)) {
    C  <- crossprod(Zf) / (nT - 1)             # correlation-scale covariance
    eg <- eigen(C, symmetric = TRUE)
    V  <- eg$vectors[, 1:k, drop = FALSE]
    Zhat <- (Zf %*% V) %*% t(V)               # rank-k reconstruction
    Zf[miss] <- Zhat[miss]                    # E-step: refill missing cells
    obj <- sqrt(mean((Zf - Zhat)^2))
    if (abs(prev - obj) < tol) { converged <- TRUE; break }
    prev <- obj
  }
  if (!converged && any(miss))
    warning("EM did not converge within ", maxit, " iterations.")

  # final decomposition on the filled panel; unit-variance scores
  C   <- crossprod(Zf) / (nT - 1)
  eg  <- eigen(C, symmetric = TRUE)
  lam <- eg$values[1:k]
  V   <- eg$vectors[, 1:k, drop = FALSE]
  Fs  <- (Zf %*% V) %*% diag(1 / sqrt(lam), k, k)   # T x k, unit variance
  L   <- V %*% diag(sqrt(lam), k, k)                # N x k, Stata convention in original GSS code

  # --- 1 factor: overall surprise ---------------------------------------------
  if (k == 1) {
    b <- unname(stats::coef(stats::lm(X[, norm_anchor[1]] ~ Fs[, 1]))[2])
    factors  <- Fs * b
    loadings <- L / b
    colnames(factors) <- "overall"
    norm_coef <- stats::setNames(b, "overall")
  }

  # --- 2 factors: GSS rotation, path has zero loading on series 1 ------------
  if (k == 2) {
    # rotation weights, mirroring GSS Stata file (a's from series-1 loadings,
    # b's orthogonal with zero loading on series 1; valid since var(f)=1)
    a1 <- L[1, 1] / (L[1, 1] + L[1, 2])
    a2 <- L[1, 2] / (L[1, 1] + L[1, 2])
    b1 <- -a2 / (a1 - a2)
    b2 <-  a1 / (a1 - a2)

    target <- a1 * Fs[, 1] + a2 * Fs[, 2]
    path   <- b1 * Fs[, 1] + b2 * Fs[, 2]

    # rotated loadings: Z ~ target * l_t' + path * l_p'
    R    <- rbind(c(a1, a2), c(b1, b2))       # (target, path)' = R (f1, f2)'
    Lrot <- L %*% solve(R)                    # N x 2

    # --- GSS-style normalisation  --------------------------
    # target: rescale so coef(anchor1 ~ target) == 1 (into anchor1's units)
    b_t <- unname(stats::coef(stats::lm(X[, norm_anchor[1]] ~ target))[2])
    target <- target * b_t;  Lrot[, 1] <- Lrot[, 1] / b_t

    # path: rescale so its effect on anchor2 matches target's effect on it
    b_target_on_a2 <- unname(stats::coef(stats::lm(X[, norm_anchor[2]] ~ target))[2])
    b_path_on_a2   <- unname(stats::coef(stats::lm(X[, norm_anchor[2]] ~ path))[2])
    c_path <- b_path_on_a2 / b_target_on_a2
    path <- path * c_path;  Lrot[, 2] <- Lrot[, 2] / c_path

    factors  <- cbind(target = target, path = path)
    loadings <- Lrot
    norm_coef <- stats::setNames(c(1, b_target_on_a2), c("target", "path"))
  }

  rownames(factors)  <- rownames(X)
  dimnames(loadings) <- list(sn, colnames(factors))

  structure(list(
    factors       = factors,
    loadings      = loadings,
    var_explained = eg$values / sum(eg$values),
    restriction   = if (k == 2) unname(loadings[1, 2]) else NULL,
    norm_anchor   = stats::setNames(sn[norm_anchor], colnames(factors)),
    norm_coef     = norm_coef,
    mu = mu, sd = sdev,
    X_filled      = Zf,
    n_missing     = sum(miss),
    iters         = it,
    converged     = converged || !any(miss),
    call          = cl
  ), class = "gss_pca")
}

#' @rdname gss_pca
#' @param x An object of class `gss_pca`.
#' @param ... Ignored.
#' @export
print.gss_pca <- function(x, ...) {
  k <- ncol(x$factors)
  cat(sprintf("gss_pca: %d factor(s) [%s], %d events x %d series, %d missing cells\n",
              k, paste(colnames(x$factors), collapse = ", "),
              nrow(x$factors), nrow(x$loadings), x$n_missing))
  cat(sprintf("var explained (PC1..PC%d): %s\n", k,
              paste(sprintf("%.1f%%", 100 * x$var_explained[1:k]), collapse = ", ")))
  if (k == 2)
    cat(sprintf("path loading on '%s' (target restriction, should be 0): %.2e\n",
                rownames(x$loadings)[1], x$restriction))
  cat("normalisation: ", paste(sprintf("%s -> coef %.3f on '%s'",
      names(x$norm_anchor), x$norm_coef, x$norm_anchor), collapse = "; "), "\n", sep = "")
  cat("loadings:\n"); print(round(x$loadings, 3))
  invisible(x)
}
