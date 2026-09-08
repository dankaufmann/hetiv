#' Bias correction of impulse responses (Herbst & Johannsen, 2024)
#'
#' Ported from https://github.com/ckwolf92/lp_var_nberma/tree/main
#' Reference: Montiel Olea, José Luis (https://www.joseluismontielolea.com),
#'   Mikkel Plagborg-Møller (https://www.mikkelpm.com),
#'   Eric Qian (https://www.eric-qian.com) and Christian K. Wolf
#'   (https://www.christiankwolf.com/) (2025),
#'   "Local Projections or VARs? A Primer for Macroeconomists"
#'   (https://arxiv.org/abs/2503.17144).
#'
#' @param irs Numeric vector of length H (1 x H). Impulse response; the first
#'   element is the response on impact.
#' @param w   T x k matrix of control variables used in the local projections
#'
#' @return Numeric vector of length H: the bias-corrected impulse response.
#' @export
biascorr <- function(irs, w) {

  irs     <- as.numeric(irs)
  w       <- as.matrix(w)
  irs_hor <- length(irs) - 1L
  T       <- nrow(w)

  # ACF term in bias correction
  acf_corr <- rep(NA_real_, irs_hor)
  w        <- sweep(w, 2, colMeans(w))        # de-mean (column-wise)
  Sigma_0  <- stats::cov(w)                   # var-cov (k x k)
  for (j in seq_len(irs_hor)) {
    # j-th autocovariance: (w(1:end-j,:)' * w(j+1:end,:)) / (T-j-1)
    Sigma_j     <- crossprod(w[1:(T - j), , drop = FALSE],
                             w[(j + 1):T, , drop = FALSE]) / (T - j - 1)
    acf_corr[j] <- 1 + sum(diag(solve(Sigma_0, Sigma_j)))  # 1 + trace(Sigma_0 \ Sigma_j)
  }

  # Iterate on bias correction
  irs_corr <- irs
  for (h in seq_len(irs_hor)) {
    # acf_corr(1:h) * irs_corr(h:-1:1)' : dot product with reversed IRF
    irs_corr[h + 1] <- irs[h + 1] +
      (1 / (T - h)) * sum(acf_corr[1:h] * rev(irs_corr[1:h]))
  }

  irs_corr
}
