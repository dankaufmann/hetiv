#' Extract monetary policy shocks via Kalman filter
#'
#' Extracts structural shocks from reduced-form residuals using the Kalman
#' filter formula. Given the impact matrix `Psi` and the covariance matrices
#' of residuals on event and control days, the function recovers the latent
#' shocks for each observation. Optionally rescales shocks to unit variance.
#'
#' @param Sig Numeric matrix (N x N). Covariance matrix of reduced-form
#'   residuals on policy event days (e.g. `et` from [hetiv()]).
#' @param SigR Numeric matrix (N x N). Covariance matrix of reduced-form
#'   residuals on control (non-event) days. Used to back out the implied
#'   shock variances when `scale = TRUE`.
#' @param Psi Numeric matrix (N x E). Impact matrix, i.e. the contemporaneous
#'   responses of all N variables to the E structural shocks (e.g. `Psi` from
#'   [hetiv()]).
#' @param et Numeric matrix or data frame (T x N). Reduced-form residuals on
#'   event days (e.g. `et` from [hetiv()]). Rows correspond to time
#'   periods, columns to variables.
#' @param tol Numeric scalar. Tolerance for the generalized inverse
#'   ([MASS::ginv()]). Clipped from below at `sqrt(.Machine$double.eps)`.
#' @param scale Logical. If `TRUE` (default), shocks are rescaled to unit
#'   variance using the implied shock variances recovered from `Sig` and
#'   `SigR`. This scaling assumes that structural shock variances are diagonal,
#'   so that `Sig - SigR = sum_i sigma_i Psi_i Psi_i'`. If `FALSE`, the raw
#'   Kalman filter projection is returned.
#'
#' @return A numeric matrix (T x E) of extracted structural shocks, with the
#'   same number of rows as `et` and one column per shock dimension.
#'
#' @references
#' Burri, M. and Kaufmann, D. (2026). Measuring monetary policy shocks.
#' IRENE Working Paper 24-03, IRENE Institute of Economic Research, University of Neuchâtel.
#'
#' @importFrom MASS ginv
#'
#' @examples
#' Sig <- diag(2)
#' SigR <- diag(c(0.5, 0.5))
#' Psi <- matrix(c(1, 0.5), nrow = 2)
#' et <- matrix(rnorm(20), ncol = 2)
#' kfpredict(Sig = Sig, SigR = SigR, Psi = Psi, et = et)
#'
#' @export
kfpredict <- function(Sig, SigR, Psi, et, tol = sqrt(.Machine$double.eps), scale = TRUE){

  Sig <- .as_numeric_matrix(Sig, "Sig", square = TRUE, allow_na = FALSE)
  scale <- .check_logical_scalar(scale, "scale")
  tol <- .check_numeric_scalar(tol, "tol", min = 0)
  tol <- max(tol, sqrt(.Machine$double.eps))
  Psi <- .as_numeric_matrix(Psi, "Psi", nrow = nrow(Sig), allow_na = FALSE)
  et <- .as_numeric_matrix(et, "et", ncol = nrow(Sig))

  E <- ncol(Psi)

  if(scale == TRUE){
    SigR <- .as_numeric_matrix(SigR, "SigR", square = TRUE, allow_na = FALSE)
    if (!identical(dim(Sig), dim(SigR)))
      stop("Sig and SigR must have the same dimensions.", call. = FALSE)
    # Back out implied shock variances from the difference in covariance matrices
    # between event and control days (Sig - SigR = Psi * SigEps * Psi')
    # Unit impact normalization on Psi means we need to recover the scale of
    # the shocks to renormalize them to unit variance
    q      <- as.vector(Sig - SigR)
    A      <- sapply(seq_len(ncol(Psi)), function(i) c(Psi[, i] %*% t(Psi[, i])))
    sig    <- c(MASS::ginv(A) %*% q)

    if(any(sig < 0))
      warning("Some estimated shock variances are negative (weak heteroskedasticity); ",
              "using absolute values for scaling.")

    myScale <- diag(sqrt(abs(sig)), nrow = E, ncol = E)

  }else{
    myScale <- diag(E)
  }

  # Pre-compute projection matrix: scale * Psi' * Sig^{-1}
  proj <- myScale %*% t(Psi) %*% MASS::ginv(Sig, tol = tol)

  # Kalman filter projection: eps_t = proj * e_t  (vectorised over T)
  eps <- t(proj %*% t(et))

  return(eps)
}
