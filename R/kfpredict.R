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
#'   `SigR`. If `FALSE`, the raw Kalman filter projection is returned.
#'
#' @return A numeric matrix (T x E) of extracted structural shocks, with the
#'   same number of rows as `et` and one column per shock dimension.
#'
#' @references
#' Burri, M. and Kaufmann, D. (2024). Measuring monetary policy shocks.
#' IRENE Working Paper 24-03, IRENE Institute of Economic Research, University of Neuchâtel.
#'
#' @importFrom MASS ginv
#' @importFrom matrixcalc vec
#'
#' @export
kfpredict <- function(Sig, SigR, Psi, et, tol, scale = TRUE){

  # Clip tolerance from below to avoid numerical issues in ginv
  tol <- max(sqrt(.Machine$double.eps), tol)
  if(is.na(tol)){
    tol <- sqrt(.Machine$double.eps)
  }

  eps <- as.matrix(et[, 1:dim(Psi)[2]])
  eps[, ] <- NA

  if(scale == TRUE){
    # Back out implied shock variances from the difference in covariance matrices
    # between event and control days (Sig - SigR = Psi * SigEps * Psi')
    # Unit impact normalization on Psi means we need to recover the scale of
    # the shocks to renormalize them to unit variance
    q      <- matrixcalc::vec(Sig - SigR)
    A      <- sapply(1:ncol(Psi), function(i) c(Psi[, i] %*% t(Psi[, i])))
    sig    <- c(MASS::ginv(A) %*% q)

    if(length(sig) > 1){
      SigEps   <- diag(sig)
      myScale  <- solve(diag(sqrt(1/sig)))
      if(any(is.na(myScale))){
        myScale <- diag(abs(sig))
      }
    }else{
      SigEps  <- sig
      myScale <- 1/sqrt(1/sig)
      if(any(is.na(myScale))){
        myScale <- abs(sig)
      }
    }

  }else{
    myScale <- 1
  }

  # Kalman filter projection: eps_t = scale * Psi' * Sig^{-1} * e_t
  for(t in 1:dim(eps)[1]){
    eps[t, ] <- myScale %*% t(Psi) %*% MASS::ginv(Sig, tol = tol) %*% as.matrix(et[t, ])
  }

  return(eps)
}
