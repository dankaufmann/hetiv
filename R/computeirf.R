#' Compute impulse response functions for a VAR(P)
#'
#' Recursively computes IRFs given an impact matrix and VAR coefficient array.
#' Optionally cumulates responses for selected variables.
#'
#' @param Psi Impact matrix, dimension `N x R`.
#' @param Phi Array of VAR coefficient matrices, dimension `N x N x P`.
#' @param H Horizon (number of periods) for which to compute IRFs.
#' @param cum Logical scalar or logical vector of length `N`. If `TRUE` for
#'   variable `i`, the IRF is cumulated via `cumsum`. A single value is applied to all variables.
#'
#' @return An array:
#' \describe{
#'   \item{irf}{Array of impulse responses, dimension `H x N x R`. Row names
#'     are labelled `0` to `H-1`.}
#' }
#'
#' @export
computeirf <- function(Psi, Phi, H, cum) {
  if (dim(Phi)[1] != dim(Phi)[2])
    stop("Each Phi slice must be a square matrix (N x N).")
  if (dim(Psi)[1] != dim(Phi)[1])
    stop("Psi and Phi must have the same first dimension (N).")
  if (H < 1)
    stop("H must be a positive integer.")

  N <- dim(Phi)[1]
  P <- dim(Phi)[3]
  R <- dim(Psi)[2]

  if (!is.logical(cum) || anyNA(cum) || !(length(cum) %in% c(1, N)))
    stop("cum must be a non-missing logical scalar or a logical vector of length nrow(Psi).")
  if (length(cum) == 1) {
    cum <- rep(cum, N)
  }

  irf <- array(0, dim = c(H, N, R))

  for (h in seq_len(H)) {
    if (h == 1) {
      irf[h, , ] <- Psi
    } else {
      for (p in seq_len(min(P, h - 1))) {
        irf[h, , ] <- irf[h, , ] + Phi[, , p] %*% irf[h - p, , ]
      }
    }
  }
  dimnames(irf)[[1]] <- 0:(H - 1)

  # Cumulate responses for selected variables
  for (i in seq_len(N)) {
    if (cum[i] == TRUE) {
      for (r in seq_len(R)) {
        irf[, i, r] <- cumsum(irf[, i, r])
      }
    }
  }

  return(irf)
}
