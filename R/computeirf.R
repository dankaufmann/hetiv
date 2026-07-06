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
#' @examples
#' Psi <- matrix(c(1, 0.5), nrow = 2)
#' Phi <- array(0, dim = c(2, 2, 1))
#' Phi[, , 1] <- matrix(c(0.4, 0.1, 0, 0.3), 2, 2)
#' computeirf(Psi = Psi, Phi = Phi, H = 4, cum = FALSE)
#'
#' @export
computeirf <- function(Psi, Phi, H, cum) {
  Phi <- .as_numeric_array3(Phi, "Phi", allow_na = FALSE)
  if (dim(Phi)[1] != dim(Phi)[2]) {
    stop("Each Phi slice must be a square matrix (N x N).", call. = FALSE)
  }
  Psi <- .as_numeric_matrix(Psi, "Psi", nrow = dim(Phi)[1], allow_na = FALSE)
  H <- .check_integerish_scalar(H, "H", min = 1)

  N <- dim(Phi)[1]
  P <- dim(Phi)[3]
  R <- dim(Psi)[2]

  cum <- .check_cum(cum, N)

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

  irf
}
