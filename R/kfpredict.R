#' Deprecated: extract shocks via Kalman filter
#'
#' `kfpredict()` has been renamed to [predict_shocks()]. It is kept for
#' backward compatibility, issues a deprecation warning and will be removed in
#' a future release.
#'
#' @inheritParams predict_shocks
#'
#' @return See [predict_shocks()].
#'
#' @keywords internal
#' @export
kfpredict <- function(Sig, SigR, Psi, et, tol = sqrt(.Machine$double.eps), scale = TRUE) {
  .Deprecated("predict_shocks", package = "hetiv")
  predict_shocks(Sig = Sig, SigR = SigR, Psi = Psi, et = et, tol = tol, scale = scale)
}
