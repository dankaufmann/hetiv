#' Plot and compare impulse responses from two estimation approaches
#'
#' Produces a panel of IRF plots comparing two sets of impulse responses with
#' normal confidence bands at level `ci`. This function is kept for backward
#' compatibility and is a thin wrapper around [plotirf()], which accepts up to
#' three sets of impulse responses as well as bootstrap draws.
#'
#' @param IRF1 Array (H x N x E) of impulse responses for the first approach.
#'   Dimnames on the first dimension must be the horizon labels (numeric).
#' @param IRF1se Array (H x N x E) of standard errors for the first approach
#'   (set to 0 to suppress confidence intervals).
#' @param IRF2 Array (H x N x E) of impulse responses for the second approach.
#' @param IRF2se Array (H x N x E) of standard errors for the second approach
#'   (set to 0 to suppress confidence intervals).
#' @param HTick Integer. Step size for x-axis tick marks.
#' @param Labels Character vector of length N. Variable names used as panel titles.
#' @param ci Numeric scalar in (0, 1). Confidence level for the shaded bands.
#'   Defaults to `0.90` (90% CI).
#'
#' @return A list of `ggplot` objects, one per variable-shock combination,
#'   ordered by shock (outer loop) then variable (inner loop).
#'
#' @seealso [plotirf()]
#'
#' @examples
#' irf1 <- array(c(1, 0.5, 0.2, 0.1), dim = c(4, 1, 1))
#' dimnames(irf1)[[1]] <- 0:3
#' irf2 <- irf1 * 0.8
#' se <- array(0.1, dim = dim(irf1), dimnames = dimnames(irf1))
#' plot2irf(irf1, se, irf2, se * 0, HTick = 1, Labels = "Output")
#'
#' @export
plot2irf <- function(IRF1, IRF1se, IRF2, IRF2se, HTick, Labels, ci = 0.90) {
  IRF1 <- .validate_irf_array(IRF1, "IRF1")
  IRF2 <- .validate_irf_array(IRF2, "IRF2")
  IRF1se <- .validate_irf_array(IRF1se, "IRF1se")
  IRF2se <- .validate_irf_array(IRF2se, "IRF2se")
  if (
    !identical(dim(IRF1), dim(IRF2)) ||
      !identical(dim(IRF1), dim(IRF1se)) ||
      !identical(dim(IRF1), dim(IRF2se))
  ) {
    stop("IRF1, IRF1se, IRF2, and IRF2se must have matching dimensions.",
      call. = FALSE
    )
  }
  .check_numeric_scalar(ci, "ci",
    min = .Machine$double.eps,
    max = 1 - .Machine$double.eps
  )

  plotirf(
    IRF = list(
      list(irf = IRF1, se = IRF1se),
      list(irf = IRF2, se = IRF2se)
    ),
    HTick = HTick, Labels = Labels, ci = ci,
    names = c("Approach 1", "Approach 2"),
    colors = c("steelblue", "darkred"),
    legend = FALSE
  )
}
