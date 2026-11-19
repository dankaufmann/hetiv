#' Squared correlation of the predicted shocks with the true shocks
#'
#' Computes the population squared correlation between the shocks predicted
#' by [predict_shocks()] and the true underlying shocks, as a measure of the
#' precision of the prediction (Burri and Kaufmann, 2026a).
#'
#' @details Let the residuals on event days be `e_t = Psi eps_t + v_t`, with
#'   unit-variance shocks `eps_t` and non-policy noise `v_t` with covariance
#'   `SigR`. The linear prediction of the shocks has covariance with the true
#'   shocks `S (I + S)^{-1}`, with the signal-to-noise matrix
#'   `S = Psi' SigR^{-1} Psi`. The diagonal of this matrix is the share of the
#'   shock variance explained by the prediction (one minus the relative mean
#'   squared error), that is, the squared correlation between the predicted
#'   and the true shock.
#'
#'   The correlation itself does not depend on the units of the shocks, but
#'   the formula does: `S` is the signal-to-noise ratio only if `Psi` is the
#'   response to a one-standard-deviation shock. The estimators return `Psi`
#'   normalised to a unit impact response instead, which leaves the variance
#'   of the shocks out of `S`.
#'   If `Sig` is supplied, the columns of `Psi` are therefore rescaled by the
#'   implied standard deviations of the shocks, recovered from `Sig - SigR`
#'   in the same way as in [predict_shocks()]. If `Sig` is `NULL`, `Psi` is
#'   used as is.
#'
#'   The measure is a population quantity evaluated at the estimated
#'   parameters. It does not account for estimation error in `Psi`, `Sig` and
#'   `SigR`, and it is informative only if the responses in `Psi` are
#'   consistently estimated, that is, if the identifying assumptions of the
#'   estimator hold.
#'
#'   The measure also requires a nonsingular covariance `SigR`. If fewer
#'   non-policy shocks than variables drive the data, `SigR` is (nearly)
#'   singular, the implied signal-to-noise ratio becomes very large and the
#'   measure approaches one whatever the actual precision. A warning is
#'   issued in this case.
#'
#' @param Psi Numeric matrix (M x E) of responses of the residuals to the E
#'   structural shocks (e.g. `Psi` from [hetiv()]).
#' @param SigR Numeric matrix (M x M). Covariance matrix of the residuals on
#'   control (non-event) days.
#' @param Sig Numeric matrix (M x M) or `NULL`. Covariance matrix of the
#'   residuals on event days, used to rescale `Psi` to one-standard-deviation
#'   shocks. If `NULL` (default), `Psi` is assumed to be in these units
#'   already.
#' @param squared Logical. If `TRUE` (default), the squared correlation is
#'   returned. If `FALSE`, its square root (the correlation).
#'
#' @return A numeric vector of length E with the squared correlation between
#'   each predicted shock and the true shock (the diagonal of
#'   `S (I + S)^{-1}`), or the correlation if `squared = FALSE`.
#'
#' @references
#' Burri, M. and Kaufmann, D. (2026a). Measuring monetary policy shocks.
#' IRENE Working Papers 24-03, IRENE Institute of Economic Research,
#' University of Neuchatel.
#'
#' @seealso [predict_shocks()]
#'
#' @importFrom MASS ginv
#'
#' @examples
#' # One shock with unit variance, two variables
#' Psi <- matrix(c(1, 0.5), nrow = 2)
#' SigR <- diag(c(0.5, 0.5))
#' prediction_corr(Psi, SigR)
#'
#' # Psi normalised to unit impact; the shock has variance 4
#' Sig <- SigR + 4 * Psi %*% t(Psi)
#' prediction_corr(Psi, SigR, Sig = Sig)
#'
#' @export
prediction_corr <- function(Psi, SigR, Sig = NULL, squared = TRUE) {
  SigR <- .as_numeric_matrix(SigR, "SigR", square = TRUE, allow_na = FALSE)
  Psi <- .as_numeric_matrix(Psi, "Psi", nrow = nrow(SigR), allow_na = FALSE)
  squared <- .check_logical_scalar(squared, "squared")
  E <- ncol(Psi)

  # Rescale the unit-impact responses to one-standard-deviation shocks
  if (!is.null(Sig)) {
    Sig <- .as_numeric_matrix(Sig, "Sig", square = TRUE, allow_na = FALSE)
    if (!identical(dim(Sig), dim(SigR))) {
      stop("Sig and SigR must have the same dimensions.", call. = FALSE)
    }
    sig <- .shock_variances(Sig, SigR, Psi)
    if (any(sig < 0)) {
      warning(
        "Some estimated shock variances are negative (weak heteroskedasticity); ",
        "using absolute values for scaling."
      )
    }
    Psi <- Psi %*% diag(sqrt(abs(sig)), nrow = E, ncol = E)
  }

  # Flag a (nearly) singular noise covariance: the signal-to-noise ratio then
  # becomes very large and the measure approaches one whatever the actual
  # precision. This happens when fewer non-policy shocks than variables drive
  # the data. The threshold is stricter than in predict_shocks(), because
  # SigR enters only through its inverse.
  ev <- eigen(stats::cov2cor(SigR), symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < 1e-3 * max(ev)) {
    warning(
      "SigR is (nearly) singular: the smallest eigenvalue of its correlation ",
      "matrix is ", format(min(ev) / max(ev), digits = 2), " relative to the ",
      "largest. The implied squared correlation is then close to one and ",
      "overstates the precision of the prediction. This happens when fewer ",
      "non-policy shocks than variables drive the data.",
      call. = FALSE
    )
  }

  # Signal-to-noise matrix and covariance of predicted and true shocks
  S <- t(Psi) %*% MASS::ginv(SigR) %*% Psi
  pred_cov <- S %*% MASS::ginv(diag(E) + S)

  # Diagonal: squared correlation of each predicted shock with the true shock
  r2 <- pmin(pmax(diag(pred_cov), 0), 1)
  if (squared) r2 else sqrt(r2)
}
