#' Rescale a series to match the mean and variance of another series
#'
#' Applies the linear transformation to `y` that matches the mean and the
#' standard deviation of `x` on their overlapping (jointly observed) sample.
#' Useful to compare measures of the same shock that are expressed in
#' different units, for example a predicted shock and a proxy.
#'
#' @param x Numeric vector. Reference series (mean and variance target).
#' @param y Numeric vector of the same length as `x`. Series to rescale.
#'
#' @return Numeric vector of the same length as `y`: `a + b * y`, with `a` and
#'   `b` chosen such that the result has the mean and standard deviation of
#'   `x` on the overlapping sample. Missing values of `y` remain missing.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(50, mean = 2, sd = 3)
#' y <- 0.1 * x + rnorm(50, sd = 0.1)
#' y2 <- rescale_to_match(x, y)
#' c(mean(y2), sd(y2))
#'
#' @export
rescale_to_match <- function(x, y) {
  if (!is.numeric(x) || !is.numeric(y) || length(x) != length(y)) {
    stop("x and y must be numeric vectors of the same length.", call. = FALSE)
  }

  overlap <- !is.na(x) & !is.na(y)
  if (sum(overlap) < 2L) stop("Need at least 2 overlapping observations.", call. = FALSE)

  mu_x <- mean(x[overlap])
  mu_y <- mean(y[overlap])
  sd_x <- stats::sd(x[overlap])
  sd_y <- stats::sd(y[overlap])
  if (sd_y == 0) {
    stop("`y` has zero variance on the overlap; cannot rescale.", call. = FALSE)
  }

  b <- sd_x / sd_y
  a <- mu_x - b * mu_y
  a + b * y
}
