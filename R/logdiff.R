#' Log-difference transformation
#'
#' Computes the first log-difference of a numeric vector or time series,
#' i.e. `log(x_t) - log(x_{t-1})`. Gaps in dates are ignored; the difference
#' is always taken with respect to the immediately preceding observation.
#'
#' @param TS A numeric vector or time series.
#'
#' @return A numeric vector of the same length as `TS`, with `NA` in the
#'   first position.
#'
#' @examples
#' logdiff(c(1, 2, 4))
#'
#' @export
logdiff <- function(TS) {
  if (!is.numeric(TS)) stop("TS must be numeric.", call. = FALSE)
  if (any(TS <= 0, na.rm = TRUE))
    warning("Non-positive values in TS; log-differences will contain NaN.")
  log(TS) - dplyr::lag(log(TS), 1)
}
