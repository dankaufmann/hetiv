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
#' @export
logdiff <- function(TS) {
  log(TS) - dplyr::lag(log(TS), 1)
}
