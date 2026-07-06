#' First-difference transformation
#'
#' Computes the first difference of a numeric vector or time series,
#' i.e. `x_t - x_{t-1}`. Gaps in dates are ignored; the difference is always
#' taken with respect to the immediately preceding observation.
#'
#' @param TS A numeric vector or time series.
#'
#' @return A numeric vector of the same length as `TS`, with `NA` in the
#'   first position.
#'
#' @examples
#' firstdiff(c(1, 3, 6))
#'
#' @export
firstdiff <- function(TS) {
  if (!is.numeric(TS)) stop("TS must be numeric.", call. = FALSE)
  TS - dplyr::lag(TS, 1)
}
