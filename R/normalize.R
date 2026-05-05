#' Normalize a time series
#'
#' Standardizes a numeric vector to zero mean and unit standard deviation
#' (z-score normalization), ignoring missing values.
#'
#' @param x A numeric vector or time series.
#'
#' @return A numeric vector of the same length as `x`, with mean 0 and
#'   standard deviation 1.
#'
#' @examples
#' normalize(c(1, 2, 3, 4, 5))
#'
#' @export
normalize <- function(x) {
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}
