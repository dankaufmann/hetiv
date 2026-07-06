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
  if (!is.numeric(x)) stop("x must be numeric.", call. = FALSE)
  if (is.na(sd(x, na.rm = TRUE)) || sd(x, na.rm = TRUE) == 0) {
    stop("x must contain at least two non-missing, non-constant values.",
         call. = FALSE)
  }
  (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}
