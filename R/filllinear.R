#' Linearly interpolate gaps in a data frame
#'
#' Fills missing values (`NA`) in every column of a data frame using linear
#' interpolation, up to a maximum of `gap` consecutive missing observations.
#' Gaps longer than `gap` are left as `NA`.
#'
#' @param DF A data frame or matrix with numeric columns.
#' @param gap Integer. Maximum number of consecutive `NA`s to fill.
#'
#' @return The data frame `DF` with short gaps linearly interpolated.
#'
#' @importFrom zoo na.approx
#'
#' @examples
#' filllinear(data.frame(x = c(1, NA, 3), y = c(2, NA, 4)), gap = 1)
#'
#' @export
filllinear <- function(DF, gap) {
  if (!is.data.frame(DF) && !is.matrix(DF)) {
    stop("DF must be a data frame or matrix.", call. = FALSE)
  }
  if (!all(vapply(as.data.frame(DF), is.numeric, logical(1)))) {
    stop("All columns in DF must be numeric.", call. = FALSE)
  }
  gap <- .check_integerish_scalar(gap, "gap", min = 0)
  for (Stats in colnames(DF)) {
    DF[, Stats] <- zoo::na.approx(DF[, Stats], na.rm = FALSE, maxgap = gap)
  }
  return(DF)
}
