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
#' @export
filllinear <- function(DF, gap) {
  for (Stats in colnames(DF)) {
    DF[, Stats] <- zoo::na.approx(DF[, Stats], na.rm = FALSE, maxgap = gap)
  }
  return(DF)
}
