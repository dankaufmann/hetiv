#' Arrange impulse response plots into one figure per shock
#'
#' Arranges a list of IRF panels, as returned by [plotirf()] or [plot2irf()],
#' into a grid with one figure per shock dimension and a common x-axis label.
#' Optionally saves each figure as a PDF.
#'
#' @param myGraphs List of `ggplot` objects of length `E * N`, ordered by
#'   shock (outer loop) and then by variable (inner loop), as returned by
#'   [plotirf()] and [plot2irf()].
#' @param E Integer. Number of shock dimensions (one figure per shock).
#' @param N Integer. Number of variables (panels per figure).
#' @param nRows Integer. Number of rows of panels in each figure. The number
#'   of columns is `ceiling(N / nRows)`.
#' @param xLab Character. Common x-axis label printed below the panels.
#' @param file Character or `NA`. If not `NA`, file path stem; figure `e` is
#'   saved as `paste0(file, "_", e, ".pdf")`. If `NA` (default), nothing is
#'   saved.
#' @param figScaleW Numeric. Width of one panel column in inches.
#' @param figScaleH Numeric. Height of one panel row in inches. The saved
#'   figure has height `0.9 * nRows * figScaleH`.
#'
#' @return Invisibly, a list of length `E` of arranged figures (`gtable`
#'   objects), which can be redrawn with [grid::grid.draw()]. Each figure is
#'   also drawn on the current graphics device as a side effect.
#'
#' @importFrom gridExtra grid.arrange
#' @importFrom grid textGrob gpar
#' @importFrom ggplot2 ggsave
#'
#' @examples
#' irf <- array(c(1, 0.5, 0.2, 0.1, 0, 0.3, 0.2, 0.1), dim = c(4, 2, 1))
#' dimnames(irf)[[1]] <- 0:3
#' se <- array(0.1, dim = dim(irf), dimnames = dimnames(irf))
#' g <- plotirf(irf, se, HTick = 1, Labels = c("Output", "Prices"))
#' figs <- arrangeirf(g, E = 1, N = 2, nRows = 1)
#'
#' \donttest{
#' # Save to PDF (file is written as <stem>_1.pdf)
#' stem <- file.path(tempdir(), "irf")
#' arrangeirf(g, E = 1, N = 2, nRows = 1, file = stem)
#' }
#'
#' @export
arrangeirf <- function(myGraphs, E, N, nRows, xLab = "Horizon", file = NA,
                       figScaleW = 3.5, figScaleH = 2.5) {
  # Validate inputs
  E <- .check_integerish_scalar(E, "E", min = 1)
  N <- .check_integerish_scalar(N, "N", min = 1)
  nRows <- .check_integerish_scalar(nRows, "nRows", min = 1, max = N)
  figScaleW <- .check_numeric_scalar(figScaleW, "figScaleW", min = 0)
  figScaleH <- .check_numeric_scalar(figScaleH, "figScaleH", min = 0)
  if (!is.list(myGraphs) || length(myGraphs) != E * N) {
    stop("myGraphs must be a list of length E * N.", call. = FALSE)
  }
  if (!(length(file) == 1 && (is.na(file) || is.character(file)))) {
    stop("file must be NA or a single character string.", call. = FALSE)
  }

  nCols <- ceiling(N / nRows)
  myXLab <- grid::textGrob(xLab, gp = grid::gpar(fontsize = 9.5), vjust = -1)
  figs <- vector("list", E)

  # One figure per shock dimension
  for (e in 1:E) {
    figs[[e]] <- gridExtra::grid.arrange(
      grobs = myGraphs[((e - 1) * N + 1):(e * N)],
      nrow = nRows, ncol = nCols, bottom = myXLab
    )

    if (!is.na(file)) {
      ggplot2::ggsave(
        filename = paste0(file, "_", e, ".pdf"), plot = figs[[e]],
        width = nCols * figScaleW, height = nRows * 0.9 * figScaleH
      )
    }
  }

  invisible(figs)
}
