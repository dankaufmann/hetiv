#' Plot p-values across horizons
#'
#' Produces a panel of p-value plots across IRF horizons, with optional
#' horizontal reference lines at user-supplied significance levels. Returns a
#' list of `ggplot` objects (one per variable-shock combination) that can be
#' arranged with e.g. `cowplot::plot_grid()`.
#'
#' @param pvals Array (H x N x E) of p-values. Dimnames on the first dimension
#'   must be the horizon labels (numeric).
#' @param HTick Integer. Step size for x-axis tick marks.
#' @param Labels Character vector of length N. Variable names used as panel titles.
#' @param sigLevels Numeric vector of significance levels at which to draw
#'   horizontal reference lines. The first value is drawn as a solid line, the
#'   second as dashed, and any further values as dotted. Defaults to
#'   `c(0.05, 0.10)`.
#'
#' @return A list of `ggplot` objects, one per variable-shock combination,
#'   ordered by shock (outer loop) then variable (inner loop).
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_hline scale_x_continuous
#' @importFrom ggplot2 theme_minimal theme element_text element_rect
#' @importFrom ggplot2 element_blank unit xlab ylab ggtitle element_line
#'
#' @examples
#' pvals <- array(c(0.2, 0.08, 0.04), dim = c(3, 1, 1))
#' dimnames(pvals)[[1]] <- 0:2
#' plotpval(pvals, HTick = 1, Labels = "Output")
#'
#' @export
plotpval <- function(pvals, HTick, Labels, sigLevels = c(0.05, 0.10)) {
  pvals <- .validate_irf_array(pvals, "pvals")
  if (any(pvals < 0 | pvals > 1, na.rm = TRUE)) {
    stop("pvals must contain values between 0 and 1.", call. = FALSE)
  }
  if (
    !is.numeric(sigLevels) || length(sigLevels) < 1 || anyNA(sigLevels) ||
      any(sigLevels < 0 | sigLevels > 1)
  ) {
    stop("sigLevels must contain values between 0 and 1.", call. = FALSE)
  }
  linetypes <- c("solid", "dashed", "dotted")

  myGraphs <- list()
  noDims <- length(dim(pvals))
  HNum <- dim(pvals)[1]
  HSeries <- .check_horizon_labels(pvals, "pvals")

  N <- dim(pvals)[2]
  .validate_plot_common(HSeries, HTick, Labels, N)
  n <- 1

  # Handle 2-dimensional input (single shock, no E dimension)
  if (noDims == 3) {
    E <- dim(pvals)[3]
  } else {
    E <- 1
    dim(pvals) <- c(HNum, N, 1)
  }

  for (j in 1:E) {
    for (i in 1:N) {
      myPval <- data.frame(Horizon = HSeries, pval = pvals[, i, j])

      g1 <- ggplot2::ggplot(myPval, ggplot2::aes(x = Horizon)) +
        ggplot2::theme_minimal() +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::ggtitle(Labels[i]) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 10)) +
        ggplot2::scale_x_continuous(breaks = seq(HSeries[1], max(HSeries), HTick)) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::geom_line(ggplot2::aes(y = pval), colour = "steelblue", linewidth = 0.7)

      # Draw one reference line per significance level
      for (k in seq_along(sigLevels)) {
        lt <- linetypes[min(k, length(linetypes))]
        g1 <- g1 + ggplot2::geom_hline(
          yintercept = sigLevels[k],
          colour     = "black",
          linetype   = lt,
          linewidth  = 0.4
        )
      }

      g1 <- g1 +
        ggplot2::theme(
          panel.grid = ggplot2::element_line(
            color = "gray", linewidth = 0.2, linetype = "dotted"
          ),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(
            color = "black", fill = NA, linewidth = 0.2, linetype = "solid"
          )
        ) +
        ggplot2::theme(plot.margin = ggplot2::unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

      myGraphs[[n]] <- g1
      n <- n + 1
    }
  }

  return(myGraphs)
}
