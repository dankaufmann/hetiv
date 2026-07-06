#' Plot and compare impulse responses from two estimation approaches
#'
#' Produces a panel of IRF plots comparing two sets of impulse responses.
#' Each panel shows the point estimates and confidence bands
#' at level `ci` for both approaches. Returns a list of `ggplot` objects
#' (one per variable-shock combination) that can be arranged with e.g.
#' `cowplot::plot_grid()`.
#'
#' @param IRF1 Array (H x N x E) of impulse responses for the first approach.
#'   Dimnames on the first dimension must be the horizon labels (numeric).
#' @param IRF1se Array (H x N x E) of standard errors for the first approach
#'   (set to 0 to suppress confidence intervals).
#' @param IRF2 Array (H x N x E) of impulse responses for the second approach.
#' @param IRF2se Array (H x N x E) of standard errors for the second approach
#'   (set to 0 to suppress confidence intervals).
#' @param HTick Integer. Step size for x-axis tick marks.
#' @param Labels Character vector of length N. Variable names used as panel titles.
#' @param ci Numeric scalar in (0, 1). Confidence level for the shaded bands.
#'   Defaults to `0.90` (90% CI).
#'
#' @return A list of `ggplot` objects, one per variable-shock combination,
#'   ordered by shock (outer loop) then variable (inner loop).
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_hline
#' @importFrom ggplot2 scale_x_continuous scale_colour_manual
#' @importFrom ggplot2 scale_linetype_manual theme_minimal theme element_text
#' @importFrom ggplot2 element_rect element_blank xlab ylab ggtitle
#' @importFrom ggplot2 element_line labs
#' @importFrom grid unit
#'
#' @examples
#' irf1 <- array(c(1, 0.5, 0.2, 0.1), dim = c(4, 1, 1))
#' dimnames(irf1)[[1]] <- 0:3
#' irf2 <- irf1 * 0.8
#' se <- array(0.1, dim = dim(irf1), dimnames = dimnames(irf1))
#' plot2irf(irf1, se, irf2, se * 0, HTick = 1, Labels = "Output")
#'
#' @export
plot2irf <- function(IRF1, IRF1se, IRF2, IRF2se, HTick, Labels, ci = 0.90) {
  IRF1 <- .validate_irf_array(IRF1, "IRF1")
  IRF2 <- .validate_irf_array(IRF2, "IRF2")
  IRF1se <- .validate_irf_array(IRF1se, "IRF1se")
  IRF2se <- .validate_irf_array(IRF2se, "IRF2se")
  if (
    !identical(dim(IRF1), dim(IRF2)) ||
      !identical(dim(IRF1), dim(IRF1se)) ||
      !identical(dim(IRF1), dim(IRF2se))
  ) {
    stop("IRF1, IRF1se, IRF2, and IRF2se must have matching dimensions.",
      call. = FALSE
    )
  }
  myGraphs <- list()
  noDims <- length(dim(IRF1))
  HNum <- dim(IRF1)[1]
  HSeries <- .check_horizon_labels(IRF1, "IRF1")

  N <- dim(IRF1)[2]
  .validate_plot_common(HSeries, HTick, Labels, N)
  .check_numeric_scalar(ci, "ci",
    min = .Machine$double.eps,
    max = 1 - .Machine$double.eps
  )
  n <- 1

  # Handle 2-dimensional input (single shock, no E dimension)
  if (noDims == 3) {
    E <- dim(IRF1)[3]
  } else {
    E <- 1
    dim(IRF1) <- c(HNum, N, 1)
    dim(IRF2) <- c(HNum, N, 1)
    dim(IRF1se) <- c(HNum, N, 1)
    dim(IRF2se) <- c(HNum, N, 1)
  }

  z <- qnorm((1 + ci) / 2)
  upper1 <- IRF1 + z * IRF1se
  lower1 <- IRF1 - z * IRF1se
  upper2 <- IRF2 + z * IRF2se
  lower2 <- IRF2 - z * IRF2se

  for (j in 1:E) {
    for (i in 1:N) {
      myIRF <- data.frame(
        Horizon = HSeries,
        IRF1    = IRF1[, i, j],
        IRF2    = IRF2[, i, j],
        upper1  = upper1[, i, j],
        lower1  = lower1[, i, j],
        upper2  = upper2[, i, j],
        lower2  = lower2[, i, j]
      )

      g1 <- ggplot2::ggplot(myIRF, ggplot2::aes(x = Horizon)) +
        ggplot2::theme_minimal() +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::ggtitle(Labels[i]) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 10)) +
        ggplot2::scale_x_continuous(breaks = seq(HSeries[1], max(HSeries), HTick)) +
        ggplot2::geom_hline(yintercept = 0, linewidth = 0.2) +
        ggplot2::geom_line(
          ggplot2::aes(y = IRF1, colour = "Approach 1", linetype = "Approach 1"),
          linewidth = 0.7
        ) +
        ggplot2::geom_line(
          ggplot2::aes(y = IRF2, colour = "Approach 2", linetype = "Approach 2"),
          linewidth = 0.7
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lower1, ymax = upper1),
          fill = "steelblue", alpha = 0.1
        ) +
        ggplot2::geom_ribbon(
          ggplot2::aes(ymin = lower2, ymax = upper2),
          fill = "darkred", alpha = 0.1
        ) +
        ggplot2::scale_colour_manual(
          values = c("Approach 1" = "steelblue", "Approach 2" = "darkred")
        ) +
        ggplot2::scale_linetype_manual(
          values = c("Approach 1" = "solid", "Approach 2" = "dotted")
        ) +
        ggplot2::labs(colour = NULL, linetype = NULL) +
        ggplot2::theme(
          panel.grid = ggplot2::element_line(
            color = "gray", linewidth = 0.2, linetype = "dotted"
          ),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(
            color = "black", fill = NA, linewidth = 0.2, linetype = "solid"
          )
        ) +
        ggplot2::theme(plot.margin = grid::unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
        ggplot2::theme(legend.position = "none")


      myGraphs[[n]] <- g1
      n <- n + 1
    }
  }

  return(myGraphs)
}
