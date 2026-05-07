#' Plot and compare impulse responses from two estimation approaches
#'
#' Produces a panel of IRF plots comparing two sets of impulse responses
#' side-by-side. Each panel shows the point estimates and 90% confidence
#' bands for both approaches. Returns a list of `ggplot` objects (one per
#' variable-shock combination) that can be arranged with e.g. `gridExtra::grid.arrange()`.
#'
#' @param IRF1 Array (H x N x E) of impulse responses for the first approach.
#'   Dimnames on the first dimension must be the horizon labels (numeric).
#' @param IRF1se Array (H x N x E) of standard errors for the first approach.
#' @param IRF2 Array (H x N x E) of impulse responses for the second approach.
#' @param IRF2se Array (H x N x E) of standard errors for the second approach.
#' @param HTick Integer. Step size for x-axis tick marks.
#' @param Labels Character vector of length N. Variable names used as panel titles.
#' @param ci Numeric scalar in (0, 1). Confidence level for the shaded bands.
#'   Defaults to `0.90` (90% CI).
#'
#' @return A list of `ggplot` objects, one per variable-shock combination,
#'   ordered by shock (outer loop) then variable (inner loop).
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_hline scale_x_continuous theme_minimal theme element_text element_rect element_blank unit xlab ylab ggtitle element_line
#'
#' @export
plot2irf <- function(IRF1, IRF1se, IRF2, IRF2se, HTick, Labels, ci = 0.90){

  myGraphs <- list()
  noDims   <- length(dim(IRF1))
  HNum     <- dim(IRF1)[1]
  HSeries  <- as.numeric(dimnames(IRF1)[[1]])

  N <- dim(IRF1)[2]
  n <- 1

  # Handle 2-dimensional input (single shock, no E dimension)
  if(noDims == 3){
    E <- dim(IRF1)[3]
  }else{
    E <- 1
    dim(IRF1)   <- c(HNum, N, 1)
    dim(IRF2)   <- c(HNum, N, 1)
    dim(IRF1se) <- c(HNum, N, 1)
    dim(IRF2se) <- c(HNum, N, 1)
  }

  z         <- qnorm((1 + ci) / 2)
  upper1    <- IRF1 + z * IRF1se
  lower1    <- IRF1 - z * IRF1se
  upper2    <- IRF2 + z * IRF2se
  lower2    <- IRF2 - z * IRF2se

  for(j in 1:E){
    for(i in 1:N){

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
        ggplot2::theme(legend.position = "none") +
        ggplot2::geom_line(ggplot2::aes(y = IRF1), colour = "steelblue", linewidth = 0.7) +
        ggplot2::geom_line(ggplot2::aes(y = IRF2), colour = "darkred",   linewidth = 0.7, linetype = "dotted") +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower1, ymax = upper1), fill = "steelblue", alpha = 0.1) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = lower2, ymax = upper2), fill = "darkred",   alpha = 0.1) +
        ggplot2::theme(
          panel.grid       = ggplot2::element_line(color = "gray", linewidth = 0.2, linetype = "dotted"),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border     = ggplot2::element_rect(color = "black", fill = NA, linewidth = 0.2, linetype = "solid")
        ) +
        ggplot2::theme(plot.margin = ggplot2::unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

      myGraphs[[n]] <- g1
      n <- n + 1
    }
  }

  return(myGraphs)
}
