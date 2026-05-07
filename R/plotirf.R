#' Plot impulse responses from a single estimation approach
#'
#' Produces a panel of IRF plots for a single set of impulse responses,
#' with optional 90% and 95% confidence bands. Returns a list of `ggplot`
#' objects (one per variable-shock combination) that can be arranged with
#' e.g. `cowplot::plot_grid()`.
#'
#' @param IRFest Array (H x N x E) of impulse response point estimates.
#'   Dimnames on the first dimension must be the horizon labels (numeric).
#' @param IRFse Array (H x N x E) of standard errors, or `NULL` to suppress
#'   confidence bands.
#' @param HTick Integer. Step size for x-axis tick marks.
#' @param Labels Character vector of length N. Variable names used as panel titles.
#' @param ci Numeric vector of confidence levels in (0, 1). One shaded band is
#'   drawn per level, with wider bands rendered at lower opacity. Defaults to
#'   `c(0.90, 0.95)`.
#'
#' @return A list of `ggplot` objects, one per variable-shock combination,
#'   ordered by shock (outer loop) then variable (inner loop).
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_hline scale_x_continuous theme_minimal theme element_text element_rect element_blank unit xlab ylab ggtitle element_line
#'
#' @export
plotirf <- function(IRFest, IRFse = NULL, HTick, Labels, ci = c(0.90, 0.95)){

  myGraphs <- list()
  noDims   <- length(dim(IRFest))
  HNum     <- dim(IRFest)[1]
  HSeries  <- as.numeric(dimnames(IRFest)[[1]])

  N <- dim(IRFest)[2]
  n <- 1

  # Handle 2-dimensional input (single shock, no E dimension)
  if(noDims == 3){
    E <- dim(IRFest)[3]
  }else{
    E <- 1
    dim(IRFest) <- c(HNum, N, 1)
  }

  # Compute confidence bands for each level in ci (sorted widest first for layering)
  has_se <- is.array(IRFse)
  ci     <- sort(ci, decreasing = TRUE)   # widest band drawn first (behind)
  if(has_se){
    if(noDims != 3) dim(IRFse) <- c(HNum, N, 1)
    bands <- lapply(ci, function(level){
      z <- qnorm((1 + level) / 2)
      list(upper = IRFest + z * IRFse, lower = IRFest - z * IRFse)
    })
  }

  for(j in 1:E){
    for(i in 1:N){

      myIRF <- data.frame(Horizon = HSeries, Estimate = IRFest[, i, j])
      if(has_se){
        for(k in seq_along(ci)){
          myIRF[[paste0("upper", k)]] <- bands[[k]]$upper[, i, j]
          myIRF[[paste0("lower", k)]] <- bands[[k]]$lower[, i, j]
        }
      }

      g1 <- ggplot2::ggplot(myIRF, ggplot2::aes(x = Horizon)) +
        ggplot2::theme_minimal() +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::ggtitle(Labels[i]) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 10)) +
        ggplot2::scale_x_continuous(breaks = seq(HSeries[1], max(HSeries), HTick)) +
        ggplot2::geom_hline(yintercept = 0, linewidth = 0.2) +
        ggplot2::theme(legend.position = "none") +
        ggplot2::geom_line(ggplot2::aes(y = Estimate), colour = "steelblue", linewidth = 0.7)

      if(has_se){
        # Alpha decreases for wider bands so narrower bands appear more opaque
        alphas <- seq(0.15, 0.10, length.out = length(ci))
        for(k in seq_along(ci)){
          g1 <- g1 + ggplot2::geom_ribbon(
            ggplot2::aes(ymin = .data[[paste0("lower", k)]], ymax = .data[[paste0("upper", k)]]),
            fill = "steelblue", alpha = alphas[k]
          )
        }
      }

      g1 <- g1 +
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
