#' Plot impulse responses from one or several estimation approaches
#'
#' Produces a panel of IRF plots, one per variable-shock combination, for up
#' to three sets of impulse responses drawn in the same panels. Returns a list
#' of `ggplot` objects that can be arranged with [arrangeirf()] (or e.g.
#' `cowplot::plot_grid()`).
#'
#' @details Each set of impulse responses is plotted in one of three ways,
#'   depending on what is supplied:
#'   \itemize{
#'     \item Point estimate only: a line without confidence bands.
#'     \item Point estimate with standard errors: shaded normal bands
#'       `irf +/- qnorm((1 + ci) / 2) * se` for up to two confidence levels.
#'     \item Point estimate with bootstrap draws: shaded percentile bands for
#'       up to two confidence levels. With `boot_type = "recentred"` (the
#'       default) the quantiles are computed from the deviations of the draws
#'       from their bootstrap mean and added to the point estimate, so the
#'       (possibly asymmetric) band is centred on the reported point estimate.
#'       With `boot_type = "percentile"` the raw quantiles of the draws are
#'       used; such bands can be off-centre when the bootstrap distribution
#'       is biased relative to the point estimate.
#'   }
#'   Bootstrap draws take precedence over standard errors if both are present
#'   (as for [svariv()] with `nboot > 0`). Standard errors that are all zero
#'   suppress the bands.
#'
#' @param IRF Impulse responses to plot. One of
#'   \itemize{
#'     \item a numeric array (H x N x E, or H x N for a single shock) of point
#'       estimates, with numeric horizon labels in `dimnames()[[1]]`;
#'     \item an estimator result as returned by [hetiv()], [proxyiv()],
#'       [recols()] or [svariv()], i.e. a list with elements `irf` and
#'       optionally `se`, `boot` (array H x N x E x B of bootstrap draws) and
#'       `Method`;
#'     \item an (unnamed) list of up to three such objects, which are overlaid
#'       in the same panels.
#'   }
#' @param IRFse Array of standard errors with the same dimensions as `IRF`, or
#'   `NULL`. Only used when `IRF` is a numeric array (kept for backward
#'   compatibility). For list inputs, supply the standard errors as element
#'   `se` of each list.
#' @param HTick Integer. Step size for x-axis tick marks.
#' @param Labels Character vector of length N. Variable names used as panel
#'   titles.
#' @param ci Numeric vector of one or two confidence levels in (0, 1). One
#'   shaded band is drawn per level, with wider bands rendered at lower
#'   opacity. Defaults to `0.90` (a single 90% band); supply two levels, e.g.
#'   `c(0.90, 0.95)`, to draw two bands.
#' @param names Character vector with one name per set of impulse responses,
#'   used in the legend. Defaults to the `Method` element of each result if
#'   available and to `"IRF 1"`, `"IRF 2"`, ... otherwise.
#' @param colors Character vector of colours, one per set of impulse responses.
#'   Defaults to `c("steelblue", "darkred", "darkgreen")`.
#' @param boot_type Character. How bootstrap draws are turned into bands:
#'   `"recentred"` (default) or `"percentile"`. See Details. Ignored when no
#'   bootstrap draws are supplied.
#' @param legend Logical. If `TRUE`, a legend with `names` is added to each
#'   panel when more than one set of impulse responses is plotted. Default
#'   `FALSE`.
#'
#' @return A list of `ggplot` objects, one per variable-shock combination,
#'   ordered by shock (outer loop) then variable (inner loop).
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_ribbon geom_hline
#' @importFrom ggplot2 scale_x_continuous scale_colour_manual scale_fill_manual
#' @importFrom ggplot2 scale_linetype_manual theme_minimal theme element_text
#' @importFrom ggplot2 element_rect element_blank xlab ylab ggtitle element_line
#' @importFrom ggplot2 labs
#' @importFrom grid unit
#'
#' @examples
#' irf <- array(c(1, 0.5, 0.2, 0.1), dim = c(4, 1, 1))
#' dimnames(irf)[[1]] <- 0:3
#' se <- array(0.1, dim = dim(irf), dimnames = dimnames(irf))
#'
#' # One IRF with normal bands at two confidence levels
#' plotirf(irf, se, HTick = 1, Labels = "Output")
#'
#' # Two IRFs: the second without bands, the first with bootstrap draws
#' set.seed(1)
#' boot <- array(irf, dim = c(dim(irf), 200)) + rnorm(4 * 200, sd = 0.1)
#' plotirf(
#'   list(list(irf = irf, boot = boot), list(irf = irf * 0.8)),
#'   HTick = 1, Labels = "Output", ci = 0.90,
#'   names = c("Bootstrap", "Alternative"), legend = TRUE
#' )
#'
#' @export
plotirf <- function(IRF, IRFse = NULL, HTick, Labels, ci = 0.90,
                    names = NULL, colors = c("steelblue", "darkred", "darkgreen"),
                    boot_type = c("recentred", "percentile"), legend = FALSE) {
  boot_type <- match.arg(boot_type)
  # --- Normalise the input into a list of IRF specs ------------------------
  # A numeric array or a named list (estimator result) is a single IRF; an
  # unnamed list holds several IRFs to overlay
  is_single <- is.array(IRF) ||
    (is.list(IRF) && !is.null(names(IRF)) && all(nzchar(names(IRF))))
  if (is_single) {
    specs <- list(.as_irf_spec(IRF, "IRF", se = IRFse))
  } else if (is.list(IRF) && length(IRF) >= 1 && length(IRF) <= 3) {
    if (!is.null(IRFse)) {
      stop("IRFse cannot be used when IRF is a list of several IRFs.",
        call. = FALSE
      )
    }
    specs <- lapply(seq_along(IRF), function(k) {
      .as_irf_spec(IRF[[k]], paste0("IRF[[", k, "]]"))
    })
  } else {
    stop("IRF must be a numeric array, an estimator result, or a list of up ",
      "to three of these.", call. = FALSE)
  }
  nIRF <- length(specs)

  # All IRFs must share the horizon grid and the N x E layout
  dims <- dim(specs[[1]]$irf)
  HSeries <- specs[[1]]$HSeries
  for (k in seq_len(nIRF)) {
    if (!identical(dim(specs[[k]]$irf), dims) ||
      !isTRUE(all.equal(specs[[k]]$HSeries, HSeries))) {
      stop("All IRFs must have the same dimensions and horizon labels.",
        call. = FALSE
      )
    }
  }
  N <- dims[2]
  E <- dims[3]
  .validate_plot_common(HSeries, HTick, Labels, N)
  ci <- .check_ci_levels(ci)
  legend <- .check_logical_scalar(legend, "legend")

  # Series names, colours and line types
  if (is.null(names)) {
    names <- vapply(seq_len(nIRF), function(k) {
      if (!is.null(specs[[k]]$name)) specs[[k]]$name else paste("IRF", k)
    }, character(1))
  }
  if (!is.character(names) || length(names) != nIRF || anyDuplicated(names)) {
    stop("names must be a character vector with one distinct entry per IRF.",
      call. = FALSE
    )
  }
  if (!is.character(colors) || length(colors) < nIRF) {
    stop("colors must be a character vector with at least one colour per IRF.",
      call. = FALSE
    )
  }
  colors <- stats::setNames(colors[seq_len(nIRF)], names)
  ltypes <- stats::setNames(c("solid", "dotted", "dashed")[seq_len(nIRF)], names)
  # Alpha decreases for wider bands so narrower bands appear more opaque;
  # overlaid IRFs use lighter bands to keep all of them readable
  alpha_max <- if (nIRF == 1) 0.2 else 0.12
  alphas <- seq(alpha_max, alpha_max / 2, length.out = length(ci))

  # --- Confidence bands -----------------------------------------------------
  bands <- lapply(specs, .irf_bands, ci = ci, boot_type = boot_type)

  # --- Build one panel per shock and variable -------------------------------
  myGraphs <- list()
  n <- 1
  for (j in seq_len(E)) {
    for (i in seq_len(N)) {
      g1 <- ggplot2::ggplot(mapping = ggplot2::aes(x = Horizon)) +
        ggplot2::theme_minimal() +
        ggplot2::xlab("") +
        ggplot2::ylab("") +
        ggplot2::ggtitle(Labels[i]) +
        ggplot2::theme(plot.title = ggplot2::element_text(size = 10)) +
        ggplot2::scale_x_continuous(breaks = seq(HSeries[1], max(HSeries), HTick)) +
        ggplot2::geom_hline(yintercept = 0, linewidth = 0.2)

      # Ribbons first (behind), widest level first, then the lines on top
      for (k in seq_len(nIRF)) {
        if (is.null(bands[[k]])) next
        for (l in seq_along(ci)) {
          band_df <- data.frame(
            Horizon = HSeries,
            lower = unname(bands[[k]]$lower[, i, j, l]),
            upper = unname(bands[[k]]$upper[, i, j, l]),
            Series = names[k]
          )
          g1 <- g1 + ggplot2::geom_ribbon(
            data = band_df,
            ggplot2::aes(ymin = lower, ymax = upper, fill = Series),
            alpha = alphas[l]
          )
        }
      }
      for (k in seq_len(nIRF)) {
        line_df <- data.frame(
          Horizon = HSeries,
          Estimate = unname(specs[[k]]$irf[, i, j]),
          Series = names[k]
        )
        g1 <- g1 + ggplot2::geom_line(
          data = line_df,
          ggplot2::aes(y = Estimate, colour = Series, linetype = Series),
          linewidth = 0.7
        )
      }

      g1 <- g1 +
        ggplot2::scale_colour_manual(values = colors, breaks = names) +
        ggplot2::scale_fill_manual(values = colors, breaks = names) +
        ggplot2::scale_linetype_manual(values = ltypes, breaks = names) +
        ggplot2::labs(colour = NULL, fill = NULL, linetype = NULL) +
        ggplot2::theme(
          panel.grid = ggplot2::element_line(
            color = "gray", linewidth = 0.2, linetype = "dotted"
          ),
          panel.grid.minor = ggplot2::element_blank(),
          panel.border = ggplot2::element_rect(
            color = "black", fill = NA, linewidth = 0.2, linetype = "solid"
          )
        ) +
        ggplot2::theme(plot.margin = grid::unit(c(0.1, 0.1, 0.1, 0.1), "cm"))

      if (legend && nIRF > 1) {
        g1 <- g1 + ggplot2::theme(
          legend.position = "bottom",
          legend.text = ggplot2::element_text(size = 8),
          legend.key.height = grid::unit(0.3, "cm")
        )
      } else {
        g1 <- g1 + ggplot2::theme(legend.position = "none")
      }

      myGraphs[[n]] <- g1
      n <- n + 1
    }
  }

  return(myGraphs)
}
