#' Plot one or several series in one panel
#'
#' Line plot of any number of series against a common x-axis variable (for
#' example time, or a parameter varied in a simulation), with one colour and
#' line type per series and a single legend. Optionally standardises the
#' series and, for two series, reports their correlation in the panel.
#'
#' @param data A data frame with one column for the x-axis and one column per
#'   series. Alternatively, a time series object (or a named list of time
#'   series objects) supported by the `tsbox` package, which is converted to
#'   a data frame with the x-axis column `time`; this requires `tsbox` to be
#'   installed.
#' @param x_var Character. Name of the x-axis column. Defaults to the first
#'   column of `data`.
#' @param series Character vector with the names of the columns to plot. If
#'   the vector is named, the names are the column names and the values are
#'   the legend labels. The order determines the legend order. Defaults to
#'   all columns except `x_var`.
#' @param colors,linetypes Vectors with at least one entry per series.
#' @param x_lab,y_lab Character or `NULL`. Axis labels.
#' @param standardize Logical. If `TRUE`, every series is demeaned and divided
#'   by its standard deviation before plotting. Default `FALSE`.
#' @param show_cor Logical. If `TRUE`, the correlation of the series on their
#'   common sample is printed in the lower left corner of the panel. Requires
#'   exactly two series. Default `FALSE`.
#' @param show_n Logical. If `TRUE` (default), the number of common
#'   observations is added to the correlation label.
#' @param digits Integer. Number of digits of the correlation label.
#' @param legend_position Character such as `"bottom"` (default) or `"none"`,
#'   or a numeric vector of length 2 with the position inside the panel in
#'   relative coordinates.
#' @param linewidth Numeric. Line width.
#' @param xlim,ylim Numeric vectors of length 2 or `NULL`. Axis limits; the
#'   plot is zoomed without dropping observations outside the limits.
#' @param file Character or `NULL`. If not `NULL`, the plot is saved to this
#'   file with [ggplot2::ggsave()].
#' @param width,height Numeric. Size of the saved figure in inches.
#'
#' @return A `ggplot` object.
#'
#' @seealso [plotirf()] for impulse responses.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(t = 1:50, shock = rnorm(50))
#' d$prediction <- d$shock + rnorm(50, sd = 0.5)
#' d$proxy <- d$shock + rnorm(50, sd = 1)
#'
#' # Two series with their correlation
#' plotseries(d, x_var = "t", series = c(shock = "True shock", prediction = "Prediction"),
#'   standardize = TRUE, show_cor = TRUE
#' )
#'
#' # Three series, legend inside the panel
#' plotseries(d, legend_position = c(0.8, 0.15))
#'
#' @export
plotseries <- function(data, x_var = NULL, series = NULL,
                       colors = c("darkred", "steelblue", "darkgreen", "black", "orange", "purple"),
                       linetypes = c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash"),
                       x_lab = NULL, y_lab = NULL,
                       standardize = FALSE, show_cor = FALSE, show_n = TRUE, digits = 2,
                       legend_position = "bottom", linewidth = 0.5,
                       xlim = NULL, ylim = NULL,
                       file = NULL, width = 7, height = 4) {
  # --- Input: data frame, or time series object(s) converted with tsbox ------
  if (!is.data.frame(data)) {
    if (!requireNamespace("tsbox", quietly = TRUE)) {
      stop("data must be a data frame; time series objects require the ",
        "'tsbox' package.", call. = FALSE)
    }
    if (is.list(data)) {
      # Named list of single series: one column per list element, joined by time
      if (is.null(names(data)) || any(!nzchar(names(data)))) {
        stop("A list of time series must be named.", call. = FALSE)
      }
      tbls <- lapply(names(data), function(nm) {
        tb <- as.data.frame(tsbox::ts_tbl(data[[nm]]))[, c("time", "value")]
        stats::setNames(tb, c("time", nm))
      })
      data <- Reduce(function(a, b) dplyr::full_join(a, b, by = "time"), tbls)
      data <- data[order(data$time), , drop = FALSE]
    } else {
      data <- as.data.frame(tsbox::ts_wide(tsbox::ts_tbl(data)))
    }
    if (is.null(x_var)) x_var <- "time"
  }
  if (ncol(data) < 2) {
    stop("data must have one x-axis column and at least one series.", call. = FALSE)
  }
  if (is.null(x_var)) x_var <- names(data)[1]
  if (!is.character(x_var) || length(x_var) != 1 || !(x_var %in% names(data))) {
    stop("x_var must be the name of one column of data.", call. = FALSE)
  }

  # --- Series: column names (ids) and legend labels --------------------------
  if (is.null(series)) series <- setdiff(names(data), x_var)
  if (!is.character(series) || length(series) < 1 || anyNA(series)) {
    stop("series must be a character vector of column names.", call. = FALSE)
  }
  ids <- if (is.null(names(series))) series else names(series)
  labels <- unname(series)
  if (!all(ids %in% names(data)) || anyDuplicated(ids) || anyDuplicated(labels)) {
    stop("series must refer to distinct columns of data and have distinct labels.",
      call. = FALSE
    )
  }
  nS <- length(ids)
  if (!all(vapply(data[ids], is.numeric, logical(1)))) {
    stop("All series must be numeric.", call. = FALSE)
  }
  if (length(colors) < nS || length(linetypes) < nS) {
    stop("colors and linetypes need at least one entry per series.", call. = FALSE)
  }
  standardize <- .check_logical_scalar(standardize, "standardize")
  show_cor <- .check_logical_scalar(show_cor, "show_cor")
  show_n <- .check_logical_scalar(show_n, "show_n")
  if (show_cor && nS != 2) {
    stop("show_cor = TRUE requires exactly two series.", call. = FALSE)
  }

  # --- Correlation on the common sample (scale invariant) --------------------
  cor_lab <- NULL
  if (show_cor) {
    ok <- stats::complete.cases(data[ids])
    rho <- stats::cor(data[[ids[1]]][ok], data[[ids[2]]][ok])
    cor_lab <- paste0(
      "Correlation: ", formatC(rho, format = "f", digits = digits),
      if (show_n) paste0("  (T = ", sum(ok), ")")
    )
  }

  # --- Long format, optionally standardised ----------------------------------
  plot_data <- dplyr::bind_rows(lapply(ids, function(id) {
    v <- data[[id]]
    if (standardize) v <- (v - mean(v, na.rm = TRUE)) / stats::sd(v, na.rm = TRUE)
    data.frame(x = data[[x_var]], series_id = id, value = v)
  }))
  plot_data$series_id <- factor(plot_data$series_id, levels = ids)

  # --- Plot -------------------------------------------------------------------
  # Identical breaks and labels in both scales merge them into one legend
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data$x, y = .data$value,
      colour = .data$series_id, linetype = .data$series_id
    )
  ) +
    ggplot2::geom_line(linewidth = linewidth, na.rm = TRUE) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(colors[seq_len(nS)], ids), breaks = ids, labels = labels
    ) +
    ggplot2::scale_linetype_manual(
      values = stats::setNames(linetypes[seq_len(nS)], ids), breaks = ids, labels = labels
    ) +
    ggplot2::labs(x = x_lab, y = y_lab, colour = NULL, linetype = NULL) +
    .theme_series()

  if (!is.null(xlim) || !is.null(ylim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim, ylim = ylim)
  }
  if (!is.null(cor_lab)) {
    p <- p + ggplot2::annotate(
      "text",
      x = min(plot_data$x, na.rm = TRUE), y = min(plot_data$value, na.rm = TRUE),
      label = cor_lab, size = 3, hjust = 0, vjust = 0
    )
  }
  if (is.numeric(legend_position)) {
    p <- p + ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = legend_position
    )
  } else {
    p <- p + ggplot2::theme(legend.position = legend_position)
  }

  if (!is.null(file)) ggplot2::ggsave(file, p, width = width, height = height)
  p
}

# Theme of plotseries(): minimal theme with a black panel border and dotted
# major grid lines, as in the panels of plotirf()
.theme_series <- function(base_size = 9) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      legend.title = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.margin = ggplot2::margin(0, 0, 0, 0),
      legend.key.width = grid::unit(1.2, "cm"),
      panel.border = ggplot2::element_rect(
        colour = "black", fill = NA, linewidth = 0.2, linetype = "solid"
      ),
      panel.grid = ggplot2::element_line(
        colour = "gray", linewidth = 0.2, linetype = "dotted"
      ),
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(size = base_size + 1)
    )
}
