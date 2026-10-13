# Internal helpers shared by plotirf() and plot2irf(): normalise the different
# ways an IRF can be supplied and compute confidence bands.

# Normalise one IRF input into a list(irf, se, boot, name).
#
# `x` may be
#   (a) a numeric array (H x N x E or H x N) of point estimates, in which case
#       `se` (same shape) gives optional standard errors, or
#   (b) an estimator result, i.e. a list with `$irf` and optionally `$se`,
#       `$boot` (H x N x E x B bootstrap draws) and `$Method`.
# Two-dimensional inputs (single shock) are reshaped to H x N x 1 together
# with `se` and `boot`, so that all downstream code can assume 3-D arrays.
# An all-zero `se` is treated as "no standard errors" (line only).
.as_irf_spec <- function(x, arg, se = NULL, name = NULL) {
  boot <- NULL
  if (is.list(x) && !is.data.frame(x)) {
    if (is.null(x$irf)) {
      stop(arg, " must be a numeric array or a list with an element 'irf'.",
        call. = FALSE
      )
    }
    if (!is.null(se)) {
      stop("IRFse cannot be combined with a list input; supply 'se' inside ",
        "the list instead.", call. = FALSE)
    }
    se <- x$se
    boot <- x$boot
    if (is.null(name) && !is.null(x$Method)) name <- x$Method
    x <- x$irf
  }
  irf <- .validate_irf_array(x, arg)
  HSeries <- .check_horizon_labels(irf, arg)
  noDims <- length(dim(irf))
  HNum <- dim(irf)[1]
  N <- dim(irf)[2]
  E <- if (noDims == 3) dim(irf)[3] else 1L
  dim3 <- c(HNum, N, E)

  if (!is.null(se)) {
    se <- .validate_irf_array(se, paste0(arg, " standard errors"))
    if (!identical(as.integer(dim(se)), as.integer(dim(irf)))) {
      stop("Standard errors of ", arg, " must have the same dimensions as the ",
        "point estimates.", call. = FALSE)
    }
    dim(se) <- dim3
    if (all(se == 0, na.rm = TRUE)) se <- NULL
  }
  if (!is.null(boot)) {
    if (!is.array(boot) || !is.numeric(boot) ||
      !(length(dim(boot)) %in% c(noDims + 1L))) {
      stop("Bootstrap draws of ", arg, " must be a numeric array with one ",
        "more dimension than the point estimates.", call. = FALSE)
    }
    if (!identical(as.integer(dim(boot)[seq_len(noDims)]), as.integer(dim(irf)))) {
      stop("Bootstrap draws of ", arg, " must match the dimensions of the ",
        "point estimates in all but the last dimension.", call. = FALSE)
    }
    dim(boot) <- c(dim3, dim(boot)[noDims + 1L])
  }
  dim(irf) <- dim3
  dimnames(irf) <- list(as.character(HSeries), NULL, NULL)

  list(irf = irf, se = se, boot = boot, name = name, HSeries = HSeries)
}

# Confidence bands for one IRF spec at the levels in `ci`.
#
# Returns NULL if the spec has neither bootstrap draws nor standard errors.
# Otherwise returns list(lower, upper), arrays of dimension
# H x N x E x length(ci).
#
# Bootstrap draws take precedence over standard errors. With
# boot_type = "recentred" (default) the bands are recentred on the point
# estimate: quantiles are computed from the deviations of the draws from their
# bootstrap mean and added to the point estimate, so that the band is
# asymmetric but centred on the reported IRF. With boot_type = "percentile"
# the raw quantiles of the draws are used.
.irf_bands <- function(spec, ci, boot_type = c("recentred", "percentile")) {
  boot_type <- match.arg(boot_type)
  irf <- spec$irf
  dims <- c(dim(irf), length(ci))
  lower <- array(NA_real_, dim = dims)
  upper <- array(NA_real_, dim = dims)

  if (!is.null(spec$boot)) {
    boot <- spec$boot
    if (boot_type == "recentred") {
      # Deviations from the bootstrap mean, added back to the point estimate
      bmean <- apply(boot, c(1, 2, 3), mean, na.rm = TRUE)
      dev <- sweep(boot, c(1, 2, 3), bmean, "-")
      centre <- irf
    } else {
      # Raw percentile interval of the draws
      dev <- boot
      centre <- 0
    }
    for (k in seq_along(ci)) {
      lower[, , , k] <- centre + apply(dev, c(1, 2, 3), stats::quantile,
        probs = (1 - ci[k]) / 2, na.rm = TRUE, names = FALSE
      )
      upper[, , , k] <- centre + apply(dev, c(1, 2, 3), stats::quantile,
        probs = (1 + ci[k]) / 2, na.rm = TRUE, names = FALSE
      )
    }
  } else if (!is.null(spec$se)) {
    for (k in seq_along(ci)) {
      z <- stats::qnorm((1 + ci[k]) / 2)
      lower[, , , k] <- irf - z * spec$se
      upper[, , , k] <- irf + z * spec$se
    }
  } else {
    return(NULL)
  }

  list(lower = lower, upper = upper)
}

# Validate confidence levels for the plotting functions: between one and two
# levels, each strictly inside (0, 1), returned sorted widest first so that
# wider bands are drawn behind narrower ones.
.check_ci_levels <- function(ci) {
  if (
    !is.numeric(ci) || length(ci) < 1 || length(ci) > 2 || anyNA(ci) ||
      any(ci <= 0 | ci >= 1)
  ) {
    stop("ci must contain one or two confidence levels between 0 and 1.",
      call. = FALSE
    )
  }
  sort(ci, decreasing = TRUE)
}
