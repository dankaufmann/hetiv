.as_numeric_matrix <- function(x, arg, nrow = NULL, ncol = NULL,
                               square = FALSE, allow_na = TRUE) {
  if (missing(x) || is.null(x)) {
    stop(arg, " must not be NULL.", call. = FALSE)
  }

  x <- as.matrix(x)
  if (!is.numeric(x)) {
    stop(arg, " must be numeric.", call. = FALSE)
  }
  if (!allow_na && anyNA(x)) {
    stop(arg, " must not contain missing values.", call. = FALSE)
  }
  if (any(!is.finite(x[!is.na(x)]))) {
    stop(arg, " must contain only finite values or NA.", call. = FALSE)
  }
  if (!is.null(nrow) && nrow(x) != nrow) {
    stop(arg, " must have ", nrow, " rows.", call. = FALSE)
  }
  if (!is.null(ncol) && ncol(x) != ncol) {
    stop(arg, " must have ", ncol, " columns.", call. = FALSE)
  }
  if (square && nrow(x) != ncol(x)) {
    stop(arg, " must be a square matrix.", call. = FALSE)
  }

  x
}

.as_numeric_array3 <- function(x, arg, dim = NULL, allow_na = TRUE) {
  if (!is.array(x) || length(dim(x)) != 3) {
    stop(arg, " must be a three-dimensional numeric array.", call. = FALSE)
  }
  if (!is.numeric(x)) {
    stop(arg, " must be numeric.", call. = FALSE)
  }
  if (!allow_na && anyNA(x)) {
    stop(arg, " must not contain missing values.", call. = FALSE)
  }
  if (any(!is.finite(x[!is.na(x)]))) {
    stop(arg, " must contain only finite values or NA.", call. = FALSE)
  }
  if (!is.null(dim) && !identical(base::dim(x), as.integer(dim))) {
    stop(arg, " must have dimensions ", paste(dim, collapse = " x "), ".",
      call. = FALSE
    )
  }

  x
}

.check_integerish_scalar <- function(x, arg, min = -Inf, max = Inf) {
  if (
    !is.numeric(x) || length(x) != 1 || is.na(x) || !is.finite(x) ||
      x != as.integer(x) || x < min || x > max
  ) {
    range_msg <- ""
    if (is.finite(min) && is.finite(max)) {
      range_msg <- paste0(" between ", min, " and ", max)
    } else if (is.finite(min)) {
      range_msg <- paste0(" at least ", min)
    } else if (is.finite(max)) {
      range_msg <- paste0(" at most ", max)
    }
    stop(arg, " must be an integer", range_msg, ".", call. = FALSE)
  }

  as.integer(x)
}

.check_numeric_scalar <- function(x, arg, min = -Inf, max = Inf,
                                  allow_na = FALSE) {
  bad <- !is.numeric(x) || length(x) != 1 ||
    (!allow_na && is.na(x)) || (!is.na(x) && !is.finite(x)) ||
    (!is.na(x) && (x < min || x > max))
  if (bad) {
    range_msg <- ""
    if (is.finite(min) && is.finite(max)) {
      range_msg <- paste0(" between ", min, " and ", max)
    } else if (is.finite(min)) {
      range_msg <- paste0(" at least ", min)
    } else if (is.finite(max)) {
      range_msg <- paste0(" at most ", max)
    }
    stop(arg, " must be a numeric scalar", range_msg, ".", call. = FALSE)
  }

  x
}

.check_logical_scalar <- function(x, arg) {
  if (!is.logical(x) || length(x) != 1 || is.na(x)) {
    stop(arg, " must be TRUE or FALSE.", call. = FALSE)
  }

  x
}

.check_choice <- function(x, arg, choices) {
  match <- pmatch(x, choices)
  if (!is.character(x) || length(x) != 1 || is.na(x) ||
    is.na(match) || match == 0) {
    stop(arg, " must be one of ",
      paste0("\"", choices, "\"", collapse = ", "), ".",
      call. = FALSE
    )
  }

  choices[[match]]
}

.check_indicator <- function(Ind, nobs) {
  if (length(Ind) != nobs) {
    stop("Ind must have the same length as nrow(y).", call. = FALSE)
  }
  if (!all(Ind %in% c(0, 1, 2, NA))) {
    stop("Ind must contain only 0, 1, 2, or NA.", call. = FALSE)
  }

  as.integer(Ind)
}

.check_cum <- function(cum, n) {
  if (!is.logical(cum) || anyNA(cum) || !(length(cum) %in% c(1, n))) {
    stop("cum must be a non-missing logical scalar or a logical vector of length ",
      n, ".",
      call. = FALSE
    )
  }
  if (length(cum) == 1) {
    cum <- rep(cum, n)
  }

  cum
}

.validate_estimator_inputs <- function(y, O, X, Ind, P, H, E, norm, cum,
                                       Hstep, cov_type) {
  y <- .as_numeric_matrix(y, "y")
  O <- .as_numeric_matrix(O, "O", nrow = nrow(y))
  if (!is.null(X)) {
    X <- .as_numeric_matrix(X, "X", nrow = nrow(y))
  }

  P <- .check_integerish_scalar(P, "P", min = 0)
  H <- .check_integerish_scalar(H, "H", min = 1)
  E <- .check_integerish_scalar(E, "E", min = 1)
  Hstep <- .check_integerish_scalar(Hstep, "Hstep", min = 1)
  norm <- .check_numeric_scalar(norm, "norm")
  Ind <- .check_indicator(Ind, nrow(y))
  cum <- .check_cum(cum, ncol(y))
  cov_type <- .check_choice(cov_type, "cov_type", c("HC0", "NW"))

  if (E > ncol(y)) {
    stop("E cannot exceed the number of columns in y.", call. = FALSE)
  }
  if (nrow(y) <= P + H) {
    stop("Not enough observations: nrow(y) must exceed P + H.", call. = FALSE)
  }
  if (sum(Ind == 1, na.rm = TRUE) == 0 || sum(Ind == 0, na.rm = TRUE) == 0) {
    stop("Ind must contain at least one event day (1) and one control day (0).",
      call. = FALSE
    )
  }

  list(
    y = y, O = O, X = X, Ind = Ind, P = P, H = H, E = E, norm = norm,
    cum = cum, Hstep = Hstep, cov_type = cov_type
  )
}

.check_horizon_labels <- function(x, arg) {
  labels <- dimnames(x)[[1]]
  if (is.null(labels) || anyNA(suppressWarnings(as.numeric(labels)))) {
    stop(arg, " must have numeric horizon labels in dimnames()[[1]].",
      call. = FALSE
    )
  }

  as.numeric(labels)
}

.validate_irf_array <- function(x, arg) {
  if (!is.array(x) || !(length(dim(x)) %in% c(2, 3))) {
    stop(arg, " must be a numeric matrix or three-dimensional array.",
      call. = FALSE
    )
  }
  if (!is.numeric(x)) {
    stop(arg, " must be numeric.", call. = FALSE)
  }
  if (any(!is.finite(x[!is.na(x)]))) {
    stop(arg, " must contain only finite values or NA.", call. = FALSE)
  }

  x
}

.validate_plot_common <- function(HSeries, HTick, Labels, N) {
  .check_integerish_scalar(HTick, "HTick", min = 1)
  if (!is.character(Labels) || length(Labels) != N || anyNA(Labels)) {
    stop("Labels must be a character vector of length ", N, ".", call. = FALSE)
  }
  if (length(HSeries) == 0 || anyNA(HSeries)) {
    stop("Horizon labels must be non-missing numeric values.", call. = FALSE)
  }

  invisible(TRUE)
}
