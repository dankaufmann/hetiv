#' Optimal lag length of a VAR
#'
#' Selects the lag length of a VAR (with a constant and, optionally, exogenous
#' variables) by an information criterion, using [vars::VAR()]. The result can
#' be used as the number of lags `P` of the information set in [hetiv()],
#' [proxyiv()], [recols()] and [svariv()].
#'
#' @param Y Numeric matrix or data frame (T x N) of endogenous variables.
#' @param X Numeric matrix or data frame (T x K) of exogenous variables, or
#'   `NULL` (default).
#' @param crit Character. Information criterion: `"AIC"` (default, following
#'   Montiel Olea et al., 2025), `"HQ"`, `"SC"` or `"FPE"`.
#' @param lag.max Integer. Largest lag length considered. Default 15.
#'
#' @return Integer. The lag length selected by the criterion.
#'
#' @details Requires the `vars` package, which is not installed automatically
#'   with `hetiv`.
#'
#' @references
#' Montiel Olea, J. L., Plagborg-Moller, M., Qian, E. and Wolf, C. K. (2025).
#' Local projections or VARs? A primer for macroeconomists. NBER
#' Macroeconomics Annual, forthcoming.
#'
#' @examples
#' if (requireNamespace("vars", quietly = TRUE)) {
#'   set.seed(1)
#'   Y <- matrix(rnorm(400), ncol = 2)
#'   for (t in 2:200) Y[t, ] <- 0.5 * Y[t - 1, ] + Y[t, ]
#'   lag_opt(Y, lag.max = 4)
#' }
#'
#' @export
lag_opt <- function(Y, X = NULL, crit = "AIC", lag.max = 15) {
  crit <- .check_choice(crit, "crit", c("AIC", "HQ", "SC", "FPE"))
  lag.max <- .check_integerish_scalar(lag.max, "lag.max", min = 1)
  if (!requireNamespace("vars", quietly = TRUE)) {
    stop("lag_opt() requires the 'vars' package. Install it with ",
      "install.packages(\"vars\").", call. = FALSE)
  }

  # vars::VAR() needs column names
  Y <- as.matrix(Y)
  if (is.null(colnames(Y))) colnames(Y) <- paste0("y", seq_len(ncol(Y)))
  if (!is.null(X)) {
    X <- as.matrix(X)
    if (is.null(colnames(X))) colnames(X) <- paste0("x", seq_len(ncol(X)))
  }

  # Lag order of a VAR (X = NULL) or a VARX chosen by the criterion
  varEst <- if (is.null(X)) {
    vars::VAR(Y, type = "const", ic = crit, lag.max = lag.max)
  } else {
    vars::VAR(Y, exogen = X, type = "const", ic = crit, lag.max = lag.max)
  }

  as.integer(varEst$p)
}
