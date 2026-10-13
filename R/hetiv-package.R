#' @keywords internal
#' @importFrom stats as.formula complete.cases integrate lm optim qchisq qnorm residuals rnorm rt sd var
#' @importFrom rlang .data
"_PACKAGE"

utils::globalVariables(c(
  ".data", "Estimate", "Horizon", "IRF1", "IRF2", "lower", "lower1", "lower2",
  "pval", "Series", "upper", "upper1", "upper2"
))
