#' Estimate impulse responses via heteroskedasticity-based IV local projections
#'
#' Estimates impulse response functions (IRFs) using recursive
#' heteroskedasticity-IV identification (Rigobon, 2003; Rigobon and Sack, 2004;
#' Lewis, 2022; Burri and Kaufmann, 2026a, 2026b) combined with local projections
#' (Jorda, 2005). Identification exploits the difference in variance between
#' policy event days and control days to construct instruments for the
#' endogenous variables.
#'
#' @details For `E > 1`, identification is recursive and order-dependent: the
#'   column order of `y` defines both the shock ordering and the normalization
#'   variable for each shock dimension.
#'
#' @param y Numeric matrix of stationary outcome variables (T x N). The effect on
#'   the first variable in each dimension is normalized to unity at horizon 0.
#'   These variables are also used to construct heteroskedasticity-based instruments
#'   and for recursive ordering to identify multiple dimensions.
#' @param O Numeric matrix of information set variables (T x M). May be
#'   identical to `y`. Included as lags 1 through `P`.
#' @param X Numeric matrix of deterministic variables (T x K). For example,
#'   time trend, seasonal dummies or other deterministic controls. Included as is (no lags).
#'   A constant is included by default.
#' @param Ind Integer vector of length T, event indicator:
#'   \itemize{
#'     \item `0` Control day (no event)
#'     \item `1` Policy day (event)
#'     \item `2` Contaminated control day (excluded from estimation)
#'   }
#' @param P Integer. Maximum lag order for the information set. Set to `0` for
#'   no lags (regression on deterministic terms only).
#' @param H Integer. Maximum horizon (in periods) up to which IRFs are estimated.
#' @param E Integer. Number of shock dimensions to identify via recursive
#'   ordering.
#' @param norm Numeric scalar. Normalize the impact response of the first
#'   variable to a specific value. Set to `1` for standard unit-effect
#'   normalization.
#' @param interact Logical. If `TRUE`, lagged information set variables are
#'   interacted with event/non-event dummies.
#' @param cum Logical vector of length N. For each variable in `y`, whether
#'   to report the cumulative impulse response instead of the level response.
#'   If only one provided, applied to all impulse responses.
#' @param Hstep Integer. Step size between horizons. The default `1` estimates
#'   all horizons 0 through H - 1. Values greater than 1 estimate only the
#'   selected horizons.
#' @param cov_type Covariance estimator for local-projection standard errors:
#'   `"HC0"` (default) for heteroskedasticity-robust standard errors or `"NW"`
#'   for Newey-West HAC standard errors.
#' @param details Logical. If `TRUE`, code saves detailed IV results, which is slightly slower.
#'   if set to `FALSE`, returns only impulse response and standard error (e.g. for bootstrap)
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{`irf`}{Array (H x N x E) of estimated impulse responses.}
#'     \item{`se`}{Array (H x N x E) of local-projection standard errors.}
#'     \item{`IVRes`}{List of `ivreg` model objects, one per horizon, variable,
#'       and shock dimension.}
#'     \item{`Obs`}{Data frame with observation counts: `Tp` (policy days),
#'       `Tc` (control days), `To` (contaminated days), `Tt` (total used).}
#'     \item{`Method`}{Character string `"Heteroscedasticity-IV"`.}
#'     \item{`et`}{Data frame of OLS residuals on event days (used for
#'       covariance estimation and shock extraction).}
#'     \item{`Sig`}{Covariance matrix of residuals on event days (used for
#'       shock extraction).}
#'     \item{`SigR`}{Covariance matrix of residuals on control days, or `NA`
#'       if unavailable (used for shock extraction).}
#'     \item{`Psi`}{Impact matrix (N x E), equal to `irf[1, , ]`. By the
#'       package's indexing convention `HSeries` starts at 1, so the first LP
#'       uses `lead(y, 0)` (the contemporaneous value) and is labelled horizon
#'       0; `irf[1, , ]` is therefore always the impact response.}
#'     \item{`WeakData`}{Data frame of endogenous variables and instruments for
#'       the Lewis-Mertens (2025) weak instrument test.}
#'   }
#'
#' @references
#' Burri, M. and D. Kaufmann (2026a). Measuring monetary policy shocks.
#' *Economics Letters*. \doi{10.1016/j.econlet.2026.113091}
#'
#' Burri, M. and D. Kaufmann (2026b). Multiple monetary policy shocks from
#' daily data: A heteroskedasticity IV approach. IRENE Working Papers 26-06,
#' IRENE Institute of Economic Research, University of Neuchatel.
#'
#' Jorda, O. (2005). Estimation and inference of impulse responses by local
#' projections. *American Economic Review*, 95(1), 161-182.
#'
#' Lewis, D. J. (2022). Robust inference in models identified via
#' heteroskedasticity. *Review of Economics and Statistics*, 104(3), 510-524.
#'
#' Lewis, D. J. and Mertens, K. (2025). A robust test for weak instruments for 2SLS with multiple
#' endogenous regressors. *The Review of Economic Studies*, DOI: 10.1093/restud/rdaf103
#'
#' Rigobon, R. (2003). Identification through heteroskedasticity.
#' *Review of Economics and Statistics*, 85(4), 777-792.
#'
#' Rigobon, R. and Sack, B. (2004). The impact of monetary policy on asset
#' prices. *Journal of Monetary Economics*, 51(8), 1553-1575.
#'
#' @importFrom dplyr lag lead
#' @importFrom ivreg ivreg
#' @importFrom sandwich NeweyWest vcovHC
#'
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(80), ncol = 2)
#' Ind <- rep(0L, nrow(y))
#' Ind[seq(5, nrow(y), by = 5)] <- 1L
#' res <- hetiv(y = y, O = y, Ind = Ind, P = 1, H = 3, details = TRUE)
#' dim(res$irf)
#'
#' @export
hetiv <- function(y, O, X = NULL, Ind, P, H, E = 1, norm = 1, interact = FALSE,
                  cum = FALSE, Hstep = 1, cov_type = "HC0", details = FALSE) {
  args <- .validate_estimator_inputs(
    y = y, O = O, X = X, Ind = Ind, P = P, H = H, E = E, norm = norm,
    cum = cum, Hstep = Hstep, cov_type = cov_type
  )
  y <- args$y
  O <- args$O
  X <- args$X
  Ind <- args$Ind
  P <- args$P
  H <- args$H
  E <- args$E
  norm <- args$norm
  cum <- args$cum
  Hstep <- args$Hstep
  cov_type <- args$cov_type
  interact <- .check_logical_scalar(interact, "interact")
  details <- .check_logical_scalar(details, "details")

  # Collect various properties of the data and observations to be used
  Nobs <- dim(y)[1] # Number of observations in y
  beg <- P + 1 # Start of the sample
  end <- Nobs - H + 1 # End of the sample
  N <- dim(y)[2] # Number of variables in y
  M <- dim(O)[2] # Number of variables in O
  K <- if (!is.null(X)) dim(X)[2] else 0

  if (sum(is.na(y)) > 0) {
    warning("Missing values in y")
  }
  if (sum(is.na(O)) > 0) {
    warning("Missing values in O")
  }
  if (sum(is.na(Ind)) > 0) {
    warning("Missing values in Ind")
  }
  if (!is.null(X) && sum(is.na(X)) > 0) warning("Missing values in X")

  # Specify at which horizons IRF should be computed
  HSeries <- seq(1, H, Hstep)
  HNum <- length(HSeries)

  # Set up data set and and objects to save results
  if (is.null(X)) {
    DataM <- data.frame(y, O, Ind)
    colnames(DataM) <- c(paste0("y", 1:N), paste0("o", 1:M), "Ind")
  } else {
    DataM <- data.frame(y, O, X, Ind)
    colnames(DataM) <- c(paste0("y", 1:N), paste0("o", 1:M), paste0("x", 1:K), "Ind")
  }
  irfest <- array(NA, dim = c(HNum, N, E))
  irfse <- array(NA, dim = c(HNum, N, E))
  IVRes <- list()
  OLSRes <- list()
  ORTHRes <- list()

  # Compute lagged variables included in information set O(t-1)...O(t-P)
  infoVars <- c()
  if (P > 0) {
    for (j in 1:M) {
      for (p in 1:P) {
        if (interact == TRUE) {
          # Compute lags of variables interacted with event dummies (excluding contaminated events)
          for (e in 0:1) {
            DataM[, paste0("o", j, ".l", p, ".i", e)] <-
              dplyr::lag(DataM[, paste0("o", j)], p) * (Ind == e)
            infoVars <- c(infoVars, paste0("o", j, ".l", p, ".i", e))

            # Include interacted constant (only for e == 1 to avoid perfect multicollinearity)
            if (p == 1 && e == 1) {
              DataM[, paste0("i", e)] <- (Ind == e)
              infoVars <- c(infoVars, paste0("i", e))
            }
          }
        } else {
          # Compute lags of variables
          DataM[, paste0("o", j, ".l", p)] <- dplyr::lag(DataM[, paste0("o", j)], p)
          infoVars <- c(infoVars, paste0("o", j, ".l", p))
        }
      }
    }
  }

  # Otherwise, only regress on a constant
  if (P == 0) {
    if (interact == TRUE) {
      # Include interacted constant (only for e == 1 to avoid perfect multicollinearity)
      DataM[, paste0("i", 1)] <- (Ind == 1)
      infoVars <- c(infoVars, paste0("i", 1))
    } else {
      infoVars <- "1"
    }
  }

  # Add deterministic variables
  if (!is.null(X)) {
    for (k in 1:K) {
      infoVars <- c(infoVars, paste0("x", k))
    }
  }

  # Compute heteroskedasticity-based instrument separately for every dimension
  for (e in 1:E) {
    DataM$Ind <- Ind

    # Set up the shock variable (instrumented variable), which is the e th variable in y
    DataM$shockVar <- DataM[, paste0("y", e)]

    # Set up the instrument variable, which is also based on the e th variable in y
    DataM$Ze <- DataM[, paste0("y", e)]

    # Set up the control variables for recursive ordering in case of E>1
    recVars <- c()
    recInst <- c()
    if (e > 1) {
      for (q in 1:(e - 1)) {
        recVars <- c(recVars, paste0("y", q))
        recInst <- c(recInst, paste0("Z", q))
      }
    }

    # Set up control variables for the local projections
    controls.info <- unique(c(infoVars))
    controls.info <- controls.info[controls.info != ""]

    controls.lp <- unique(c(infoVars, recVars))
    controls.lp <- controls.lp[controls.lp != ""]

    controls.iv <- unique(c(infoVars, recInst))
    controls.iv <- controls.iv[controls.iv != ""]

    # Account for missing data in instrumental variable by setting it to an other event
    # which is not used in the estimation
    DataM$Ind[is.na(DataM$Ze)] <- 2

    # Set up event days (Policy event, control, and other event)
    DataM$Event <- (DataM$Ind == 1)
    DataM$NoEvent <- (DataM$Ind == 0)
    DataM$OthEvent <- (DataM$Ind == 2)

    # Shorten data to subset of observations without missing values at beginning or end
    DataMSub <- DataM[beg:end, ]

    # Orthogonalize the variable used to construct the instrument if we include
    # lags of dependent variable. See Lewis (2022). Note that if no control
    # variables are used, this is just a regression of Ze on a constant.
    myFormula <- paste0("Ze ~ ", paste(controls.info, collapse = "+"))
    orthModel <- lm(
      as.formula(myFormula),
      data = subset(DataMSub, Ind < 2),
      na.action = "na.exclude"
    )
    DataMSub$Ze <- residuals(orthModel, na.action = "na.exclude")

    # Compute instrument and F-Statistic (see Lewis, 2022, and Rigobon and Sachs, 2004)
    Te <- sum(DataMSub$Event)
    Tn <- sum(DataMSub$NoEvent)
    To <- sum(DataMSub$OthEvent)
    Tt <- Te + Tn
    if (Te == 0 || Tn == 0) {
      stop("At least one event day (Ind == 1) and one control day (Ind == 0) are required.")
    }
    DataMSub$Z <- (DataMSub$Event * Tt / Te - DataMSub$NoEvent * Tt / Tn) * DataMSub$Ze

    # Save instrument in original data for later use
    DataM[beg:end, paste0("Z", e)] <- DataMSub$Z
    DataM[beg:end, "Z"] <- DataMSub$Z

    # Estimate the impulse responses based on heteroskedasticity-based IV for every variable in y
    for (i in 1:N) {
      # Set up dependent variable
      DataM$depVar <- DataM[, paste0("y", i)]
      cumi <- cum[i]

      for (h_idx in seq_along(HSeries)) {
        h <- HSeries[h_idx]

        # Compute dependent variable at horizon h (level or cumulative response)
        if (cumi == TRUE) {
          for (f in 1:h) {
            if (f == 1) {
              DataM$depVar.h <- dplyr::lead(DataM$depVar, f - 1)
            } else {
              DataM$depVar.h <- DataM$depVar.h + dplyr::lead(DataM$depVar, f - 1)
            }
          }
        } else {
          DataM$depVar.h <- dplyr::lead(DataM$depVar, h - 1)
        }

        # Shorten data to subset which contains no missing values
        DataMSub <- DataM[beg:end, ]

        # LP (Jorda, 2005): instrument shockVar with Z, control for information set
        myFormula <- paste0(
          "depVar.h ~ shockVar + ",
          paste(controls.lp, collapse = "+"),
          paste0(
            "| Z + ",
            paste(controls.iv, collapse = "+")
          )
        )

        # Estimate IV regression and compute local-projection standard errors
        IV.mod <- ivreg::ivreg(as.formula(myFormula), data = subset(DataMSub, Ind < 2))
        if (cov_type == "NW") {
          IV.vcov <- sandwich::NeweyWest(IV.mod, prewhite = FALSE, adjust = TRUE)
        } else {
          IV.vcov <- sandwich::vcovHC(IV.mod, type = "HC0")
        }
        IV.se <- sqrt(diag(IV.vcov))


        if (details == TRUE) {
          # Compute OLS residuals on event and control days (once per shock, at h=1)
          # Save residuals for later computation of variance-covariance matrix
          if (e == 1 && h == 1) {
            DataMSub$depVar.h2 <- DataMSub$depVar.h
            DataMSub$depVar.h2[DataMSub$Ind == 2] <- NA

            myFormula <- paste0("depVar.h2 ~", paste(controls.info, collapse = "+"))
            OLS.mod <- lm(
              as.formula(myFormula),
              data = subset(DataMSub, Ind < 2),
              na.action = "na.exclude"
            )

            eti <- residuals(OLS.mod)
            eti[DataMSub$Event != 1] <- NA

            vti <- residuals(OLS.mod)
            vti[DataMSub$NoEvent != 1] <- NA

            # Account for missing values when saving back to original data set
            DataM$eti <- NA
            DataM$eti[beg:end] <- eti
            eti <- DataM$eti

            DataM$vti <- NA
            DataM$vti[beg:end] <- vti
            vti <- DataM$vti

            if (i == 1) {
              et <- data.frame(eti)
              vt <- data.frame(vti)
            } else {
              et <- data.frame(et, eti)
              vt <- data.frame(vt, vti)
            }
          }
        }

        # Normalizes IRFs to a specific value. Because initial response
        irfest[h_idx, i, e] <- IV.mod$coefficients["shockVar"] * norm
        irfse[h_idx, i, e] <- IV.se["shockVar"] * abs(norm)

        # Save IV results for every variable and every horizon
        if (details == TRUE) {
          IVRes[[paste0("IV.h", h, ".n", i, ".e", e)]] <- IV.mod

          if (h == 1) {
            # Save OLS and orthogonalization results for horizon 1 only to save memory
            OLSRes[[paste0("OLS.n", i)]] <- OLS.mod

            if (i == e) {
              ORTHRes[[paste0("ORTH.h", h, ".n", i, ".e", e)]] <- orthModel
            }
          }
        }
      }
    }
  }

  # Label rows of impulse responses to start at 0 (immediate response)
  dimnames(irfest)[[1]] <- HSeries - 1
  dimnames(irfse)[[1]] <- HSeries - 1

  if (details == TRUE) {
    # Compute variance-covariance matrix of residuals on event days, and impact matrix
    Sig <- var(et, use = "complete.obs")

    if (sum(!is.na(vt) > 0)) {
      SigR <- var(vt, use = "complete.obs")
    } else {
      SigR <- NA
    }
    Psi <- matrix(irfest[1, , , drop = FALSE], nrow = N, ncol = E)

    # Save data for weak instruments test by Lewis-Mertens (2025)
    if (controls.info[1] != "1") {
      WeakData <- data.frame(
        DataM[, paste0("y", 1:N)],
        DataM[, paste0("Z", 1:E)],
        DataM[, controls.info]
      )
    } else {
      WeakData <- data.frame(DataM[, paste0("y", 1:N)], DataM[, paste0("Z", 1:E)])
    }
    wk_names <- c(paste0("y", seq_len(N)), paste0("Z", seq_len(E)))
    if (controls.info[1] != "1") wk_names <- c(wk_names, controls.info)
    colnames(WeakData) <- wk_names

    # Set data to missing if indicator is equal to 2
    WeakData[DataM$Ind == 2, ] <- NA

    Obs <- data.frame(Tp = Te, Tc = Tn, To = To, Tt = Tt)
  }

  Method <- "Heteroscedasticity-IV"

  if (details == TRUE) {
    return(list(
      irf = irfest, se = irfse,
      IVRes = IVRes, OLSRes = OLSRes, ORTHRes = ORTHRes,
      Obs = Obs, Method = Method,
      et = as.matrix(et), Sig = Sig, SigR = SigR, Psi = Psi, WeakData = WeakData
    ))
  } else {
    return(list(irf = irfest, se = irfse, Method = Method))
  }
}
