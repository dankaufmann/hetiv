#' Estimate impulse responses via proxy-IV local projections
#'
#' Estimates impulse response functions (IRFs) using user-provided external
#' instruments (proxies) combined with local projections (Jordà, 2005). The
#' proxy variables serve directly as instruments for the endogenous shock
#' variables. Optionally imposes recursive zero restrictions across shock
#' dimensions and supports deterministic controls and interaction terms,
#' following the same conventions as [hetiv()].
#'
#' @param y Numeric matrix of stationary outcome variables (T x N). The effect
#'   on the first variable in each dimension is normalized to `norm` at horizon
#'   0. These variables are also used as the endogenous regressors instrumented
#'   by the columns of `Z`.
#' @param O Numeric matrix of information set variables (T x M). May be
#'   identical to `y`. Included as lags 1 through `P`.
#' @param Z Numeric matrix of external instruments (T x E). Column `e` is used
#'   as the proxy for shock dimension `e`. Missing values on control days
#'   (`Ind == 0`) are treated as contaminated and excluded from estimation;
#'   missing values on policy days (`Ind == 1`) are retained so that those
#'   observations remain available for shock prediction even when the instrument
#'   is unobserved.
#' @param X Numeric matrix of deterministic variables (T x K), or `NULL`
#'   (default). May include a constant, time trend, or seasonal dummies.
#'   Included as-is (no lags).
#' @param Ind Integer vector of length T, event indicator:
#'   \itemize{
#'     \item `0` Control day (no event)
#'     \item `1` Policy day (event)
#'     \item `2` Contaminated control day (excluded from estimation)
#'   }
#' @param P Integer. Maximum lag order for the information set. Set to `0` for
#'   no lags (regression on constant only).
#' @param H Integer. Maximum horizon (in periods) up to which IRFs are
#'   estimated.
#' @param E Integer. Number of shock dimensions to identify. Default `1`.
#' @param norm Numeric scalar. Normalize the impact response of the first
#'   variable to this value. Set to `1` for standard unit-effect normalization.
#' @param cum Logical scalar or vector of length N. If `TRUE` for variable `i`,
#'   the cumulative IRF is reported. A single value is recycled to all
#'   variables. Default `FALSE`.
#' @param Hstep Integer. Step size between horizons. The default `1` estimates
#'   all horizons 0 through H - 1. Values greater than 1 are intended only for
#'   fast testing; they are only safe when `Hstep >= H` (a single horizon is
#'   stored). For complete IRF estimation always use `Hstep = 1`. Default `1`.
#' @param recursive Logical. If `TRUE`, imposes recursive zero restrictions
#'   across shock dimensions: for shock `e > 1`, the variables and instruments
#'   from dimensions `1, ..., e-1` are added as controls. Default `FALSE`.
#' @param details Logical. If `TRUE`, returns detailed results including IV
#'   model objects, OLS residuals, and covariance matrices. If `FALSE` (default),
#'   returns only impulse responses and standard errors (faster; use for
#'   bootstrap).
#'
#' @return A named list. Always contains:
#'   \describe{
#'     \item{`irf`}{Array (H x N x E) of estimated impulse responses.}
#'     \item{`se`}{Array (H x N x E) of HC0 heteroscedasticity-robust standard
#'       errors.}
#'     \item{`Method`}{Character string `"Proxy-IV"`.}
#'   }
#'   With `details = TRUE`, additionally contains:
#'   \describe{
#'     \item{`IVRes`}{List of `ivreg` model objects, one per horizon, variable,
#'       and shock dimension.}
#'     \item{`OLSRes`}{List of OLS model objects used for residual-based
#'       covariance estimation, one per outcome variable.}
#'     \item{`Obs`}{Data frame with observation counts: `Tp` (policy days),
#'       `Tc` (control days), `To` (contaminated days), `Tt` (total used).}
#'     \item{`et`}{Data frame of OLS residuals on event days.}
#'     \item{`Sig`}{Covariance matrix of residuals on event days.}
#'     \item{`SigR`}{Covariance matrix of residuals on control days, or `NA`
#'       if unavailable.}
#'     \item{`Psi`}{Impact matrix (N x E), equal to `irf[1, , ]`. By the
#'       package's indexing convention `HSeries` starts at 1, so the first LP
#'       uses `lead(y, 0)` (the contemporaneous value) and is labelled horizon
#'       0; `irf[1, , ]` is therefore always the impact response.}
#'     \item{`WeakData`}{Data frame of endogenous variables and instruments for
#'       the Lewis-Mertens (2025) weak instrument test.}
#'   }
#'
#' @references
#' Jordà, Ò. (2005). Estimation and Inference of Impulse Responses by Local
#' Projections. *American Economic Review*, 95(1), 161–182.
#'
#' Lewis, D. J. and Mertens, K. (2025). A Robust Test for Weak Instruments for
#' 2SLS with Multiple Endogenous Regressors. *The Review of Economic Studies*,
#' DOI: 10.1093/restud/rdaf103
#'
#' Mertens, K. and Ravn, M. O. (2013). The Dynamic Effects of Personal and
#' Corporate Income Tax Changes in the United States. *American Economic
#' Review*, 103(4), 1212–1247.
#'
#' Stock, J. H. and Watson, M. W. (2018). Identification and Estimation of
#' Dynamic Causal Effects in Macroeconomics Using External Instruments.
#' *Economic Journal*, 128(610), 917–948.
#'
#' @importFrom dplyr lag lead
#' @importFrom ivreg ivreg
#' @importFrom sandwich vcovHC
#'
#' @export
proxyiv <- function(y, O, Z, X = NULL, Ind, P, H, E = 1, norm = 1,
                    cum = FALSE, Hstep = 1,
                    recursive = FALSE, details = FALSE) {

  # Collect properties of the data
  Nobs <- dim(y)[1]
  beg  <- P + 1
  end  <- Nobs - H + 1
  N    <- dim(y)[2]
  M    <- dim(O)[2]
  K    <- if (!is.null(X)) dim(X)[2] else 0

  if (sum(is.na(y)) > 0)   warning("Missing values in y")
  if (sum(is.na(O)) > 0)   warning("Missing values in O")
  if (sum(is.na(Ind)) > 0) warning("Missing values in Ind")
  if (!is.null(X) && sum(is.na(X)) > 0) warning("Missing values in X")

  # Recycle cum to all N variables if a single value is given
  if (length(cum) == 1) cum <- rep(cum, N)

  # Horizons at which IRFs are estimated
  HSeries <- seq(1, H, Hstep)
  HNum    <- length(HSeries)

  # Build main data frame
  if (is.null(X)) {
    DataM <- data.frame(y, O, Ind, Z)
    colnames(DataM) <- c(paste0("y", 1:N), paste0("o", 1:M), "Ind", paste0("z", 1:E))
  } else {
    DataM <- data.frame(y, O, X, Ind, Z)
    colnames(DataM) <- c(paste0("y", 1:N), paste0("o", 1:M), paste0("x", 1:K), "Ind", paste0("z", 1:E))
  }

  irfest <- array(NA, dim = c(HNum, N, E))
  irfse  <- array(NA, dim = c(HNum, N, E))
  IVRes  <- list()
  OLSRes <- list()

  # Build information set: lags of O
  infoVars <- c()
  if (P > 0) {
    for (j in 1:M) {
      for (p in 1:P) {
        DataM[, paste0("o", j, ".l", p)] <- dplyr::lag(DataM[, paste0("o", j)], p)
        infoVars <- c(infoVars, paste0("o", j, ".l", p))
      }
    }
  }
  if (P == 0) {
    infoVars <- "1"
  }

  # Add deterministic variables (X) without lags
  if (!is.null(X)) {
    for (k in 1:K) {
      infoVars <- c(infoVars, paste0("x", k))
    }
  }

  # Estimate proxy-IV impulse responses separately for every shock dimension
  for (e in 1:E) {

    # Endogenous variable: the e-th variable in y
    DataM$shockVar <- DataM[, paste0("y", e)]

    # Instrument: the e-th column of user-provided Z
    DataM$Ze <- DataM[, paste0("z", e)]

    # Recursive ordering: for shock e > 1, add variables and instruments from
    # previously identified dimensions as additional controls
    recVars <- c()
    recInst <- c()
    if (recursive == TRUE && e > 1) {
      for (q in 1:(e - 1)) {
        recVars <- c(recVars, paste0("y", q))
        recInst <- c(recInst, paste0("Z", q))
      }
    }

    controls.info <- unique(infoVars[infoVars != ""])
    controls.lp   <- unique(c(infoVars, recVars))
    controls.lp   <- controls.lp[controls.lp != ""]
    controls.iv   <- unique(c(infoVars, recInst))
    controls.iv   <- controls.iv[controls.iv != ""]

    DataM$Event    <- (DataM$Ind == 1)
    DataM$NoEvent  <- (DataM$Ind == 0)
    DataM$OthEvent <- (DataM$Ind == 2)

    DataMSub <- DataM[beg:end, ]

    Te <- sum(DataMSub$Event)
    Tn <- sum(DataMSub$NoEvent)
    To <- sum(DataMSub$OthEvent)
    Tt <- Te + Tn

    # Use the user-provided proxy directly as instrument (no orthogonalization;
    # 2SLS projection in ivreg handles the required partialling out)
    DataMSub$Z <- DataMSub$Ze
    DataM[beg:end, paste0("Z", e)] <- DataMSub$Z
    DataM[beg:end, "Z"]            <- DataMSub$Z

    # Estimate impulse responses for every outcome variable in y
    for (i in 1:N) {

      DataM$depVar <- DataM[, paste0("y", i)]
      cumi <- cum[i]

      for (h in HSeries) {

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

        DataMSub <- DataM[beg:end, ]

        # Proxy-IV LP (Jordà, 2005): instrument shockVar with Z, control for information set
        myFormula <- paste0("depVar.h ~ shockVar + ",
                            paste(controls.lp, collapse = "+"),
                            " | Z + ",
                            paste(controls.iv, collapse = "+"))

        IV.mod <- suppressWarnings(ivreg::ivreg(as.formula(myFormula), data = subset(DataMSub, Ind < 2)))
        IV.se  <- sqrt(diag(sandwich::vcovHC(IV.mod, type = "HC0")))

        # Compute OLS residuals for covariance estimation (once per outcome variable, at e=1, h=1)
        if (details == TRUE && e == 1 && h == 1) {
          DataMSub$depVar.h2 <- DataMSub$depVar.h
          DataMSub$depVar.h2[DataMSub$Ind == 2] <- NA

          myFormula.ols <- paste0("depVar.h2 ~ ", paste(controls.info, collapse = "+"))
          OLS.mod <- lm(as.formula(myFormula.ols), data = DataMSub, na.action = "na.exclude")

          eti <- residuals(OLS.mod)
          eti[DataMSub$Event != 1] <- NA

          vti <- residuals(OLS.mod)
          vti[DataMSub$NoEvent != 1] <- NA

          # Pad residuals back to full sample length for consistent indexing
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

          OLSRes[[paste0("OLS.n", i)]] <- OLS.mod
        }

        irfest[h, i, e] <- IV.mod$coefficients["shockVar"] * norm
        irfse[h, i, e]  <- IV.se["shockVar"] * norm

        if (details == TRUE) {
          IVRes[[paste0("IV.h", h, ".n", i, ".e", e)]] <- IV.mod
        }
      }

    }
  }

  dimnames(irfest)[[1]] <- HSeries - 1
  dimnames(irfse)[[1]]  <- HSeries - 1

  Method <- "Proxy-IV"

  if (details == TRUE) {
    Sig  <- var(et, use = "complete.obs")
    SigR <- if (sum(!is.na(vt)) > 0) var(vt, use = "complete.obs") else NA
    Psi  <- irfest[1, , ]

    # Data for Lewis-Mertens (2025) weak instrument test
    if (controls.info[1] != "1") {
      WeakData <- data.frame(DataM[, paste0("y", 1:N)], DataM[, paste0("Z", 1:E)], DataM[, controls.info])
    } else {
      WeakData <- data.frame(DataM[, paste0("y", 1:N)], DataM[, paste0("Z", 1:E)])
    }
    if (E == 1) {
      colnames(WeakData) <- if (controls.info[1] != "1") c("y1", "Z1", controls.info) else c("y1", "Z1")
    }
    WeakData[DataM$Ind == 2, ] <- NA

    Obs <- data.frame(Tp = Te, Tc = Tn, To = To, Tt = Tt)

    return(list(irf = irfest, se = irfse,
                IVRes = IVRes, OLSRes = OLSRes,
                Obs = Obs, Method = Method,
                et = as.matrix(et), Sig = Sig, SigR = SigR, Psi = Psi,
                WeakData = WeakData))
  } else {
    return(list(irf = irfest, se = irfse, Method = Method))
  }
}
