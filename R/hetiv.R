#' Estimate impulse responses via heteroskedasticity-based IV local projections
#'
#' Estimates impulse response functions (IRFs) using recursive
#' heteroskedasticity-IV identification (Rigobon, 2003; Rigobon and Sack, 2004;
#' Lewis, 2022; Burri and Kaufmann, 2026a, 2026b) combined with local projections
#' (Jorda, 2005). Identification exploits the difference in variance between
#' policy event days and control days to construct instruments for the
#' endogenous variables. It also returns weak instrument tests as well as shock
#' predictions on event days. 
#'
#' @details For `E > 1`, identification is recursive and order-dependent: the
#'   column order of `y` defines both the shock ordering and the normalization
#'   variable for each shock dimension.
#'
#' @param y Numeric matrix of outcome variables (T x N). The effect on
#'   the first variable in each dimension is normalized to unity at horizon 0.
#'   These variables are also used to construct heteroskedasticity-based instruments
#'   and for recursive ordering to identify multiple dimensions.
#' @param O Numeric matrix of information set variables (T x M). Should include at least
#'   `y` in case of non-stationary variables. Included as lags 1 through `P`.
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
#'   variable to a specific value. Defaults to `1` for standard unit-effect
#'   normalization.
#' @param interact Logical. If `TRUE`, lagged information set variables are
#'   interacted with event/non-event dummies.
#' @param cum Logical vector of length N. For each variable in `y`, whether
#'   to report the cumulative impulse response instead of the level response.
#'   If only one provided, applied to all impulse responses.
#' @param cov_type Covariance estimator for local-projection standard errors:
#'   `"HC3"` (default) for heteroskedasticity-robust standard errors or `"NW"`
#'   for Newey-West HAC standard errors. `"HC3"` is the default because
#'   Montiel Olea et al. (2025) show that heteroskedasticity-robust standard
#'   errors suffice for local-projection impulse responses under weak
#'   conditions, even though multi-step forecast errors are typically serially
#'   correlated. `"NW"` remains available as an optional HAC robustness check.
#' @param hpredict Integer. Horizon of the responses and forecast errors used
#'   for shock prediction, on the same scale as the row labels of `irf`.
#'   Defaults to 0 (recommended; impact response and one-step-ahead residual are used)
#'   Must lie between 0 and `H - 1`.
#' @param pooled Logical. If `FALSE` (default), the shock prediction uses the
#'   forecast errors and responses at horizon `hpredict` only. If `TRUE`, the
#'   forecast errors of all horizons `0, ..., hpredict` are stacked and used
#'   jointly (pooled prediction); see [kfpredict()]. With `hpredict = 0` both
#'   options coincide.
#' @param details Logical. If `TRUE`, saves detailed IV results, weak instrument
#'   tests, and shock prediction results, which is slightly slower. If `FALSE`, 
#'   returns only impulse responses and standard errors.
#'
#' @return A named list. Always contains:
#'   \describe{
#'     \item{`irf`}{Array (H x N x E) of estimated impulse responses.}
#'     \item{`se`}{Array (H x N x E) of local-projection standard errors.}
#'     \item{`Method`}{Character string `"Heteroskedasticity-IV"`.}
#'   }
#'   With `details = TRUE`, additionally contains:
#'   \describe{
#'     \item{`irf`}{Array (H x N x E) of estimated impulse responses.}
#'     \item{`se`}{Array (H x N x E) of local-projection standard errors.}
#'     \item{`Shocks`}{Predicted structural shocks on event days from
#'       [kfpredict()], based on the responses at horizon `hpredict`.}
#'     \item{`IVRes`}{List of `ivreg` model objects, one per horizon, variable,
#'       and shock dimension.}
#'     \item{`OLSRes`}{List of the information-set regressions whose residuals
#'       enter the shock prediction, named `OLS.n<i>.h<h>` by outcome
#'       variable and horizon.}
#'     \item{`ORTHRes`}{List of OLS models used to orthogonalize the variable
#'       underlying the instrument, one per shock dimension.}
#'     \item{`FirstStage`}{For `E = 1`: data frame with one row per horizon
#'       `h` giving the number of observations `Tobs`, the homoskedastic
#'       first-stage F statistic `F_hom` (the squared t statistic of the
#'       proxy; Stock and Yogo, 2005) and the robust effective F statistic
#'       `F_rob` of Montiel Olea and Pflueger (2013), always computed with HAC
#'       robust covariance estimator. The statistics do not depend on
#'       the outcome variable. `NULL` for `E > 1`.}
#'     \item{`WeakTest`}{For `E > 1`: data frame with one row per horizon `h`
#'       and outcome variable `i` and shock dimension `e` 
#'       giving the results of [gweakivtest()] for
#'       the system with `E` endogenous policy variables and `E` instruments:
#'       `gmin` (generalized minimum eigenvalue statistic, which depends on
#'       the horizon only through the sample), `cv` and `cv_simplified`
#'       (sharp and simplified critical values, which depend on the outcome
#'       and horizon), and the Stock-Yogo statistic `sy_stat` with its
#'       critical value `sy_cv` (`NaN` unless over-identified). Rows where
#'       the outcome is a policy variable on impact are omitted. `NULL` for
#'       `E = 1`. Note that the test statistic is always based on HAC robust 
#'        covariance, the target is always the `e` endogenous variable,
#'        and the relative criterion is reported for comparability 
#'        with Montiel Olea and Pflueger (2013).}
#'     \item{`Obs`}{Data frame with observation counts: `Tp` (policy days),
#'       `Tc` (control days), `To` (contaminated days), `Tt` (total used).}
#'     \item{`Method`}{Character string `"Heteroskedasticity-IV"`.}
#'     \item{`et`}{Data frame of forecast errors on event days used for shock
#'       prediction: residuals of the regression of `y_{t+h}` on the
#'       information set dated `t - 1`, one column per variable and horizon
#'       (`y1.h0`, ..., `yN.h0`, `y1.h1`, ...). One horizon (`hpredict`) if
#'       `pooled = FALSE`, horizons `0, ..., hpredict` if `pooled = TRUE`.}
#'     \item{`Sig`}{Covariance matrix of `et` on event days.}
#'     \item{`SigR`}{Covariance matrix of the same residuals on control days,
#'       or `NA` if unavailable. A control day enters horizon h only if no
#'       event occurs in its window `t + 1, ..., t + h`.}
#'     \item{`Psi`}{Response matrix used for shock prediction, with rows
#'       matching the columns of `et`: `irf[hpredict + 1, , ]` (N x E) if
#'       `pooled = FALSE`, the rows `1, ..., hpredict + 1` of `irf` stacked
#'       (N(hpredict + 1) x E) if `pooled = TRUE`. With the default
#'       `hpredict = 0` this is the impact response.}
#'     \item{`WeakData`}{Data frame of endogenous variables and instruments for
#'       the Lewis-Mertens (2025) weak instrument test.}
#'   }
#'
#' @references
#' Burri, M. and D. Kaufmann (2026a). Measuring monetary policy shocks.
#' IRENE Working Papers 24-03, IRENE Institute of Economic Research,
#' University of Neuchatel.
#'
#' Burri, M. and D. Kaufmann (2026b). Multiple monetary policy shocks from
#' daily data: A heteroskedasticity IV approach. *Economics Letters*, 268,
#' 113091. \doi{10.1016/j.econlet.2026.113091}
#'
#' Jorda, O. (2005). Estimation and inference of impulse responses by local
#' projections. *American Economic Review*, 95(1), 161-182.
#'
#' Lewis, D. J. (2022). Robust inference in models identified via
#' heteroskedasticity. *Review of Economics and Statistics*, 104(3), 510-524.
#'
#' Lewis, D. J. and Mertens, K. (2025). A robust test for weak instruments for 2SLS with multiple
#' endogenous regressors. *Review of Economic Studies*, DOI: 10.1093/restud/rdaf103
#'
#' Montiel Olea, J. L., M. Plagborg-Moller, E. Qian, and C. K. Wolf (2025).
#' Local projections or VARs? A primer for macroeconomists. *NBER Working
#' Paper* No. 33871.
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
                  cum = FALSE, cov_type = "HC3", hpredict = 0,
                  pooled = FALSE, details = FALSE) {
  args <- .validate_estimator_inputs(
    y = y, O = O, X = X, Ind = Ind, P = P, H = H, E = E, norm = norm,
    cum = cum, cov_type = cov_type, hpredict = hpredict, pooled = pooled
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
  cov_type <- args$cov_type
  interact <- .check_logical_scalar(interact, "interact")
  hpredict <- args$hpredict
  hrow <- hpredict + 1L # row of the IRF array holding horizon hpredict
  pooled <- args$pooled
  
  # IRF rows whose residuals enter the shock prediction
  hres <- if (pooled) seq_len(hrow) else hrow
  etList <- list()
  vtList <- list()
  details <- .check_logical_scalar(details, "details")

  # Collect various properties of the data and observations to be used
  Nobs <- dim(y)[1] # Number of observations in y
  beg <- P + 1 # Start of the sample
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

  # Horizons at which IRFs are computed: h = 1..H corresponds to horizon h - 1
  HSeries <- seq_len(H)
  HNum <- H

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
  FirstStage <- NULL # first-stage F statistics by horizon (E == 1, details)
  WeakTest <- NULL # Lewis-Mertens tests by horizon and outcome (E > 1, details)

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
    end      <- Nobs             # End of the sample 
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
    
    # Only save residuals when Ind < 2, because contaminated events (Ind == 2) are not used in the estimation
    # Make sure that other are missing
    DataMSub$Ze[DataMSub$Ind < 2]  <- residuals(orthModel, na.action = "na.exclude")
    DataMSub$Ze[DataMSub$Ind == 2] <- NA
    
    # Compute instrument and F-Statistic (see Lewis, 2022, and Rigobon and Sack, 2004)
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
        end      <- Nobs - h + 1      # End of the sample (change with every horizon)
        DataMSub <- DataM[beg:end, ]

        # Local projection IV (Jorda, 2005, Rigobon and Sack, 2004, Lewis, 2022): 
        # regress the dependent variable at horizon h
        # instrument shockVar with Z, control for information set
        myFormula <- paste0(
          "depVar.h ~ shockVar + ",
          paste(controls.lp, collapse = "+"),
          paste0(
            "| Z + ",
            paste(controls.iv, collapse = "+")
          )
        )

        # Estimate IV regression and compute local-projection standard errors
        DataIV <- subset(DataMSub, Ind < 2)
        IV.mod <- ivreg::ivreg(as.formula(myFormula), data = DataIV)
        if (cov_type == "NW") {
          IV.vcov <- sandwich::NeweyWest(IV.mod, prewhite = FALSE, adjust = TRUE)
        } else {
          IV.vcov <- sandwich::vcovHC(IV.mod, type = "HC3")
        }
        IV.se <- sqrt(diag(IV.vcov))

        # Weak-instrument diagnostics on the sample of the IV regression
        if (details == TRUE) {
          ctrl <- controls.info[controls.info != "1"]
          Xctrl <- if (length(ctrl) > 0) as.matrix(DataIV[, ctrl, drop = FALSE]) else NULL
          if (E == 1 && i == 1) {
            # One instrument: homoskedastic F (= t^2) and robust effective F
            # (Montiel Olea and Pflueger, 2013); identical for every outcome
            fs <- .firststage_f(DataIV$shockVar, DataIV$Z, X = Xctrl, cov_type = "NW")
            FirstStage <- rbind(FirstStage, data.frame(
              h = h - 1, Tobs = fs$Tobs, F_hom = fs$F_hom, F_rob = fs$F_rob
            ))
          }
          if (E > 1 && !(h == 1 && i <= E)) {
            # Several instruments: generalized test of Lewis and Mertens (2025)
            # for the full system, whose critical values depend on the outcome
            # at horizon h (skipped where the outcome is itself a policy
            # variable on impact)
            WeakTest <- rbind(WeakTest, .lm_weaktest_row(
              y = DataIV$depVar.h,
              Y = as.matrix(DataIV[, paste0("y", seq_len(e)), drop = FALSE]),
              X = Xctrl,
              Z = as.matrix(DataIV[, paste0("Z", seq_len(e)), drop = FALSE]),
              h = h - 1, i = i, cov_type = "NW", crit = "rel", target = e
            ))
          }
        }

        if (details == TRUE) {
          # Residuals of the information-set regression at the horizons used
          # for shock prediction (event and control days; once, for e = 1)
          if (e == 1 && h %in% hres) {
            sr <- .shock_residuals(DataMSub, controls.info, DataM$Ind, beg, end, Nobs, h - 1)
            lab <- paste0("y", i, ".h", h - 1)
            etList[[lab]] <- sr$et
            vtList[[lab]] <- sr$vt
            OLSRes[[paste0("OLS.n", i, ".h", h - 1)]] <- sr$mod
          }
        }

        # Normalizes IRFs to a specific value. Because initial response
        irfest[h_idx, i, e] <- IV.mod$coefficients["shockVar"] * norm
        irfse[h_idx, i, e] <- IV.se["shockVar"] * abs(norm)
        

        # Save IV results for every variable and every horizon
        if (details == TRUE) {
          IVRes[[paste0("IV.h", h, ".n", i, ".e", e)]] <- IV.mod

          if (h == 1) {
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
    # Stack residuals and responses over the prediction horizons, compute the
    # covariance matrices on event and control days, and predict the shocks
    sk <- .stack_shock_inputs(etList, vtList, irfest, hres, N, E)
    et <- sk$et
    Sig <- sk$Sig
    SigR <- sk$SigR
    Psi <- sk$Psi
    Shocks <- kfpredict(Sig, SigR, Psi, et)
    
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

  Method <- "Heteroskedasticity-IV"

  if (details == TRUE) {
    return(list(
      irf = irfest, se = irfse,
      Shocks = Shocks,
      IVRes = IVRes, OLSRes = OLSRes, ORTHRes = ORTHRes,
      Obs = Obs, Method = Method,
      et = as.matrix(et), Sig = Sig, SigR = SigR, Psi = Psi, WeakData = WeakData,
      FirstStage = FirstStage, WeakTest = WeakTest
    ))
  } else {
    return(list(irf = irfest, se = irfse, Method = Method))
  }
}
