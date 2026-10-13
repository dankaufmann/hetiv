#' Estimate impulse responses via proxy-IV local projections
#'
#' Estimates impulse response functions (IRFs) using user-provided external
#' instruments (proxies) combined with local projections (Jorda, 2005). The
#' proxy variables serve directly as instruments for the endogenous shock
#' variables. Optionally imposes recursive zero restrictions across shock
#' dimensions and supports deterministic controls following the same
#' conventions as [hetiv()].
#'
#' @details For `E > 1`, identification can be order-dependent: the column order
#'   of `y` defines the endogenous shock variables and normalizations, the
#'   column order of `Z` assigns proxies to shock dimensions, and
#'   `recursive = TRUE` imposes restrictions in that order.
#'
#' @param y Numeric matrix of stationary outcome variables (T x N). The effect
#'   on the first variable in each dimension is normalized to `norm` at horizon
#'   0. These variables are also used as the endogenous regressors instrumented
#'   by the columns of `Z`.
#' @param O Numeric matrix of information set variables (T x M). May be
#'   identical to `y`. Included as lags 1 through `P`.
#' @param Z Numeric matrix of external instruments (T x E). Column `e` is used
#'   as the proxy for shock dimension `e`. Following the proxy-SVAR literature
#'   (Mertens and Ravn, 2013), the recommended convention is to set the proxy
#'   to zero on control days (`Ind == 0`): these days then remain in the IV
#'   sample and sharpen the estimates of the control coefficients, while the
#'   identifying moments are unaffected. Alternatively, missing values on
#'   control days are treated as contaminated and excluded from estimation,
#'   which reproduces the event-day sample of the event-study literature
#'   (Nakamura and Steinsson, 2018). With zeros, prefer the robust first-stage
#'   statistics (`F_rob`, `WeakTest`) over the homoskedastic F, which treats
#'   the zero observations as informative. Missing values on policy days
#'   (`Ind == 1`) remain labelled as events for residual-based shock
#'   prediction, but are necessarily dropped from IV regressions that use the
#'   missing proxy.
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
#' @param cov_type Covariance estimator for local-projection standard errors:
#'   `"HC3"` (default) for heteroskedasticity-robust standard errors or `"NW"`
#'   for Newey-West HAC standard errors. `"HC3"` is the default because
#'   Montiel Olea et al. (2025) show that heteroskedasticity-robust standard
#'   errors suffice for local-projection impulse responses under weak
#'   conditions, even though multi-step forecast errors are typically serially
#'   correlated. `"NW"` remains available as an optional HAC robustness check.
#' @param recursive Logical. If `TRUE`, imposes recursive zero restrictions
#'   across shock dimensions: for shock `e > 1`, the variables and instruments
#'   from dimensions `1, ..., e-1` are added as controls. Default `FALSE`.
#' @param hpredict Integer. Horizon of the responses and forecast errors used
#'   for shock prediction, on the same scale as the row labels of `irf`.
#'   Defaults to 0 (impact response and one-step-ahead residual). Must lie
#'   between 0 and `H - 1`.
#' @param pooled Logical. If `FALSE` (default), the shock prediction uses the
#'   forecast errors and responses at horizon `hpredict` only. If `TRUE`, the
#'   forecast errors of all horizons `0, ..., hpredict` are stacked and used
#'   jointly (pooled prediction), which weights every horizon by its
#'   signal-to-noise ratio; see [kfpredict()]. With `hpredict = 0` both
#'   options coincide.
#' @param details Logical. If `TRUE`, returns detailed results including IV
#'   model objects, OLS residuals, and covariance matrices. If `FALSE` (default),
#'   returns only impulse responses and standard errors (faster; use for
#'   bootstrap).
#'
#' @return A named list. Always contains:
#'   \describe{
#'     \item{`irf`}{Array (H x N x E) of estimated impulse responses.}
#'     \item{`se`}{Array (H x N x E) of local-projection standard errors.}
#'     \item{`Method`}{Character string `"Proxy-IV"`.}
#'   }
#'   With `details = TRUE`, additionally contains:
#'   \describe{
#'     \item{`Shocks`}{Predicted structural shocks on event days from
#'       [kfpredict()], based on the responses at horizon `hpredict`.}
#'     \item{`IVRes`}{List of `ivreg` model objects, one per horizon, variable,
#'       and shock dimension.}
#'     \item{`OLSRes`}{List of the information-set regressions whose residuals
#'       enter the shock prediction, named `OLS.n<i>.h<h>` by outcome
#'       variable and horizon.}
#'     \item{`Obs`}{Data frame with observation counts: `Tp` (policy days),
#'       `Tc` (control days), `To` (contaminated days), `Tt` (labelled event or
#'       control days), and `Tiv` (complete observations available to the IV
#'       regression of the first shock at the impact horizon).}
#'     \item{`FirstStage`}{For `E = 1`: data frame with one row per horizon
#'       `h` giving the number of observations `Tobs`, the homoskedastic
#'       first-stage F statistic `F_hom` (the squared t statistic of the
#'       proxy; Stock and Yogo, 2005) and the robust effective F statistic
#'       `F_rob` of Montiel Olea and Pflueger (2013), computed with the
#'       covariance estimator in `cov_type`. The statistics do not depend on
#'       the outcome variable. `NULL` for `E > 1`.}
#'     \item{`WeakTest`}{For `E > 1`: data frame with one row per horizon `h`
#'       and outcome variable `i` giving the results of [gweakivtest()] for
#'       the system with `E` endogenous policy variables and `E` proxies:
#'       `gmin` (generalized minimum eigenvalue statistic, which depends on
#'       the horizon only through the sample), `cv` and `cv_simplified`
#'       (sharp and simplified critical values, which depend on the outcome
#'       and horizon), and the Stock-Yogo statistic `sy_stat` with its
#'       critical value `sy_cv` (`NaN` unless over-identified). Rows where
#'       the outcome is a policy variable on impact are omitted. `NULL` for
#'       `E = 1`.}
#'     \item{`et`}{Data frame of forecast errors on event days used for shock
#'       prediction, one column per variable and horizon (`y1.h0`, ...,
#'       `yN.h0`, `y1.h1`, ...): horizon `hpredict` only if `pooled = FALSE`,
#'       horizons `0, ..., hpredict` if `pooled = TRUE`.}
#'     \item{`Sig`}{Covariance matrix of `et` on event days.}
#'     \item{`SigR`}{Covariance matrix of the same residuals on control days,
#'       or `NA` if unavailable. A control day enters horizon h only if no
#'       event occurs in its window `t + 1, ..., t + h`.}
#'     \item{`Psi`}{Response matrix used for shock prediction, with rows
#'       matching the columns of `et`: `irf[hpredict + 1, , ]` (N x E) if
#'       `pooled = FALSE`, the rows `1, ..., hpredict + 1` of `irf` stacked
#'       if `pooled = TRUE`. With the default `hpredict = 0` this is the
#'       impact response.}
#'     \item{`WeakData`}{Data frame of endogenous variables and instruments for
#'       the Lewis-Mertens (2025) weak instrument test.}
#'   }
#'
#' @references
#' Jorda, O. (2005). Estimation and inference of impulse responses by local
#' projections. *American Economic Review*, 95(1), 161-182.
#'
#' Lewis, D. J. and Mertens, K. (2025). A robust test for weak instruments for
#' 2SLS with multiple endogenous regressors. *Review of Economic Studies*,
#' DOI: 10.1093/restud/rdaf103
#'
#' Mertens, K. and Ravn, M. O. (2013). The dynamic effects of personal and
#' corporate income tax changes in the United States. *American Economic
#' Review*, 103(4), 1212-1247.
#'
#' Montiel Olea, J. L., M. Plagborg-Moller, E. Qian, and C. K. Wolf (2025).
#' Local projections or VARs? A primer for macroeconomists. *NBER Working
#' Paper* No. 33871.
#'
#' Nakamura, E. and Steinsson, J. (2018). High-frequency identification of
#' monetary non-neutrality: The information effect. *Quarterly Journal of
#' Economics*, 133(3), 1283-1330.
#'
#' Stock, J. H. and Watson, M. W. (2018). Identification and estimation of
#' dynamic causal effects in macroeconomics using external instruments.
#' *Economic Journal*, 128(610), 917-948.
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
#' Z <- matrix(Ind * y[, 1] + rnorm(nrow(y)), ncol = 1)
#' res <- proxyiv(
#'   y = y, O = y, Z = Z, Ind = Ind, P = 1, H = 3,
#'   details = TRUE
#' )
#' dim(res$irf)
#'
#' @export
proxyiv <- function(y, O, Z, X = NULL, Ind, P, H, E = 1, norm = 1,
                    cum = FALSE, cov_type = "HC3",
                    recursive = FALSE, hpredict = 0,
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
  Z <- .as_numeric_matrix(Z, "Z", nrow = nrow(y))
  recursive <- .check_logical_scalar(recursive, "recursive")
  hpredict <- args$hpredict
  hrow <- hpredict + 1L # row of the IRF array holding horizon hpredict
  pooled <- args$pooled
  # IRF rows whose residuals enter the shock prediction
  hres <- if (pooled) seq_len(hrow) else hrow
  etList <- list()
  vtList <- list()
  details <- .check_logical_scalar(details, "details")

  # Collect properties of the data
  Nobs <- dim(y)[1]
  beg <- P + 1
  #end <- Nobs - H + 1 # End of the sample (moved. Change with every horizon)
  N <- dim(y)[2]
  M <- dim(O)[2]
  K <- if (!is.null(X)) dim(X)[2] else 0

  if (ncol(Z) < E) {
    stop("Z must have at least E columns (one instrument per shock dimension).",
      call. = FALSE
    )
  }

  if (sum(is.na(y)) > 0) warning("Missing values in y")
  if (sum(is.na(O)) > 0) warning("Missing values in O")
  if (sum(is.na(Ind)) > 0) warning("Missing values in Ind")
  if (!is.null(X) && sum(is.na(X)) > 0) warning("Missing values in X")

  # Horizons at which IRFs are computed: h = 1..H corresponds to horizon h - 1
  HSeries <- seq_len(H)
  HNum <- H

  # Build main data frame
  if (is.null(X)) {
    DataM <- data.frame(y, O, Ind, Z)
    colnames(DataM) <- c(paste0("y", 1:N), paste0("o", 1:M), "Ind", paste0("z", 1:E))
  } else {
    DataM <- data.frame(y, O, X, Ind, Z)
    colnames(DataM) <- c(
      paste0("y", 1:N), paste0("o", 1:M), paste0("x", 1:K),
      "Ind", paste0("z", 1:E)
    )
  }

  irfest <- array(NA, dim = c(HNum, N, E))
  irfse <- array(NA, dim = c(HNum, N, E))
  IVRes <- list()
  OLSRes <- list()
  FirstStage <- NULL # first-stage F statistics by horizon (E == 1, details)
  WeakTest <- NULL # Lewis-Mertens tests by horizon and outcome (E > 1, details)
  Tiv <- NA

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
    DataM$Ind <- Ind

    # Endogenous variable: the e-th variable in y
    DataM$shockVar <- DataM[, paste0("y", e)]

    # Instrument: the e-th column of user-provided Z
    DataM$Ze <- DataM[, paste0("z", e)]
    DataM$Ind[DataM$Ind == 0 & is.na(DataM$Ze)] <- 2

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
    controls.lp <- unique(c(infoVars, recVars))
    controls.lp <- controls.lp[controls.lp != ""]
    controls.iv <- unique(c(infoVars, recInst))
    controls.iv <- controls.iv[controls.iv != ""]

    DataM$Event <- (DataM$Ind == 1)
    DataM$NoEvent <- (DataM$Ind == 0)
    DataM$OthEvent <- (DataM$Ind == 2)

    end      <- Nobs
    DataMSub <- DataM[beg:end, ]

    Te <- sum(DataMSub$Event)
    Tn <- sum(DataMSub$NoEvent)
    To <- sum(DataMSub$OthEvent)
    Tt <- Te + Tn
    if (Te == 0) {
      stop("At least one event day (Ind == 1) is required.")
    }

    # Use the user-provided proxy directly as instrument (no orthogonalization;
    # 2SLS projection in ivreg handles the required partialling out)
    DataMSub$Z <- DataMSub$Ze
    DataM[beg:end, paste0("Z", e)] <- DataMSub$Z
    DataM[beg:end, "Z"] <- DataMSub$Z

    # Estimate impulse responses for every outcome variable in y
    for (i in 1:N) {
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

        end      <- Nobs - h + 1       # End of the sample (moved. Change with every horizon)
        DataMSub <- DataM[beg:end, ]

        # Proxy-IV LP (Jorda, 2005): instrument shockVar with Z, control for information set
        myFormula <- paste0(
          "depVar.h ~ shockVar + ",
          paste(controls.lp, collapse = "+"),
          " | Z + ",
          paste(controls.iv, collapse = "+")
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
          if (e == 1 && i == 1 && h == 1) {
            # Complete observations available to the IV regression on impact
            iv_vars <- unique(c("depVar.h", "shockVar", "Z", controls.lp, controls.iv))
            iv_vars <- iv_vars[iv_vars != "1"]
            Tiv <- sum(complete.cases(DataIV[, iv_vars, drop = FALSE]))
          }
          if (E == 1 && i == 1) {
            # One instrument: homoskedastic F (= t^2) and robust effective F
            # (Montiel Olea and Pflueger, 2013); identical for every outcome
            fs <- .firststage_f(DataIV$shockVar, DataIV$Z, X = Xctrl, cov_type = cov_type)
            FirstStage <- rbind(FirstStage, data.frame(
              h = h - 1, Tobs = fs$Tobs, F_hom = fs$F_hom, F_rob = fs$F_rob
            ))
          }
          if (E > 1 && e == E && !(h == 1 && i <= E)) {
            # Several instruments: generalized test of Lewis and Mertens (2025)
            # for the full system, whose critical values depend on the outcome
            # at horizon h (skipped where the outcome is itself a policy
            # variable on impact)
            WeakTest <- rbind(WeakTest, .lm_weaktest_row(
              y = DataIV$depVar.h,
              Y = as.matrix(DataIV[, paste0("y", seq_len(E)), drop = FALSE]),
              X = Xctrl,
              Z = as.matrix(DataIV[, paste0("Z", seq_len(E)), drop = FALSE]),
              h = h - 1, i = i, cov_type = cov_type
            ))
          }
        }

        # Residuals of the information-set regression at the horizons used for
        # shock prediction (event and control days; once, for e = 1)
        if (details == TRUE && e == 1 && h %in% hres) {
          sr <- .shock_residuals(DataMSub, controls.info, DataM$Ind, beg, end, Nobs, h - 1)
          lab <- paste0("y", i, ".h", h - 1)
          etList[[lab]] <- sr$et
          vtList[[lab]] <- sr$vt
          OLSRes[[paste0("OLS.n", i, ".h", h - 1)]] <- sr$mod
        }

        irfest[h_idx, i, e] <- IV.mod$coefficients["shockVar"] * norm
        irfse[h_idx, i, e] <- IV.se["shockVar"] * abs(norm)

        if (details == TRUE) {
          IVRes[[paste0("IV.h", h, ".n", i, ".e", e)]] <- IV.mod
        }
      }
    }
  }

  dimnames(irfest)[[1]] <- HSeries - 1
  dimnames(irfse)[[1]] <- HSeries - 1
  
  Method <- "Proxy-IV"

  if (details == TRUE) {
    # Stack residuals and responses over the prediction horizons, compute the
    # covariance matrices on event and control days, and predict the shocks
    sk <- .stack_shock_inputs(etList, vtList, irfest, hres, N, E)
    et <- sk$et
    Sig <- sk$Sig
    SigR <- sk$SigR
    Psi <- sk$Psi
    Shocks <- kfpredict(Sig, SigR, Psi, et)
    
    # Data for Lewis-Mertens (2025) weak instrument test
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
    WeakData[DataM$Ind == 2, ] <- NA

    Obs <- data.frame(Tp = Te, Tc = Tn, To = To, Tt = Tt, Tiv = Tiv)

    return(list(
      irf = irfest, se = irfse,
      Shocks = Shocks,
      IVRes = IVRes, OLSRes = OLSRes,
      Obs = Obs, Method = Method,
      et = as.matrix(et), Sig = Sig, SigR = SigR, Psi = Psi,
      WeakData = WeakData,
      FirstStage = FirstStage, WeakTest = WeakTest
    ))
  } else {
    return(list(irf = irfest, se = irfse, Method = Method))
  }
}
