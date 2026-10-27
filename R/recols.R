#' Estimate impulse responses via recursive (Cholesky) local projections
#'
#' Estimates impulse response functions (IRFs) using recursive zero
#' restrictions combined with local projections (Jorda, 2005). The shock of
#' dimension `e` is identified by regressing the outcome at horizon `h` by OLS
#' on the contemporaneous value of the `e`-th variable in `y`, controlling for
#' the contemporaneous values of variables `1, ..., e - 1` in `y`, lags of the
#' information set, and deterministic terms. This is the local-projection
#' analogue of a Cholesky decomposition and serves as a benchmark for the
#' heteroskedasticity-based estimator [hetiv()].
#'
#' @details For `E > 1`, identification is recursive and order-dependent: the
#'   column order of `y` defines both the shock ordering and the normalization
#'   variable for each shock dimension. By construction, the impact response
#'   of variable `e` to shock `e` equals `norm`, and the impact responses of
#'   variables `1, ..., e - 1` to shock `e` are zero.
#'
#'   Unlike [hetiv()], the distinction between policy days (`Ind == 1`) and
#'   control days (`Ind == 0`) does not affect the IRF estimates: both enter
#'   the local projections, and only contaminated days (`Ind == 2`) are
#'   dropped. The distinction matters only for the shock extraction returned
#'   when `details = TRUE`.
#'
#' @param y Numeric matrix of stationary outcome variables (T x N). The effect
#'   on the `e`-th variable in each dimension `e` is normalized to `norm` at
#'   horizon 0. The first `E` variables are also used to impose the recursive
#'   zero restrictions.
#' @param O Numeric matrix of information set variables (T x M). May be
#'   identical to `y`. Included as lags 1 through `P`.
#' @param X Numeric matrix of deterministic variables (T x K). For example,
#'   time trend, seasonal dummies or other deterministic controls. Included as
#'   is (no lags). A constant is included by default.
#' @param Ind Integer vector of length T, event indicator:
#'   \itemize{
#'     \item `0` Control day (no event)
#'     \item `1` Policy day (event)
#'     \item `2` Contaminated control day (excluded from estimation)
#'   }
#'   At least one policy day is required.
#' @param P Integer. Maximum lag order for the information set. Set to `0` for
#'   no lags (regression on deterministic terms only).
#' @param H Integer. Maximum horizon (in periods) up to which IRFs are estimated.
#' @param E Integer. Number of shock dimensions to identify via recursive
#'   ordering.
#' @param norm Numeric scalar. Normalize the impact response of the `e`-th
#'   variable to shock `e` to a specific value. Set to `1` for standard
#'   unit-effect normalization.
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
#'   Defaults to 0 (impact response and one-step-ahead residual). Must lie
#'   between 0 and `H - 1`.
#' @param pooled Logical. If `FALSE` (default), the shock prediction uses the
#'   forecast errors and responses at horizon `hpredict` only. If `TRUE`, the
#'   forecast errors of all horizons `0, ..., hpredict` are stacked and used
#'   jointly (pooled prediction), which weights every horizon by its
#'   signal-to-noise ratio; see [kfpredict()]. With `hpredict = 0` both
#'   options coincide.
#' @param details Logical. If `TRUE`, saves detailed OLS results and extracts
#'   shocks, which is slightly slower. If `FALSE`, returns only impulse
#'   responses and standard errors (e.g. for bootstrap).
#'
#' @return A named list. If `details = FALSE`, it contains `irf`, `se`, and
#'   `Method`. If `details = TRUE`, it contains:
#'   \describe{
#'     \item{`irf`}{Array (H x N x E) of estimated impulse responses.}
#'     \item{`se`}{Array (H x N x E) of local-projection standard errors.}
#'     \item{`Shocks`}{Output of [kfpredict()]: predicted structural shocks.}
#'     \item{`OLSRes`}{List of `lm` model objects, one per horizon, variable,
#'       and shock dimension.}
#'     \item{`Obs`}{Data frame with observation counts: `Tp` (policy days),
#'       `Tc` (control days), `To` (contaminated days), `Tt` (total used).}
#'     \item{`Method`}{Character string `"Recursive-OLS"`.}
#'     \item{`et`}{Matrix of forecast errors on event days used for shock
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
#'   }
#'
#' @seealso [hetiv()] for heteroskedasticity-based identification and
#'   [proxyiv()] for proxy-based identification.
#'
#' @references
#' Jorda, O. (2005). Estimation and inference of impulse responses by local
#' projections. *American Economic Review*, 95(1), 161-182.
#'
#' Montiel Olea, J. L., M. Plagborg-Moller, E. Qian, and C. K. Wolf (2025).
#' Local projections or VARs? A primer for macroeconomists. *NBER Working
#' Paper* No. 33871.
#'
#' Plagborg-Moller, M. and C. K. Wolf (2021). Local projections and VARs
#' estimate the same impulse responses. *Econometrica*, 89(2), 955-980.
#'
#' @importFrom dplyr lag lead
#' @importFrom sandwich NeweyWest vcovHC
#'
#' @examples
#' set.seed(1)
#' y <- matrix(rnorm(80), ncol = 2)
#' Ind <- rep(0L, nrow(y))
#' Ind[seq(5, nrow(y), by = 5)] <- 1L
#' res <- recols(y = y, O = y, Ind = Ind, P = 1, H = 3, E = 2)
#' res$irf[1, , ] # impact responses: unit diagonal, zero above it
#'
#' @export
recols <- function(y, O, X = NULL, Ind, P, H, E = 1, norm = 1,
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
  beg <- P + 1 # Start of the sample (end changes with every horizon)
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

  # Set up data set and objects to save results
  if (is.null(X)) {
    DataM <- data.frame(y, O, Ind)
    colnames(DataM) <- c(paste0("y", 1:N), paste0("o", 1:M), "Ind")
  } else {
    DataM <- data.frame(y, O, X, Ind)
    colnames(DataM) <- c(paste0("y", 1:N), paste0("o", 1:M), paste0("x", 1:K), "Ind")
  }
  irfest <- array(NA, dim = c(HNum, N, E))
  irfse <- array(NA, dim = c(HNum, N, E))
  OLSRes <- list()

  # Compute lagged variables included in information set O(t-1)...O(t-P)
  infoVars <- c()
  if (P > 0) {
    for (j in 1:M) {
      for (p in 1:P) {
        DataM[, paste0("o", j, ".l", p)] <- dplyr::lag(DataM[, paste0("o", j)], p)
        infoVars <- c(infoVars, paste0("o", j, ".l", p))
      }
    }
  }

  # Otherwise, only regress on a constant
  if (P == 0) {
    infoVars <- "1"
  }

  # Add deterministic variables
  if (!is.null(X)) {
    for (k in 1:K) {
      infoVars <- c(infoVars, paste0("x", k))
    }
  }

  # Set up event days (policy event, control, and contaminated event)
  DataM$Event <- (DataM$Ind == 1)
  DataM$NoEvent <- (DataM$Ind == 0)
  DataM$OthEvent <- (DataM$Ind == 2)

  # Observation counts on the horizon-0 estimation sample
  DataMSub <- DataM[beg:Nobs, ]
  Te <- sum(DataMSub$Event, na.rm = TRUE)
  Tn <- sum(DataMSub$NoEvent, na.rm = TRUE)
  To <- sum(DataMSub$OthEvent, na.rm = TRUE)
  Tt <- Te + Tn

  # Estimate the IRFs separately for every shock dimension
  for (e in 1:E) {
    # The shock variable is the e-th variable in y
    DataM$shockVar <- DataM[, paste0("y", e)]

    # Recursive zero restrictions: control for contemporaneous y_1, ..., y_(e-1)
    recVars <- c()
    if (e > 1) {
      for (q in 1:(e - 1)) {
        recVars <- c(recVars, paste0("y", q))
      }
    }

    # Control variables for the residual regression and for the local projections
    controls.info <- unique(c(infoVars))
    controls.info <- controls.info[controls.info != ""]

    controls.lp <- unique(c(infoVars, recVars))
    controls.lp <- controls.lp[controls.lp != ""]

    # Estimate the impulse responses by OLS for every variable in y
    for (i in 1:N) {
      DataM$depVar <- DataM[, paste0("y", i)]
      cumi <- cum[i]

      for (h_idx in seq_along(HSeries)) {
        h <- HSeries[h_idx]

        # Dependent variable at horizon h (level or cumulative response)
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

        # Shorten data to the subset without missing values at this horizon
        end <- Nobs - h + 1
        DataMSub <- DataM[beg:end, ]

        # LP (Jorda, 2005), including recursive control variables for e > 1
        myFormula <- paste0(
          "depVar.h ~ shockVar + ",
          paste(controls.lp, collapse = "+")
        )

        # Estimate OLS regression (excluding contaminated days) and LP standard errors
        # At h = 0 the regressions of y_1, ..., y_e on the shock and recursive
        # controls fit perfectly by construction, so the "essentially perfect
        # fit" warning from summary.lm() is expected and muffled
        OLS.mod <- lm(as.formula(myFormula), data = subset(DataMSub, Ind < 2))
        OLS.vcov <- withCallingHandlers(
          if (cov_type == "NW") {
            sandwich::NeweyWest(OLS.mod, prewhite = FALSE, adjust = TRUE)
          } else {
            sandwich::vcovHC(OLS.mod, type = "HC3")
          },
          warning = function(w) {
            if (grepl("essentially perfect fit", conditionMessage(w))) {
              invokeRestart("muffleWarning")
            }
          }
        )
        OLS.se <- sqrt(diag(OLS.vcov))

        # Normalize IRFs (impact response of variable e to shock e equals norm)
        irfest[h_idx, i, e] <- OLS.mod$coefficients["shockVar"] * norm
        irfse[h_idx, i, e] <- OLS.se["shockVar"] * abs(norm)

        if (details == TRUE) {
          # Save LP results for every variable, horizon, and shock dimension
          OLSRes[[paste0("OLS.h", h, ".n", i, ".e", e)]] <- OLS.mod

          # Residuals of the information-set regression at the horizons used
          # for shock prediction (event and control days; once, for e = 1)
          if (e == 1 && h %in% hres) {
            sr <- .shock_residuals(DataMSub, controls.info, DataM$Ind, beg, end, Nobs, h - 1)
            lab <- paste0("y", i, ".h", h - 1)
            etList[[lab]] <- sr$et
            vtList[[lab]] <- sr$vt
          }
        }
      }
    }
  }

  # Label rows of impulse responses to start at 0 (immediate response)
  dimnames(irfest)[[1]] <- HSeries - 1
  dimnames(irfse)[[1]] <- HSeries - 1

  Method <- "Recursive-OLS"

  if (details == TRUE) {
    # Stack residuals and responses over the prediction horizons, compute the
    # covariance matrices on event and control days, and predict the shocks
    sk <- .stack_shock_inputs(etList, vtList, irfest, hres, N, E)
    et <- sk$et
    Sig <- sk$Sig
    SigR <- sk$SigR
    Psi <- sk$Psi
    
    # Check whether SigR is NA (not NA if no control days). In this case, do
    # not scale the shock using control days
     if(any(is.na(SigR))){
       SigR <- 0
      scale <- FALSE
     }else{
      scale <- TRUE
    }
    Shocks <- kfpredict(Sig, SigR, Psi, et, scale = scale)
    
    Obs <- data.frame(Tp = Te, Tc = Tn, To = To, Tt = Tt)

    return(list(
      irf = irfest, se = irfse,
      Shocks = Shocks,
      OLSRes = OLSRes,
      Obs = Obs, Method = Method,
      et = as.matrix(et), Sig = Sig, SigR = SigR, Psi = Psi
    ))
  } else {
    return(list(irf = irfest, se = irfse, Method = Method))
  }
}
