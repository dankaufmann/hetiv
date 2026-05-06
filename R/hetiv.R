#' Estimate impulse responses via heteroskedasticity-based IV local projections
#'
#' Estimates impulse response functions (IRFs) using recursive
#' heteroskedasticity-IV identification (Rigobon, 2003; Rigobon and Sachs, 2004;
#' Lewis, 2022; Burri and Kaufmann, 2026) combined with local projections 
#' (Jordà, 2005). Identification exploits the difference in variance between 
#' policy event days and control days to construct instruments for the 
#' endogenous variables.
#'
#' @param y Numeric matrix of stationary outcome variables (T x N). The effect on 
#'   the first variable in each dimension is normalized to unity at horizon 0. 
#'   These variables are also used to construct heteroskedasticity-based instruments 
#'   and for recursive ordering to identify multiple dimensions.
#' @param O Numeric matrix of information set variables (T x M). May be
#'   identical to `y`. Included as lags 1 through `P`.
#' @param X Numeric matrix of deterministic variables (T x K). May include a constant, 
#'   time trend, seasonal dummies or other deterministic controls. Included as is (no lags).
#' @param Ind Integer vector of length T, event indicator:
#'   \itemize{
#'     \item `0` Control day (no event)
#'     \item `1` Policy day (event)
#'     \item `2` Contaminated control day (excluded from estimation)
#'   }
#' @param P Integer. Maximum lag order for the information set. Set to `0` for
#'   no lags (regression on constant only).
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
#' @param Hstep Integer. Step size between horizons. If `> 1`, only every
#'   `Hstep`-th horizon is estimated starting from `H = 0`.
#' @param details Logical. If `TRUE`, code saves detailed IV results, which is slower.
#'   if set to `FALSE`, returns only impulse response and standard error (e.g. for bootstrap)
#'
#' @return A named list with the following elements:
#'   \describe{
#'     \item{`irf`}{Array (H x N x E) of estimated impulse responses.}
#'     \item{`se`}{Array (H x N x E) of HC0 heteroscedasticity-robust standard
#'       errors.}
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
#'     \item{`Psi`}{Impact matrix (N x E), i.e. `irf[1, , ]`.}
#'     \item{`WeakData`}{Data frame of endogenous variables and instruments for
#'       the Lewis-Mertens (2025) weak instrument test.}
#'   }
#'
#' @references
#' Burri, M. and D. Kaufmann (2026). Measuring monetary policy shocks.
#' IRENE Working Papers 24-03, IRENE Institute of Economic Research,
#' University of Neuchâtel.
#'
#' Burri, M. and D. Kaufmann (2026). Multiple monetary policy shocks from
#' daily data: A heteroskedasticity IV approach. IRENE Working Papers 26-06,
#' IRENE Institute of Economic Research, University of Neuchâtel.
#'
#' Jordà, Ò. (2005). Estimation and Inference of Impulse Responses by Local
#' Projections. *American Economic Review*, 95(1), 161–182.
#'
#' Lewis, D. J. (2022). Robust Inference in Models Identified via
#' Heteroskedasticity. *Review of Economics and Statistics*, 104(3), 510–524.
#' 
#' Lewis, D. J. and Mertens, K. (2025). A Robust Test for Weak Instruments for 2SLS with Multiple
#' Endogenous Regressors. *The Review of Economic Studies*, DOI: 10.1093/restud/rdaf103
#'
#' Rigobon, R. (2003). Identification Through Heteroskedasticity.
#' *Review of Economics and Statistics*, 85(4), 777–792.
#'
#' Rigobon, R. and Sack, B. (2004). The impact of monetary policy on asset
#' prices. *Journal of Monetary Economics*, 51(8), 1553–1575.
#'
#' @importFrom dplyr lag lead
#' @importFrom ivreg ivreg
#' @importFrom sandwich vcovHC
#'
#' @export
hetiv <- function(y, O, X = NULL, Ind, P, H, E = 1, norm = 1, interact = FALSE, cum = FALSE, Hstep = 1, details = FALSE){

  # Collect various properties of the data and observations to be used
  Nobs <- dim(y)[1]     # Number of observations in y
  beg  <- max(P+1)      # Start of the sample
  end  <- Nobs - H+1    # End of the sample
  N    <- dim(y)[2]     # Number of variables in y
  M    <- dim(O)[2]     # Number of variables in O
  K    <- dim(X)[2]     # Number of deterministic variables in X

  if(sum(is.na(y))>0){
    warning("Missing values in y")
  }
  if(sum(is.na(O))>0){
    warning("Missing values in O")
  }
  if(sum(is.na(Ind))>0){
    warning("Missing values in Ind")
  }
  if(sum(is.na(X))>0){
    warning("Missing values in X")
  }

  # Modify cumulation indicator if only one provided
  if(length(cum) == 1){
    cum <- rep(cum, N)
  }

  # Specify at which horizons IRF should be computed
  HSeries <- seq(1, H, Hstep)
  HNum    <- length(HSeries)

  # Set up data set and and objects to save results
  if(is.null(X)){
    DataM     <- data.frame(y, O, Ind)
    colnames(DataM) <- c(paste0("y", 1:N), paste0("o", 1:M), "Ind")
  }else{
    DataM     <- data.frame(y, O, X, Ind)
    colnames(DataM) <- c(paste0("y", 1:N), paste0("o", 1:M), paste0("x", 1:K), "Ind")
  }
  irfest    <- array(NA, dim = c(HNum, N, E))
  irfse     <- array(NA, dim = c(HNum, N, E))
  IVRes     <- list()
  OLSRes    <- list()
  ORTHRes   <- list()

  # Compute lagged variables included in information set O(t-1)...O(t-P)
  infoVars  <- c()
  if(P > 0){
    for(j in 1:M){
      for (p in 1:P){
        if(interact == TRUE){
          # Compute lags of variables interacted with event dummies (excluding contaminated events)
          for(e in 0:1){
            DataM[, paste0("Event", e, ".o", j, ".l", p)] <- dplyr::lag(DataM[, paste0("o", j)], p) * (Ind == e)
            infoVars <- c(infoVars, paste0("Event", e, ".o", j, ".l", p))

            # Include interacted constant
            DataM[, paste0("Event", e)] <- (Ind == e)
            infoVars <- c(infoVars, paste0("Event", e))
          }
        }else{
          # Compute lags of variables
          DataM[, paste0("o", j, ".l", p)] <- dplyr::lag(DataM[, paste0("o", j)], p)
          infoVars <- c(infoVars, paste0("o", j, ".l", p))
        }
      }
    }
  }

  # Otherwise, only regress on a constant
  if(P == 0){
    if(interact == TRUE){
      for(e in 0:1){
        # Include interacted constant
        DataM[, paste0("Event", e)] <- (Ind == e)
        infoVars <- c(infoVars, paste0("Event", e))
      }
    }else{
      infoVars <- "1"
    }
  }

  # Add deterministic variables
  if(!is.null(X)){
    for(k in 1:K){
      DataM[, paste0("x", k)] <- DataM[, paste0("x", k)]
      infoVars <- c(infoVars, paste0("x", k))
    }
  }

  # Compute heteroskedasticity-based instrument separately for every dimension
  for(e in 1:E){

    # Set up the shock variable (instrumented variable), which is the e th variable in y
    DataM$shockVar <- DataM[, paste0("y", e)]

    # Set up the instrument variable, which is also based on the e th variable in y
    DataM$Ze <- DataM[,paste0("y", e)]

    # Set up the control variables for recursive ordering in case of E>1
    recVars <- c()
    recInst <- c()
    if(e>1){
      for(q in 1:(e-1)){
        recVars <- c(recVars, paste0("y", q))
        recInst <- c(recInst, paste0("Z", q))
      }
    }

    # Set up control variables for the local projections
    controls.info <- unique(c(infoVars))
    controls.info <- controls.info[controls.info !=""]

    controls.lp <- unique(c(infoVars, recVars))
    controls.lp <- controls.lp[controls.lp !=""]

    controls.iv <- unique(c(infoVars, recInst))
    controls.iv <- controls.iv[controls.iv !=""]

    # Account for missing data in instrumental variable by setting it to an other event
    # which is not used in the estimation
    DataM$Ind[is.na(DataM$Ze)] = 2

    # Set up event days (Policy event, control, and other event)
    DataM$Event    <- (DataM$Ind == 1)
    DataM$NoEvent  <- (DataM$Ind == 0)
    DataM$OthEvent <- (DataM$Ind == 2)

    # Shorten data to subset of observations without missing values at beginning or end
    DataMSub <- DataM[beg:end, ]

    # Orthogonalize the variable used to construct the instrument if we include lags of dependent variable
    # See Lewis (2022). Note that if no control variables, this is just a regression of Ze on a constant
    myFormula = paste0("Ze ~ ", paste(controls.info, collapse = "+"))
    orthModel <- lm(as.formula(myFormula), data = subset(DataMSub, Ind < 2), na.action="na.exclude")
    DataMSub$Ze <- residuals(orthModel, na.action="na.exclude")

    # Compute instrument and F-Statistic (see Lewis, 2022, and Rigobon and Sachs, 2004)
    Te <- sum(DataMSub$Event)
    Tn <- sum(DataMSub$NoEvent)
    To <- sum(DataMSub$OthEvent)
    Tt <- Te+Tn
    DataMSub$Z <- (DataMSub$Event*Tt/Te - DataMSub$NoEvent*Tt/Tn)*DataMSub$Ze

    # Save instrument in original data for later use
    DataM[beg:end, paste0("Z", e)] <- DataMSub$Z
    DataM[beg:end, "Z"] <- DataMSub$Z

    # Estimate the impulse responses based on heteroskedasticity-based IV for every variable in y
    for (i in 1:N)
    {
      # Set up dependent variable
      DataM$depVar  <- DataM[, paste0("y", i)]
      cumi          <- cum[i]

      for(h in HSeries){

        # Compute dependent variable at horizon h (level or cumulative response)
        if(cumi == T){
          for(f in 1:h){
            if(f == 1){
              DataM$depVar.h <- dplyr::lead(DataM$depVar, f-1)
            }else{
              DataM$depVar.h <- DataM$depVar.h + dplyr::lead(DataM$depVar, f-1)
            }
          }
        }else{
          DataM$depVar.h <- dplyr::lead(DataM$depVar, h-1)
        }

        # Shorten data to subset which contains no missing values
        DataMSub <- DataM[beg:end, ]

        # LP (Jorda, 2005): instrument shockVar with Z, control for information set
        myFormula <- paste0("depVar.h ~ shockVar + ",
                            paste(controls.lp, collapse = "+"),
                            paste0("| Z + ",
                                   paste(controls.iv, collapse = "+")))

        IV.mod  <- ivreg::ivreg(as.formula(myFormula), data = subset(DataMSub, Ind < 2))
        IV.se   <- sqrt(diag(sandwich::vcovHC(IV.mod, type = "HC0")))
        IV.sum  <- summary(IV.mod, vcov = sandwich::vcovHC(IV.mod, type = "HC0"), df = Inf, diagnostics = TRUE)

        
        if(details == TRUE){
          # Compute OLS residuals on event and control days (once per shock, at h=1)
          # Save residuals for later computation of variance-covariance matrix
          if (e == 1 & h == 1){
            DataMSub$depVar.h2 <- DataMSub$depVar.h
            DataMSub$depVar.h2[DataMSub$Ind == 2] <- NA

            myFormula <- paste0("depVar.h2 ~", paste(controls.info, collapse = "+"))
            OLS.mod  <- lm(as.formula(myFormula), data = subset(DataMSub, Ind < 2), na.action="na.exclude")

            eti  <- residuals(OLS.mod)
            eti[DataMSub$Event != 1] <- NA

            vti  <- residuals(OLS.mod)
            vti[DataMSub$NoEvent != 1] <- NA

            # Account for missing values when saving back to original data set
            DataM$eti <- NA
            DataM$eti[beg:end] <- eti
            eti <- DataM$eti

            DataM$vti <- NA
            DataM$vti[beg:end] <- vti
            vti <- DataM$vti

            if(i == 1){
              et <- data.frame(eti)
              vt <- data.frame(vti)
            }else{
              et <- data.frame(et, eti)
              vt <- data.frame(vt, vti)
            }

          }else{
            #eti  <- NA
            #vti  <- NA
            #OLS.mod <- NA
          }
        }

        # Normalizes IRFs to a specific value. Because initial response
        # normalized to unity, just multiply by NormFac
        #IRF[h, "irf"]   <- IV.mod$coefficients["shockVar"]
        #IRF[h, "se"]    <- IV.se["shockVar"]
        #irfest[h, i, e] <- IRF[h, "irf"]*NormFac
        #irfse[h, i, e]  <- IRF[h, "se"]*NormFac

        irfest[h, i, e] <- IV.mod$coefficients["shockVar"]*norm
        irfse[h, i, e]  <- IV.se["shockVar"]*norm

        # Save IV results for every variable and every horizon
        if(details == TRUE){
          IVRes[[paste0("IV.h", h, ".n", i, ".e", e)]] <- IV.mod

          if(h == 1){
            # Save OLS and orthogonalization results for horizon 1 only to save memory
            OLSRes[[paste0("OLS.h", h, ".n", i, ".e", e)]] <- OLS.mod
            ORTHRes[[paste0("ORTH.h", h, ".n", i, ".e", e)]] <- orthModel
          }
        }  
      }

      # Label rows of impulse responses to start at 0 (immediate response)
      dimnames(irfest)[[1]] <- HSeries-1
      dimnames(irfse)[[1]]  <- HSeries-1

    }
  }

  if(details == TRUE){
    # Compute variance-covariance matrix of residuals on event days, and impact matrix
    Sig  <- var(et, use = "complete.obs")

    if(sum(!is.na(vt)>0)){
      SigR <- var(vt, use = "complete.obs")
    }else{
      SigR <- NA
    }
    Psi  <- irfest[1 , ,]
  
    # Save data for weak instruments test by Lewis-Mertens (2025)
    if(controls.info[1] != "1"){
      WeakData <- data.frame(DataM[, paste0("y", 1:E)], DataM[, paste0("Z", 1:E)], DataM[, controls.info])
    }else{
      WeakData <- data.frame(DataM[, paste0("y", 1:E)], DataM[, paste0("Z", 1:E)])
    }
    if(E == 1){
      if(controls.info[1] != "1"){
        colnames(WeakData) <- c("y1", "Z1", controls.info)
      }else{
        colnames(WeakData) <- c("y1", "Z1")
      }
    }

    # Set data to missing if indicator is equal to 2
    WeakData[DataM$Ind == 2, ] <- NA

    Obs <- data.frame(Tp = Te, Tc = Tn, To = To, Tt = Tt)
  }
  
  Method = "Heteroscedasticity-IV"
  
  if(details == TRUE){

    return(list(irf = irfest, se = irfse,
              IVRes = IVRes, OLSRes = OLSRes, ORTHRes = ORTHRes,
              Obs = Obs, Method = Method,
              et = et, Sig = Sig, SigR = SigR, Psi = Psi, WeakData = WeakData))

  }else{
    return(list(irf = irfest, se = irfse, Method = Method))
  }
}
