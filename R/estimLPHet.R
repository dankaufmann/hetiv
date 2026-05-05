#' Estimate impulse responses via heteroscedasticity-based IV local projections
#'
#' Estimates impulse response functions (IRFs) using recursive
#' heteroscedasticity-IV identification (Rigobon, 2003; Lewis, 2022) combined
#' with local projections (Jordà, 2005). Identification exploits the difference
#' in variance between policy event days and control days to construct
#' instruments for the endogenous variables.
#'
#' @param y Numeric matrix of outcome variables (T x N). The first variable in
#'   each dimension is normalized to unity at horizon 0. These variables are
#'   also used as instruments and for recursive ordering.
#' @param O Numeric matrix of information set variables (T x M). May be
#'   identical to `y`. Included as lags 1 through `P`.
#' @param Ind Integer vector of length T, event indicator:
#'   \itemize{
#'     \item `0` Control day (no event)
#'     \item `1` Policy day (event)
#'     \item `2` Contaminated control day (excluded from estimation)
#'   }
#' @param Interact Logical. If `TRUE`, lagged information set variables are
#'   interacted with event/non-event dummies.
#' @param E Integer. Number of shock dimensions to identify via recursive
#'   heteroscedasticity ordering.
#' @param P Integer. Maximum lag order for the information set. Set to `0` for
#'   no lags (regression on constant only).
#' @param cumRep Logical vector of length N. For each variable in `y`, whether
#'   to report the cumulative impulse response instead of the level response.
#' @param H Integer. Maximum horizon (in periods) up to which IRFs are computed.
#' @param HStep Integer. Step size between horizons. If `> 1`, only every
#'   `HStep`-th horizon is estimated starting from `H = 1`.
#' @param NormFac Numeric scalar. Multiplies all IRFs to normalize the impact
#'   response of the first variable to a specific value. Set to `1` for no
#'   normalization.
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
#'       covariance estimation).}
#'     \item{`Sig`}{Covariance matrix of residuals on event days.}
#'     \item{`SigR`}{Covariance matrix of residuals on control days, or `NA`
#'       if unavailable.}
#'     \item{`Psi`}{Impact matrix (N x E), i.e. `irf[1, , ]`.}
#'     \item{`WeakData`}{Data frame of endogenous variables and instruments for
#'       the Lewis-Mertens (2024) weak instrument test.}
#'   }
#'
#' @references
#' Jordà, Ò. (2005). Estimation and Inference of Impulse Responses by Local
#' Projections. *American Economic Review*, 95(1), 161–182.
#'
#' Lewis, D. J. (2022). Identifying Shocks via Time-Varying Volatility.
#' *Review of Economic Studies*, 89(6), 2893–2937.
#'
#' Rigobon, R. (2003). Identification Through Heteroskedasticity.
#' *Review of Economics and Statistics*, 85(4), 777–792.
#'
#' @importFrom dplyr lag lead
#' @importFrom ivreg ivreg
#' @importFrom sandwich vcovHC
#'
#' @export
estimLPHet <- function(y, O, Ind, Interact, E, P, cumRep, H, HStep, NormFac){

  # Collect various properties of the data and observations to be used
  Nobs <- dim(y)[1]     # Number of observations in y
  beg  <- max(P+1)      # Start of the sample
  end  <- Nobs - H+1    # End of the sample
  N    <- dim(y)[2]     # Number of variables in y
  M    <- dim(O)[2]     # Number of variables in O

  if(sum(is.na(y))>0){
    warning("Missing values in y")
  }
  if(sum(is.na(O))>0){
    warning("Missing values in O")
  }
  if(sum(is.na(Ind))>0){
    warning("Missing values in Ind")
  }

  # Specify at which horizons IRF should be computed
  HSeries <- seq(1, H, HStep)
  HNum    <- length(HSeries)

  # Set up data set and and objects to save results
  DataM     <- data.frame(y, O, Ind)
  colnames(DataM) <- c(paste0("y", 1:N), paste0("o", 1:M), "Ind")
  irfest    <- array(NA, dim = c(HNum, N, E))
  irfse     <- array(NA, dim = c(HNum, N, E))
  IVRes     <- list()

  # Compute lagged variables included in information set O(t-1)
  infoVars      <- c()
  if(P > 0){
    for(j in 1:M){
      for (p in 1:P){
        if(Interact == TRUE){
          for(e in 0:1){
            DataM[, paste0("Event", e, ".o", j, ".l", p)] <- dplyr::lag(DataM[, paste0("o", j)], p) * (Ind == e)
            infoVars <- c(infoVars, paste0("Event", e, ".o", j, ".l", p))
          }
        }else{
          DataM[, paste0("o", j, ".l", p)] <- dplyr::lag(DataM[, paste0("o", j)], p)
          infoVars <- c(infoVars, paste0("o", j, ".l", p))
        }

      }
    }
  }
  if(P == 0){
    # Otherwise, only regress on a constant
    if(Interact == TRUE){
      for(e in 0:1){
        DataM[, paste0("Event", e)] <- (Ind == e)
        infoVars <- c(infoVars, paste0("Event", e))
      }
    }else{
      infoVars <- "1"
    }

  }

  # Compute instrument separately for every dimension
  for(e in 1:E){

    # Set up the shock variable (instrumented variable), which is the e th variable in y
    DataM$shockVar <- DataM[, paste0("y", e)]

    # Set up the instrument variable, which is also the e th variable in y
    DataM$Ze <- DataM[,paste0("y", e)]

    # Set up the control variables for recursive ordering
    # Note that we interact the dependent variable from previous
    # equation with the event dummy
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

    # Account for missing data in instrument variable by setting it to an other event
    # which is not used in the estimation
    DataM$Ind[is.na(DataM$Ze)] = 2

    # Set up event days (Policy event, control, and other event)
    DataM$Event    <- (DataM$Ind == 1)
    DataM$NoEvent  <- (DataM$Ind == 0)
    DataM$OthEvent <- (DataM$Ind == 2)

    # Shorten data to subset of observations without missing values at beginning or end
    DataMSub <- DataM[beg:end, ]

    # Orthogonalize the instrument if we include lags of dependent variable
    # See Lewis (2022)
    if(P > 0){
      myFormula = paste0("Ze ~ ", paste(controls.info, collapse = "+"))
      orthModel <- lm(as.formula(myFormula), data = subset(DataMSub, Ind < 2), na.action="na.exclude")
      DataMSub$Ze <- residuals(orthModel, na.action="na.exclude")
    }

    # Compute instrument and F-Statistic (see Lewis, 2022, ReStat, and
    # Rigobon, 2003, ReStat)
    Te <- sum(DataMSub$Event)
    Tn <- sum(DataMSub$NoEvent)
    To <- sum(DataMSub$OthEvent)
    Tt <- Te+Tn

    # Compute the instrument (see Lewis, 2022, ReStat and Rigobon and Sack, 2004, JME)
    DataMSub$Z <- (DataMSub$Event*Tt/Te - DataMSub$NoEvent*Tt/Tn)*DataMSub$Ze

    # Save instrument in original data for later use
    DataM[beg:end, paste0("Z", e)] <- DataMSub$Z
    DataM[beg:end, "Z"] <- DataMSub$Z

    # Estimate the impulse responses for every variable in y
    for (i in 1:N)
    {
      # Set up dependent variable
      DataM$depVar  <- DataM[, paste0("y", i)]
      cumRepi       <- cumRep[i]

      # Create object to save impulse responses and standard errors
      IRF <- array(data = NA, dim = c(HNum, 7))
      colnames(IRF) <- c("H", "irf", "se", "upper95", "lower95", "upper90", "lower90")
      rownames(IRF) <- HSeries

      for(h in HSeries){

        # Compute dependent variable at horizon h (level or cumulative response)
        if(cumRepi == T){
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

        # Compute OLS residuals on event and control days (once per shock, at h=1)
        if (e == 1 & h == 1){
          DataMSub$depVar.h2 <- DataMSub$depVar.h
          DataMSub$depVar.h2[DataMSub$Ind == 2] <- NA

          myFormula <- paste0("depVar.h2 ~", paste(controls.info, collapse = "+"))
          OLS.mod  <- lm(as.formula(myFormula), data = subset(DataMSub, Ind < 2), na.action="na.exclude")

          eti  <- residuals(OLS.mod)
          eti[DataMSub$Event != 1] <- NA

          vti  <- residuals(OLS.mod)
          vti[DataMSub$NoEvent != 1] <- NA
        }else{
          eti  <- NA
          vti  <- NA
          OLS.mod <- NA
        }

        # Account for missing values when saving back to original data set
        DataM$eti <- NA
        DataM$eti[beg:end] <- eti
        eti <- DataM$eti

        DataM$vti <- NA
        DataM$vti[beg:end] <- vti
        vti <- DataM$vti

        # Normalizes IRFs to a specific value. Because initial response
        # normalized to unity, just multiply by NormFac
        IRF[h, "irf"]   <- IV.mod$coefficients["shockVar"]
        IRF[h, "se"]    <- IV.se["shockVar"]
        irfest[h, i, e] <- IRF[h, "irf"]*NormFac
        irfse[h, i, e]  <- IRF[h, "se"]*NormFac

        # Save residuals for later computation of variance-covariance matrix
        if(h == 1 & e == 1){
          if(i == 1){
            et <- data.frame(eti)
            vt <- data.frame(vti)
          }else{
            et <- data.frame(et, eti)
            vt <- data.frame(vt, vti)
          }
        }

        # Save IV results for every variable and every horizon
        IVRes[[paste0("IV.h", h, ".n", i, ".e", e)]] <- IV.mod

      }

      # Label rows of impulse responses to start at 0 (immediate response)
      dimnames(irfest)[[1]] <- HSeries-1
      dimnames(irfse)[[1]]  <- HSeries-1

    }
  }

  # Compute variance-covariance matrix of residuals on event days, and impact matrix
  Sig  <- var(et, use = "complete.obs")

  if(sum(!is.na(vt)>0)){
    SigR <- var(vt, use = "complete.obs")
  }else{
    SigR <- NA
  }

  Psi  <- irfest[1 , ,]

  # Save data for weak instruments test by Lewis-Mertens (2024)
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

  Method = "Heteroscedasticity-IV"
  Obs <- data.frame(Tp = Te, Tc = Tn, To = To, Tt = Tt)

  return(list(irf = irfest, se = irfse,
              IVRes = IVRes,
              Obs = Obs, Method = Method,
              et = et, Sig = Sig, SigR = SigR, Psi = Psi, WeakData = WeakData))

}
