#' Simulate VAR data with heteroskedastic shocks and given parameters
#'
#' Simulates data from a VAR(P) model with two types of structural shocks:
#' regular shocks (always present) and event shocks (occurring every `Nevn`
#' periods). The shock distribution can be normal, Student-t, or GARCH(1,1).
#'
#' @param Phi Array of VAR coefficient matrices, dimension `N x N x P`.
#' @param SigE Variance of the event shocks. This applies to all E shocks.
#' @param PsiE Impact matrix for event shocks, dimension `N x E`.
#' @param PsiR Impact matrix for regular shocks, dimension `N x R`.
#' @param Nobs Number of observations to retain after burn-in. The returned
#'   data have exactly `Nobs` rows.
#' @param Nbin Number of burn-in observations. These are simulated to
#'   initialise the VAR but are discarded before returning.
#' @param N Number of variables in the VAR.
#' @param R Number of regular shocks.
#' @param E Number of event shocks. For GARCH shock distributions (`eDist = c(alpha, beta)`), only `E = 1` is supported.
#' @param Nevn Event frequency: an event shock occurs every `Nevn` periods.
#'   Set to `0` to suppress event shocks entirely (no heteroskedasticity).
#' @param P VAR lag order.
#' @param eDist Shock distribution. Use `0` for standard normal; a positive
#'   integer for Student-t with that many degrees of freedom; or a numeric
#'   vector `c(alpha, beta)` for GARCH(1,1) with ARCH effect `alpha` and
#'   persistence `beta`.
#' @param seed Integer seed passed to [set.seed()] for reproducibility.
#'   Use `NA` to skip seeding.
#'
#' @return A named list with components:
#' \describe{
#'   \item{y}{Simulated VAR data, dimension `Nobs x N` (burn-in discarded).}
#'   \item{IndE}{Event indicator matrix, dimension `Nobs x 1`.}
#'   \item{eR}{Simulated regular shocks, dimension `Nobs x R`.}
#'   \item{eE}{Simulated event shocks, dimension `Nobs x E`.}
#'   \item{e}{Composite structural shocks, dimension `Nobs x N`.}
#'   \item{Phi}{VAR coefficient array (returned unchanged).}
#'   \item{PsiE}{Event shock impact matrix (returned unchanged).}
#'   \item{PsiR}{Regular shock impact matrix (returned unchanged).}
#'   \item{SigE}{Event shock variance (returned unchanged).}
#' }
#'
#' @export
simulatedata <- function(Phi, SigE, PsiE, PsiR, Nobs, Nbin, N, R, E, Nevn, P, eDist, seed){
  # Function to simulate VAR(P) data with heteroskedastic shocks
  # Note that this only works for E = 1. For E > 1, GARCH shocks are not independent
  # and SigE is the scale of the shocks for all shocks
  
  # Initialize random number generator to get always the same simulations
  if(!is.na(seed)){
    set.seed(seed)
  }

  # Simulate shock series (iid N(0, 1))
  if(length(eDist) == 1){
    if(eDist == 0){
      eR <- matrix(rnorm((Nobs+Nbin)*R, mean=0, sd=1), (Nobs+Nbin), R)
      eE <- matrix(rnorm((Nobs+Nbin)*E, mean=0, sd=sqrt(SigE)), (Nobs+Nbin), E)
    }
    if(eDist > 0){

      # Scale required to get a variance of SigE and 1 for the two shocks
      scaleE <- sqrt(SigE * (eDist - 2) / eDist)
      scaleR <- sqrt(1 * (eDist - 2) / eDist)

      eR <- matrix(scaleR*rt((Nobs+Nbin)*R, df = eDist), (Nobs+Nbin), R)
      eE <- matrix(scaleE*rt((Nobs+Nbin)*E, df = eDist), (Nobs+Nbin), E)

    }
  }else{
    # GARCH(1,1) shocks: eDist = c(alpha, beta)
    if (E > 1) stop("GARCH shocks (eDist = c(alpha, beta)) are only supported for E = 1.")
    omega <- 0.1          # constant term
    alpha <- eDist[1]     # ARCH effect (impact of past shocks)
    beta  <- eDist[2]     # GARCH effect (persistence of volatility)

    # Initialize
    hE  <- numeric((Nobs+Nbin)*E)      # conditional variance
    eE  <- numeric((Nobs+Nbin)*E)      # shocks

    hR  <- numeric((Nobs+Nbin)*R)      # conditional variance
    eR  <- numeric((Nobs+Nbin)*R)      # shocks

    # Set initial variance
    hE[1] <- omega / (1 - alpha - beta)  # unconditional variance
    eE[1] <- rnorm(1, mean = 0, sd = sqrt(hE[1]))

    hR[1] <- omega / (1 - alpha - beta)  # unconditional variance
    eR[1] <- rnorm(1, mean = 0, sd = sqrt(hR[1]))

    # Generate the series
    for(t in 2:((Nobs+Nbin)*E)) {
      hE[t] <- omega + alpha * eE[t-1]^2 + beta * hE[t-1]
      eE[t] <- rnorm(1, mean = 0, sd = sqrt(hE[t]))
    }
    for(t in 2:((Nobs+Nbin)*R)) {
      hR[t] <- omega + alpha * eR[t-1]^2 + beta * hR[t-1]
      eR[t] <- rnorm(1, mean = 0, sd = sqrt(hR[t]))
    }

    eR <- matrix(eR, (Nobs+Nbin), R)
    eE <- matrix(eE, (Nobs+Nbin), E)
    eE <- eE * sqrt(SigE)
  }


  # Simulate data with VAR(P)
  # y(t) = Phi*y(t-1)+...+PsiR*eR         (if no event)
  # y(t) = Phi*y(t-1)+...+PsiR*eR+PsiE*eE (if event)
  y    <- array(data = 0, dim = c(Nobs+Nbin, N))
  e    <- array(data = 0, dim = c(Nobs+Nbin, N))
  IndE <- array(data = 0, dim = c(Nobs+Nbin, 1))
  for (t in (P+1):(Nobs+Nbin)){

    # Accumulate lags (VAR prediction)
    for (p in 1:P){
      y[t, ]  <- y[t, ] + Phi[, , p] %*% y[t-p, ]
    }

    # Add shocks for periods with an event (if set to 0, then no heteroskedasticity)
    if(Nevn != 0 && t %% Nevn == 0){
      e[t, ]  <-  PsiR %*% eR[t, ] + PsiE %*% eE[t, ]
      y[t, ]  <-  y[t, ] + e[t, ]
      IndE[t] <- 1

    }
    # Periods without an event
    else{
      e[t, ]  <-  PsiR %*% eR[t, ]
      y[t, ] <- y[t, ] + e[t, ]
    }

  }

  # Discard burn-in rows and return post-burn-in observations only
  keep <- (Nbin + 1):(Nobs + Nbin)
  return(list(y = y[keep, ], IndE = IndE[keep, , drop = FALSE],
              eR = eR[keep, ], eE = eE[keep, ], e = e[keep, ],
              Phi = Phi, PsiE = PsiE, PsiR = PsiR, SigE = SigE))

}
