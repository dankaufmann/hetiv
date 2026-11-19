# Internal helpers for the residuals used in shock prediction
#
# The estimators predict the structural shocks on event days from the
# forecast errors of the local projections at horizon `hpredict` (single
# horizon) or at all horizons 0, ..., hpredict jointly (pooled). These helpers
# compute the residuals at one horizon, flag control days whose window
# contains an event, and stack the residuals and responses across horizons
# in the layout expected by predict_shocks().

# TRUE for every t whose forward window t+1, ..., t+h contains an event or a
# contaminated day (Ind != 0). Used to drop such control days from the
# horizon-h control residuals, since their forecast error contains the
# response to the later event.
.window_has_event <- function(Ind, h) {
  n <- length(Ind)
  flag <- rep(FALSE, n)
  if (h < 1) return(flag)
  ev <- !is.na(Ind) & Ind != 0
  for (k in seq_len(h)) {
    lead_ev <- c(ev[(k + 1):n], rep(FALSE, k))[seq_len(n)]
    flag <- flag | lead_ev
  }
  flag
}

# Residuals of the information-set regression of depVar.h on the controls at
# horizon h (0-based), split into event-day and control-day residuals and
# padded to the full sample length Nobs.
#
# DataMSub must contain depVar.h, Ind, Event, NoEvent and the controls; Ind
# is the full-sample indicator (length Nobs); rows beg:end of the full sample
# correspond to DataMSub.
.shock_residuals <- function(DataMSub, controls.info, Ind, beg, end, Nobs, h) {
  # Exclude contaminated events (Ind == 2) from the regression
  DataMSub$depVar.h2 <- DataMSub$depVar.h
  DataMSub$depVar.h2[DataMSub$Ind == 2] <- NA

  # Formula in the base environment, so that a stored model does not keep the
  # calling function's frame alive (and is not saved with it)
  f <- stats::as.formula(
    paste0("depVar.h2 ~ ", paste(controls.info, collapse = "+")),
    env = baseenv()
  )
  mod <- stats::lm(f, data = DataMSub, na.action = "na.exclude")
  r <- stats::residuals(mod)

  eti <- r
  eti[DataMSub$Event != 1] <- NA

  vti <- r
  vti[DataMSub$NoEvent != 1] <- NA
  # A control day at horizon h must have no event inside its window
  if (h > 0) {
    vti[.window_has_event(Ind, h)[beg:end]] <- NA
  }

  et <- rep(NA_real_, Nobs)
  et[beg:end] <- eti
  vt <- rep(NA_real_, Nobs)
  vt[beg:end] <- vti

  list(et = et, vt = vt, mod = mod)
}

# Assemble the inputs of predict_shocks() from the residual lists and the IRF
# array. `hres` holds the (1-based) IRF rows used for prediction: a single
# row, or rows 1:(hpredict + 1) when pooled. Columns of et/vt and rows of Psi
# are ordered horizon-major: y1.h0, ..., yN.h0, y1.h1, ...
.stack_shock_inputs <- function(etList, vtList, irfest, hres, N, E) {
  labs <- as.vector(outer(paste0("y", seq_len(N)), paste0(".h", hres - 1), paste0))
  et <- as.data.frame(etList[labs])
  vt <- as.data.frame(vtList[labs])
  names(et) <- labs
  names(vt) <- labs

  Psi <- do.call(rbind, lapply(hres, function(hh) {
    matrix(irfest[hh, , , drop = FALSE], nrow = N, ncol = E)
  }))
  rownames(Psi) <- labs

  Sig <- stats::var(et, use = "complete.obs")
  SigR <- if (sum(stats::complete.cases(vt)) > 1) stats::var(vt, use = "complete.obs") else NA

  list(et = et, vt = vt, Psi = Psi, Sig = Sig, SigR = SigR)
}

# Implied variances of the E structural shocks, backed out from the difference
# of the residual covariance on event and control days under the unit-impact
# normalisation of Psi: Sig - SigR = sum_i sigma_i Psi_i Psi_i' (least squares
# on the vectorised matrices). Shared by predict_shocks() and
# prediction_corr(). Estimates can be negative under weak heteroskedasticity.
.shock_variances <- function(Sig, SigR, Psi) {
  q <- as.vector(Sig - SigR)
  A <- sapply(seq_len(ncol(Psi)), function(i) c(Psi[, i] %*% t(Psi[, i])))
  c(MASS::ginv(A) %*% q)
}

# Collect the coefficients of the horizon-0 information-set regressions (the
# regressions whose residuals are the one-step forecast errors) into arrays.
#
# coefs: list of length N with the named coefficient vector of the regression
#   of y_i on the P lags of the M variables in O (names o<j>.l<p>) and the K
#   exogenous variables (names x<k>), including the intercept.
# interact: TRUE if the lags are interacted with the event indicator
#   (hetiv(interact = TRUE); names o<j>.l<p>.i0 / .i1 and i1TRUE for the
#   event-day shift of the intercept).
#
# Returns Phi (N x M x P), Theta (N x K, NULL if K = 0) and const (length N).
# With interact = TRUE, Phi and const refer to control days, and PhiE and
# constE (event days) are returned in addition.
.info_coefficients <- function(coefs, N, M, P, K, interact = FALSE) {
  get <- function(cf, nm) if (nm %in% names(cf)) unname(cf[nm]) else NA_real_
  # rep() rather than paste0() alone so that P = 0 gives an empty dimension
  dn <- list(
    paste0("y", seq_len(N)), paste0("o", seq_len(M)),
    paste0(rep("l", P), seq_len(P))
  )
  Phi <- array(NA_real_, dim = c(N, M, P), dimnames = dn)
  PhiE <- Phi
  Theta <- if (K > 0) {
    matrix(NA_real_, N, K, dimnames = list(dn[[1]], paste0("x", seq_len(K))))
  } else {
    NULL
  }
  const <- stats::setNames(rep(NA_real_, N), dn[[1]])
  constE <- const

  for (i in seq_len(N)) {
    cf <- coefs[[i]]
    if (is.null(cf)) next
    const[i] <- get(cf, "(Intercept)")
    if (interact) constE[i] <- const[i] + get(cf, "i1TRUE")
    for (j in seq_len(M)) {
      for (p in seq_len(P)) {
        nm <- paste0("o", j, ".l", p)
        if (interact) {
          Phi[i, j, p] <- get(cf, paste0(nm, ".i0"))
          PhiE[i, j, p] <- get(cf, paste0(nm, ".i1"))
        } else {
          Phi[i, j, p] <- get(cf, nm)
        }
      }
    }
    for (k in seq_len(K)) Theta[i, k] <- get(cf, paste0("x", k))
  }

  out <- list(Phi = Phi, Theta = Theta, const = const)
  if (interact) {
    out$PhiE <- PhiE
    out$constE <- constE
  }
  out
}
