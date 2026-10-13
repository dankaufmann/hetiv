# Internal helpers for the residuals used in shock prediction
#
# The estimators predict the structural shocks on event days from the
# forecast errors of the local projections at horizon `hpredict` (single
# horizon) or at all horizons 0, ..., hpredict jointly (pooled). These helpers
# compute the residuals at one horizon, flag control days whose window
# contains an event, and stack the residuals and responses across horizons
# in the layout expected by kfpredict().

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

  f <- stats::as.formula(paste0("depVar.h2 ~ ", paste(controls.info, collapse = "+")))
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

# Assemble the inputs of kfpredict() from the residual lists and the IRF
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
