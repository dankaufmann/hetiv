# Tests for the coefficients of the one-step regressions (Phi, Theta, const)
# and for the correlation of the predicted with the true shocks (PredR2)

sim_var <- function(Nobs = 3000, SigE = 4, seed = 42) {
  N <- 3
  Phi <- array(0, dim = c(N, N, 2))
  Phi[, , 1] <- matrix(c(0.5, 0.1, 0, 0, 0.4, 0.1, 0.1, 0, 0.3), N, N)
  Phi[, , 2] <- 0.1 * diag(N)
  simulatedata(
    Phi = Phi, SigE = SigE, PsiE = matrix(c(1, 0.8, 0.5), N, 1), PsiR = diag(N),
    Nobs = Nobs, Nbin = 100, N = N, R = N, E = 1, Nevn = 5, P = 2, eDist = 0,
    seed = seed
  )
}

test_that("Phi, Theta and const match the one-step OLS regressions", {
  sim <- sim_var(Nobs = 400)
  y <- sim$y
  Ind <- as.vector(sim$IndE)
  set.seed(42)
  X <- matrix(rnorm(2 * nrow(y)), ncol = 2)
  Z <- y[, 1] * (Ind == 1) + rnorm(nrow(y), sd = 0.1) * (Ind == 1)

  # Reference: equation-by-equation OLS of y(t) on two lags of y and X
  lag1 <- rbind(NA, y[-nrow(y), ])
  lag2 <- rbind(NA, NA, y[-c(nrow(y) - 1, nrow(y)), ])
  ref <- lapply(1:3, function(i) stats::coef(stats::lm(y[, i] ~ lag1 + lag2 + X)))

  check <- function(r) {
    expect_equal(dim(r$Phi), c(3L, 3L, 2L))
    expect_equal(dim(r$Theta), c(3L, 2L))
    expect_length(r$const, 3)
    for (i in 1:3) {
      expect_equal(r$const[[i]], ref[[i]][[1]])
      expect_equal(unname(r$Phi[i, , 1]), unname(ref[[i]][2:4]))
      expect_equal(unname(r$Phi[i, , 2]), unname(ref[[i]][5:7]))
      expect_equal(unname(r$Theta[i, ]), unname(ref[[i]][8:9]))
    }
  }

  check(hetiv(y, O = y, X = X, Ind = Ind, P = 2, H = 3, details = TRUE))
  check(recols(y, O = y, X = X, Ind = Ind, P = 2, H = 3, details = TRUE))
  check(proxyiv(y, O = y, Z = Z, X = X, Ind = Ind, P = 2, H = 3, details = TRUE))

  # The coefficients always refer to horizon 0, whatever horizon is used for
  # the shock prediction
  r0 <- hetiv(y, O = y, X = X, Ind = Ind, P = 2, H = 3, details = TRUE)
  r2 <- hetiv(y, O = y, X = X, Ind = Ind, P = 2, H = 3, hpredict = 2, details = TRUE)
  expect_equal(r2$Phi, r0$Phi)
  expect_equal(r2$Theta, r0$Theta)
  expect_equal(r2$et, hetiv(y, O = y, X = X, Ind = Ind, P = 2, H = 3,
    hpredict = 2, details = TRUE)$et)

  # Without exogenous variables and lags
  r <- hetiv(y, O = y, Ind = Ind, P = 0, H = 2, details = TRUE)
  expect_null(r$Theta)
  expect_equal(dim(r$Phi), c(3L, 3L, 0L))
  expect_equal(unname(r$const), unname(colMeans(y)))
  expect_null(r$PhiE)
})

test_that("Phi recovers the VAR coefficients of the data generating process", {
  sim <- sim_var(Nobs = 20000)
  r <- hetiv(sim$y, O = sim$y, Ind = as.vector(sim$IndE), P = 2, H = 2, details = TRUE)
  expect_equal(unname(r$Phi), unname(sim$Phi), tolerance = 0.05)
})

test_that("hetiv returns regime-specific coefficients with interact = TRUE", {
  sim <- sim_var(Nobs = 600)
  y <- sim$y
  Ind <- as.vector(sim$IndE)
  r <- hetiv(y, O = y, Ind = Ind, P = 1, H = 2, interact = TRUE, details = TRUE)
  expect_equal(dim(r$Phi), c(3L, 3L, 1L))
  expect_equal(dim(r$PhiE), c(3L, 3L, 1L))

  # Reference: separate OLS regressions on control and event days
  lag1 <- rbind(NA, y[-nrow(y), ])
  for (i in 1:3) {
    c0 <- stats::coef(stats::lm(y[, i] ~ lag1, subset = Ind == 0))
    c1 <- stats::coef(stats::lm(y[, i] ~ lag1, subset = Ind == 1))
    expect_equal(unname(r$Phi[i, , 1]), unname(c0[2:4]))
    expect_equal(unname(r$PhiE[i, , 1]), unname(c1[2:4]))
    expect_equal(r$const[[i]], c0[[1]])
    expect_equal(r$constE[[i]], c1[[1]])
  }
})

test_that("prediction_corr gives the squared correlation with the true shocks", {
  # Known population value: one unit-variance shock, S = Psi' SigR^-1 Psi
  Psi <- matrix(c(1, 0.5), nrow = 2)
  SigR <- diag(c(0.5, 2))
  S <- 1 / 0.5 + 0.25 / 2
  expect_equal(prediction_corr(Psi, SigR), S / (1 + S))
  expect_equal(prediction_corr(Psi, SigR, squared = FALSE), sqrt(S / (1 + S)))

  # Unit-impact Psi and a shock with variance 4: rescaled through Sig
  Sig <- SigR + 4 * Psi %*% t(Psi)
  expect_equal(prediction_corr(Psi, SigR, Sig = Sig), 4 * S / (1 + 4 * S))
  # Without Sig, the unit-impact Psi ignores the variance of the shock
  expect_equal(prediction_corr(Psi, SigR), S / (1 + S))

  # Two shocks
  Psi2 <- cbind(c(1, 0.5, 0.2), c(0, 1, 0.5))
  out <- prediction_corr(Psi2, diag(3))
  expect_length(out, 2)
  expect_true(all(out > 0 & out < 1))

  expect_error(prediction_corr(Psi, diag(3)), "Psi")

  # Nearly singular noise covariance: measure close to one, with a warning
  SigR_sing <- matrix(c(1, 0.99999, 0.99999, 1), 2)
  expect_warning(r2 <- prediction_corr(Psi, SigR_sing), "nearly\\) singular")
  expect_gt(r2, 0.99)
  expect_error(prediction_corr(Psi, SigR, Sig = diag(3)), "same dimensions")
})

test_that("PredR2 is close to the realised squared correlation in large samples", {
  # Holds when the responses are consistently estimated (HET-IV and proxy-IV
  # in this design, but not recursive OLS, which is biased here)
  for (SigE in c(0.5, 8)) {
    sim <- sim_var(Nobs = 20000, SigE = SigE)
    Ind <- as.vector(sim$IndE)
    ev <- Ind == 1
    set.seed(42)
    Z <- (sim$eE[, 1] + rnorm(nrow(sim$y), sd = 0.5)) * ev
    for (r in list(
      hetiv(sim$y, O = sim$y, Ind = Ind, P = 2, H = 2, details = TRUE),
      proxyiv(sim$y, O = sim$y, Z = Z, Ind = Ind, P = 2, H = 2, details = TRUE)
    )) {
      realised <- stats::cor(r$Shocks[ev, 1], sim$eE[ev, 1], use = "complete.obs")
      expect_length(r$PredR2, 1)
      expect_equal(r$PredR2, realised^2, tolerance = 0.03)
    }
  }
})
