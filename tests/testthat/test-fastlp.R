# The estimators compute the local projections by matrix algebra (.lp_fit())
# unless the model objects are requested (details = TRUE, models = TRUE), in
# which case ivreg()/lm() and sandwich are used. Both paths must give the same
# impulse responses and standard errors.

make_lp_data <- function(E = 1, Nobs = 300, seed = 42) {
  N <- 4
  Phi <- array(0, dim = c(N, N, 2))
  Phi[, , 1] <- 0.4 * diag(N) + 0.05
  Phi[, , 2] <- 0.1 * diag(N)
  PsiE <- if (E == 1) {
    matrix(c(1, 0.8, 0.5, 0.2), N, 1)
  } else {
    cbind(c(1, 0.2, 0.5, 0.3), c(0, 1, 0.4, -0.3))
  }
  sim <- simulatedata(
    Phi = Phi, SigE = 2, PsiE = PsiE, PsiR = diag(N), Nobs = Nobs, Nbin = 50,
    N = N, R = N, E = E, Nevn = 5, P = 2, eDist = 0, seed = seed
  )
  Ind <- as.vector(sim$IndE)
  set.seed(seed)
  Z <- (sim$eE + matrix(rnorm(length(sim$eE)), ncol = E)) * (Ind == 1)
  X <- matrix(rnorm(2 * Nobs), ncol = 2)
  list(y = sim$y, Ind = Ind, Z = Z, X = X)
}

# Same call with the fast path and with the model path
both_paths <- function(f, ...) {
  fast <- suppressWarnings(f(..., details = TRUE, models = FALSE))
  slow <- suppressWarnings(f(..., details = TRUE, models = TRUE))
  list(fast = fast, slow = slow)
}

expect_same_irf <- function(r) {
  expect_equal(r$fast$irf, r$slow$irf, tolerance = 1e-8)
  expect_equal(r$fast$se, r$slow$se, tolerance = 1e-8)
}

test_that("fast and model paths agree for hetiv", {
  d <- make_lp_data()
  for (ct in c("HC3", "NW")) {
    expect_same_irf(both_paths(hetiv, y = d$y, O = d$y, Ind = d$Ind, P = 2, H = 4, cov_type = ct))
    expect_same_irf(both_paths(hetiv, y = d$y, O = d$y, X = d$X, Ind = d$Ind, P = 2, H = 4,
      cov_type = ct, cum = c(TRUE, FALSE, TRUE, FALSE)))
    expect_same_irf(both_paths(hetiv, y = d$y, O = d$y, Ind = d$Ind, P = 0, H = 3, cov_type = ct))
    expect_same_irf(both_paths(hetiv, y = d$y, O = d$y, Ind = d$Ind, P = 1, H = 3,
      cov_type = ct, interact = TRUE))
  }
  # Two shocks (recursive controls and instruments), other information set
  d2 <- make_lp_data(E = 2)
  for (ct in c("HC3", "NW")) {
    expect_same_irf(both_paths(hetiv, y = d2$y, O = d2$y[, 1:2], X = d2$X, Ind = d2$Ind,
      P = 2, H = 3, E = 2, cov_type = ct, norm = 0.25))
  }
})

test_that("fast and model paths agree for proxyiv", {
  d <- make_lp_data()
  d2 <- make_lp_data(E = 2)
  for (ct in c("HC3", "NW")) {
    expect_same_irf(both_paths(proxyiv, y = d$y, O = d$y, Z = d$Z, X = d$X, Ind = d$Ind,
      P = 2, H = 4, cov_type = ct, cum = TRUE))
    expect_same_irf(both_paths(proxyiv, y = d$y, O = d$y, Z = d$Z, Ind = d$Ind,
      P = 0, H = 3, cov_type = ct))
    for (rec in c(FALSE, TRUE)) {
      expect_same_irf(both_paths(proxyiv, y = d2$y, O = d2$y, Z = d2$Z, Ind = d2$Ind,
        P = 1, H = 3, E = 2, cov_type = ct, recursive = rec))
    }
  }
})

test_that("fast and model paths agree for recols", {
  d <- make_lp_data()
  d2 <- make_lp_data(E = 2)
  for (ct in c("HC3", "NW")) {
    expect_same_irf(both_paths(recols, y = d$y, O = d$y, X = d$X, Ind = d$Ind,
      P = 2, H = 4, cov_type = ct, cum = c(FALSE, TRUE, FALSE, FALSE)))
    expect_same_irf(both_paths(recols, y = d$y, O = d$y, Ind = d$Ind, P = 0, H = 3, cov_type = ct))
    expect_same_irf(both_paths(recols, y = d2$y, O = d2$y, Ind = d2$Ind,
      P = 1, H = 3, E = 2, cov_type = ct))
  }
})

test_that("fast and model paths agree with missing values and contaminated days", {
  d <- make_lp_data()
  y <- d$y
  y[c(40, 41, 120), 3] <- NA
  y[200, 1] <- NA
  Ind <- d$Ind
  Ind[c(15, 85, 150)] <- 2L
  Z <- d$Z
  Z[which(Ind == 1)[c(3, 10)], 1] <- NA
  for (ct in c("HC3", "NW")) {
    expect_same_irf(both_paths(hetiv, y = y, O = y, Ind = Ind, P = 2, H = 4, cov_type = ct))
    expect_same_irf(both_paths(proxyiv, y = y, O = y, Z = Z, Ind = Ind, P = 2, H = 4, cov_type = ct))
    expect_same_irf(both_paths(recols, y = y, O = y, Ind = Ind, P = 2, H = 4, cov_type = ct))
  }
})

test_that("collinear controls fall back to the model path", {
  d <- make_lp_data()
  X <- cbind(d$X, d$X[, 1] + d$X[, 2]) # third column is collinear
  expect_null(hetiv:::.lp_fit(d$y[, 2], d$y[, 1], W = X))
  expect_null(hetiv:::.lp_fit(d$y[, 2], d$y[, 1], W = d$X, z = d$Z[, 1], Wz = X))
  for (f in list(hetiv, recols)) {
    r <- both_paths(f, y = d$y, O = d$y, X = X, Ind = d$Ind, P = 1, H = 3)
    expect_same_irf(r)
    expect_false(anyNA(r$fast$irf))
  }
  r <- both_paths(proxyiv, y = d$y, O = d$y, Z = d$Z, X = X, Ind = d$Ind, P = 1, H = 3)
  expect_same_irf(r)
})

test_that(".lp_fit reproduces ivreg and lm with sandwich covariances", {
  d <- make_lp_data()
  dat <- data.frame(y = d$y[, 2], x = d$y[, 1], z = d$Z[, 1], w1 = d$X[, 1], w2 = d$X[, 2])
  dat$y[c(5, 50)] <- NA
  W <- as.matrix(dat[, c("w1", "w2")])

  iv <- ivreg::ivreg(y ~ x + w1 + w2 | z + w1 + w2, data = dat)
  ols <- stats::lm(y ~ x + w1 + w2, data = dat)
  for (ct in c("HC3", "NW")) {
    vc <- function(m) {
      if (ct == "NW") sandwich::NeweyWest(m, prewhite = FALSE, adjust = TRUE) else sandwich::vcovHC(m, type = "HC3")
    }
    f_iv <- hetiv:::.lp_fit(dat$y, dat$x, W = W, z = dat$z, cov_type = ct)
    expect_equal(f_iv$coef, unname(coef(iv)["x"]))
    expect_equal(f_iv$se, unname(sqrt(diag(vc(iv)))["x"]))
    f_ols <- hetiv:::.lp_fit(dat$y, dat$x, W = W, cov_type = ct)
    expect_equal(f_ols$coef, unname(coef(ols)["x"]))
    expect_equal(f_ols$se, unname(sqrt(diag(vc(ols)))["x"]))
  }
  # No controls: constant only
  f0 <- hetiv:::.lp_fit(dat$y, dat$x)
  expect_equal(f0$coef, unname(coef(stats::lm(y ~ x, data = dat))["x"]))
})

test_that("model objects are returned only with models = TRUE", {
  d <- make_lp_data()
  r <- both_paths(hetiv, y = d$y, O = d$y, Ind = d$Ind, P = 1, H = 3)
  expect_false(any(c("IVRes", "OLSRes", "ORTHRes") %in% names(r$fast)))
  expect_true(all(c("IVRes", "OLSRes", "ORTHRes") %in% names(r$slow)))
  expect_length(r$slow$IVRes, 3 * 4)
  # Everything else is identical
  shared <- setdiff(names(r$slow), c("IVRes", "OLSRes", "ORTHRes"))
  expect_equal(r$fast[shared], r$slow[shared], tolerance = 1e-8)

  p <- both_paths(proxyiv, y = d$y, O = d$y, Z = d$Z, Ind = d$Ind, P = 1, H = 3)
  expect_false(any(c("IVRes", "OLSRes") %in% names(p$fast)))
  expect_true(all(c("IVRes", "OLSRes") %in% names(p$slow)))
  o <- both_paths(recols, y = d$y, O = d$y, Ind = d$Ind, P = 1, H = 3)
  expect_false("OLSRes" %in% names(o$fast))
  expect_true("OLSRes" %in% names(o$slow))

  # models has no effect without details
  expect_named(hetiv(y = d$y, O = d$y, Ind = d$Ind, P = 1, H = 3, models = TRUE),
    c("irf", "se", "Method"))

  # Stored models do not carry the function environment: the saved result
  # stays small, and the models remain usable
  size_fast <- length(serialize(r$fast, NULL))
  size_slow <- length(serialize(r$slow, NULL))
  expect_lt(size_fast, 2e5)
  expect_lt(size_slow, 12 * 3e5)
  m <- r$slow$IVRes[[1]]
  expect_identical(environment(m$formula), baseenv())
  expect_equal(dim(stats::vcov(m)), c(6L, 6L))
  expect_length(stats::predict(m), stats::nobs(m))
})
