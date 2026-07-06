test_that("simulatedata covers alternative shock distributions and errors", {
  N <- 2
  R <- 2
  E <- 1
  P <- 1
  Phi <- array(0, dim = c(N, N, P))
  Phi[, , 1] <- diag(c(0.2, 0.1))
  PsiE <- matrix(c(1, 0.4), N, E)
  PsiR <- diag(N)

  t_sim <- simulatedata(
    Phi = Phi, SigE = 3, PsiE = PsiE, PsiR = PsiR,
    Nobs = 25, Nbin = 5, N = N, R = R, E = E, Nevn = 0,
    P = P, eDist = 5, seed = NA_real_
  )
  expect_equal(sum(t_sim$IndE), 0)
  expect_equal(as.numeric(t_sim$eE), rep(0, 25))

  garch_sim <- simulatedata(
    Phi = Phi, SigE = 2, PsiE = PsiE, PsiR = PsiR,
    Nobs = 25, Nbin = 5, N = N, R = R, E = E, Nevn = 4,
    P = P, eDist = c(0.1, 0.3), seed = 11
  )
  expect_equal(dim(garch_sim$y), c(25, N))
  expect_true(any(garch_sim$IndE == 1))

  expect_error(
    simulatedata(Phi, 1, PsiE, PsiR, 10, 1, N, R, E, 5, P, 2, 1),
    "Student-t shocks require"
  )
  expect_error(
    simulatedata(Phi, 1, PsiE, PsiR, 10, 1, N, R, E, 5, P, -1, 1),
    "eDist must be"
  )
  expect_error(
    simulatedata(Phi, 1, PsiE, PsiR, 10, 1, N, R, E, 5, P,
                 c(0.8, 0.3), 1),
    "GARCH parameters"
  )
  expect_error(
    simulatedata(Phi, 1, cbind(PsiE, PsiE), PsiR, 10, 1, N, R, 2, 5, P,
                 c(0.1, 0.3), 1),
    "GARCH shocks"
  )
  expect_error(
    simulatedata(array(0, dim = c(N, N, 20)), 1, PsiE, PsiR, 10, 1,
                 N, R, E, 5, 20, 0, 1),
    "P must be smaller"
  )
  expect_error(
    simulatedata(Phi, 1, PsiE, PsiR, 10, 1, N, R, E, 5, P, NA, 1),
    "eDist must be"
  )
  expect_error(
    simulatedata(Phi, 1, PsiE, PsiR, 10, 1, N, R, E, 5, P, 0, "seed"),
    "seed must be"
  )
})

test_that("exported functions surface validation errors clearly", {
  y <- matrix(rnorm(40), ncol = 2)
  Ind <- rep(c(0L, 1L), length.out = nrow(y))

  expect_error(hetiv(y = "bad", O = y, Ind = Ind, P = 1, H = 2),
               "y must be numeric")
  expect_error(hetiv(y = y, O = y[-1, ], Ind = Ind, P = 1, H = 2),
               "O must have")
  expect_error(hetiv(y = y, O = y, X = y[-1, 1, drop = FALSE], Ind = Ind,
                     P = 1, H = 2),
               "X must have")
  expect_error(hetiv(y = y, O = y, Ind = Ind, P = 1.5, H = 2),
               "P must be an integer")
  expect_error(hetiv(y = y, O = y, Ind = Ind, P = 1, H = 2, E = 3),
               "E cannot exceed")
  expect_error(hetiv(y = y[1:3, ], O = y[1:3, ], Ind = Ind[1:3], P = 1, H = 2),
               "Not enough observations")
  expect_error(hetiv(y = y, O = y, Ind = Ind, P = 1, H = 2, cov_type = "bad"),
               "'arg' should be")
  expect_error(hetiv(y = y, O = y, Ind = Ind, P = 1, H = 2, details = NA),
               "details must be")

  expect_error(proxyiv(y = y, O = y, Z = matrix(rnorm(nrow(y)), ncol = 1),
                       Ind = Ind, P = 1, H = 2, E = 2),
               "Z must have at least")
  expect_error(proxyiv(y = y, O = y, Z = matrix(rnorm(nrow(y)), ncol = 1),
                       Ind = Ind, P = 1, H = 2, recursive = NA),
               "recursive must be")

  expect_error(computeirf(Psi = matrix(1, 2, 1), Phi = matrix(0, 2, 2),
                          H = 2, cum = FALSE),
               "Phi must be")
  expect_error(computeirf(Psi = matrix(1, 3, 1),
                          Phi = array(0, dim = c(2, 2, 1)), H = 2,
                          cum = FALSE),
               "Psi must have")
  expect_error(computeirf(Psi = matrix(1, 2, 1),
                          Phi = array(0, dim = c(2, 3, 1)), H = 2,
                          cum = FALSE),
               "square")

  expect_warning(logdiff(c(1, 0, 2)), "Non-positive")
  expect_error(filllinear("bad", gap = 1), "DF must be")
  expect_error(filllinear(data.frame(x = "bad"), gap = 1), "All columns")
  expect_error(filllinear(data.frame(x = c(1, NA, 3)), gap = 1.5),
               "gap must be")
  expect_error(normalize("bad"), "x must be numeric")
})

test_that("estimators cover optional controls, cumulation, interactions, and HAC covariance", {
  set.seed(22)
  y <- matrix(rnorm(120), ncol = 2)
  Ind <- rep(0L, nrow(y))
  Ind[seq(4, nrow(y), by = 4)] <- 1L
  X <- matrix(seq_len(nrow(y)) / nrow(y), ncol = 1)
  Z <- matrix(Ind * y[, 1] + rnorm(nrow(y)), ncol = 1)

  het <- hetiv(y = y, O = y, X = X, Ind = Ind, P = 0, H = 3,
               interact = TRUE, cum = c(TRUE, FALSE), cov_type = "NW")
  prox <- proxyiv(y = y, O = y, X = X, Z = Z, Ind = Ind, P = 0, H = 3,
                  cum = c(TRUE, FALSE), cov_type = "NW")

  expect_equal(dim(het$irf), c(3, 2, 1))
  expect_equal(dim(prox$irf), c(3, 2, 1))
  expect_true(all(is.finite(het$irf)))
  expect_true(all(is.finite(prox$irf)))
})

test_that("kfpredict covers unscaled and weak-heteroskedasticity branches", {
  Sig <- diag(2)
  Psi <- matrix(c(1, 0.5), nrow = 2)
  et <- matrix(rnorm(20), ncol = 2)

  raw <- kfpredict(Sig = Sig, SigR = diag(2), Psi = Psi, et = et,
                   scale = FALSE, tol = 0)
  expect_equal(dim(raw), c(10, 1))

  expect_warning(
    scaled <- kfpredict(Sig = diag(c(0.5, 0.5)), SigR = diag(2),
                        Psi = Psi, et = et, scale = TRUE),
    "negative"
  )
  expect_equal(dim(scaled), c(10, 1))
})

test_that("gweakivtest covers HAC, relative criterion, targeted beta, and helpers", {
  expect_error(hetiv:::.mat_pow(matrix(1, 2, 3), 0.5), "square")
  expect_error(hetiv:::.mat_pow(diag(c(1, -1)), 0.5), "positive semidefinite")
  expect_error(hetiv:::.mat_pow(matrix(c(1, 0, 0, 0), 2, 2), -0.5),
               "positive definite")

  spd <- hetiv:::.nearestSPD(matrix(c(1, 2, 2, 1), 2, 2))
  expect_true(all(eigen(spd, symmetric = TRUE, only.values = TRUE)$values > 0))

  q <- hetiv:::.myQR(matrix(c(1, 2, 3, 4, 5, 7), nrow = 3), 2)
  expect_equal(crossprod(q), diag(2), tolerance = 1e-10)

  opt <- hetiv:::.OptStiefelGBB(
    X = matrix(c(1, 0, 0, 1), nrow = 2),
    fun = function(x) list(fval = -sum(diag(x)), gradient = -diag(2)),
    opts = list(mxitr = 2)
  )
  expect_equal(crossprod(opt$X), diag(2), tolerance = 1e-8)

  set.seed(123)
  n <- 80
  Z <- cbind(rnorm(n), rnorm(n), rnorm(n))
  X <- matrix(numeric(0), n, 0)
  Y <- matrix(0.7 * Z[, 1] - 0.2 * Z[, 2] + 0.1 * Z[, 3] + rnorm(n), ncol = 1)
  y <- 1.4 * Y[, 1] + rnorm(n)
  y[1] <- NA

  out <- gweakivtest(y = y, Y = Y, X = X, Z = Z, cov_type = "NW",
                     crit = "rel", target = 1, points = 1, seed = 3)
  expect_equal(out$nobs, n - 1)
  expect_equal(out$target, "beta_1")
  expect_equal(out$criterion, "rel")
  expect_true(is.finite(out$gmin_generalized_critical_value))

  expect_error(gweakivtest(y, cbind(Y, Y), X, Z[, 1, drop = FALSE]),
               "Not identified")
  expect_error(gweakivtest(y, Y, X, Z, target = 2), "target cannot exceed")
  expect_error(gweakivtest(y, Y, X, Z, alfa = 1), "alfa must be")
  expect_error(gweakivtest(y, Y, X, Z, tau = 0), "tau must be")
  expect_error(gweakivtest(y, Y, X, Z, points = 0), "points must be")
  expect_error(gweakivtest(y, Y, X, Z, seed = "bad"), "seed must be")
  expect_error(gweakivtest(y, Y, X, Z, target = 0), "target must be")
})
