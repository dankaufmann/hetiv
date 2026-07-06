test_that("simulatedata preserves matrix dimensions and supports P = 0", {
  N <- 2
  R <- 1
  E <- 1
  P <- 0
  sim <- simulatedata(
    Phi = array(0, dim = c(N, N, P)),
    SigE = 2,
    PsiE = matrix(c(1, 0.5), N, E),
    PsiR = matrix(c(1, 0.2), N, R),
    Nobs = 30,
    Nbin = 5,
    N = N,
    R = R,
    E = E,
    Nevn = 5,
    P = P,
    eDist = 0,
    seed = 1
  )

  expect_equal(dim(sim$y), c(30, N))
  expect_equal(dim(sim$eR), c(30, R))
  expect_equal(dim(sim$eE), c(30, E))
})

test_that("computeirf validates malformed cum", {
  expect_error(
    computeirf(
      Psi = diag(2),
      Phi = array(0, dim = c(2, 2, 1)),
      H = 3,
      cum = c(TRUE, FALSE, TRUE)
    ),
    "cum must be"
  )
})

test_that("kfpredict accepts one-shock impact vectors", {
  out <- kfpredict(
    Sig = diag(2),
    SigR = diag(c(0.5, 0.5)),
    Psi = c(1, 0.5),
    et = matrix(rnorm(20), ncol = 2),
    scale = TRUE
  )

  expect_equal(dim(out), c(10, 1))
})

test_that("Hstep stores selected horizons by row position", {
  set.seed(2)
  Tobs <- 80
  y <- matrix(rnorm(Tobs * 2), Tobs, 2)
  Ind <- rep(0L, Tobs)
  Ind[seq(5, Tobs, by = 5)] <- 1L
  Z <- matrix(Ind * y[, 1] + rnorm(Tobs), Tobs, 1)

  het <- hetiv(y = y, O = y, Ind = Ind, P = 1, H = 5, Hstep = 2, details = TRUE)
  prox <- proxyiv(y = y, O = y, Z = Z, Ind = Ind, P = 1, H = 5, Hstep = 2, details = TRUE)

  expect_equal(dim(het$irf), c(3, 2, 1))
  expect_equal(dim(prox$irf), c(3, 2, 1))
  expect_equal(dim(het$Psi), c(2, 1))
  expect_equal(dim(prox$Psi), c(2, 1))
  expect_equal(dimnames(het$irf)[[1]], c("0", "2", "4"))
  expect_equal(dimnames(prox$irf)[[1]], c("0", "2", "4"))
})

test_that("negative normalization leaves standard errors non-negative", {
  set.seed(3)
  Tobs <- 80
  y <- matrix(rnorm(Tobs * 2), Tobs, 2)
  Ind <- rep(0L, Tobs)
  Ind[seq(5, Tobs, by = 5)] <- 1L
  Z <- matrix(Ind * y[, 1] + rnorm(Tobs), Tobs, 1)

  het <- hetiv(y = y, O = y, Ind = Ind, P = 1, H = 2, norm = -1)
  prox <- proxyiv(y = y, O = y, Z = Z, Ind = Ind, P = 1, H = 2, norm = -1)

  expect_true(all(het$se >= 0, na.rm = TRUE))
  expect_true(all(prox$se >= 0, na.rm = TRUE))
})

test_that("proxyiv marks missing control-day proxies as contaminated", {
  set.seed(4)
  Tobs <- 80
  y <- matrix(rnorm(Tobs * 2), Tobs, 2)
  Ind <- rep(0L, Tobs)
  Ind[seq(5, Tobs, by = 5)] <- 1L
  Z <- matrix(Ind * y[, 1] + rnorm(Tobs), Tobs, 1)
  Z[2:4, 1] <- NA
  Z[5, 1] <- NA

  prox <- proxyiv(y = y, O = y, Z = Z, Ind = Ind, P = 1, H = 2, details = TRUE)

  expect_equal(prox$Obs$To, 3)
  expect_lt(prox$Obs$Tiv, prox$Obs$Tt)
})
