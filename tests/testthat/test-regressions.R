make_lp_data <- function(Tobs = 90) {
  set.seed(101)
  y <- matrix(rnorm(Tobs * 2), Tobs, 2)
  Ind <- rep(0L, Tobs)
  Ind[seq(5, Tobs, by = 5)] <- 1L
  Z <- matrix(Ind * y[, 1] + rnorm(Tobs), Tobs, 1)
  list(y = y, O = y, Ind = Ind, Z = Z)
}

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

test_that("simulatedata rejects malformed inputs", {
  expect_error(
    simulatedata(
      Phi = array(0, dim = c(2, 2, 1)),
      SigE = 1,
      PsiE = matrix(1, 2, 1),
      PsiR = matrix(1, 2, 1),
      Nobs = 10,
      Nbin = 1,
      N = 2,
      R = 1,
      E = 2,
      Nevn = 5,
      P = 1,
      eDist = 0,
      seed = 1
    ),
    "PsiE must have"
  )
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

test_that("computeirf computes cumulated responses", {
  Psi <- matrix(c(1, 2), nrow = 2)
  Phi <- array(0, dim = c(2, 2, 1))
  Phi[, , 1] <- diag(c(0.5, 0.25))

  out <- computeirf(Psi, Phi, H = 3, cum = c(TRUE, FALSE))

  expect_equal(dim(out), c(3, 2, 1))
  expect_equal(as.numeric(out[, 1, 1]), c(1, 1.5, 1.75))
  expect_equal(as.numeric(out[, 2, 1]), c(2, 0.5, 0.125))
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

test_that("kfpredict validates dimensions", {
  expect_error(
    kfpredict(
      Sig = diag(2),
      SigR = diag(2),
      Psi = matrix(1, 3, 1),
      et = matrix(rnorm(20), ncol = 2)
    ),
    "Psi must have 2 rows"
  )
})

test_that("Hstep stores selected horizons by row position", {
  dat <- make_lp_data(80)

  het <- hetiv(
    y = dat$y, O = dat$O, Ind = dat$Ind, P = 1, H = 5,
    Hstep = 2, details = TRUE
  )
  prox <- proxyiv(
    y = dat$y, O = dat$O, Z = dat$Z, Ind = dat$Ind, P = 1,
    H = 5, Hstep = 2, details = TRUE
  )

  expect_equal(dim(het$irf), c(3, 2, 1))
  expect_equal(dim(prox$irf), c(3, 2, 1))
  expect_equal(dim(het$Psi), c(2, 1))
  expect_equal(dim(prox$Psi), c(2, 1))
  expect_equal(dimnames(het$irf)[[1]], c("0", "2", "4"))
  expect_equal(dimnames(prox$irf)[[1]], c("0", "2", "4"))
})

test_that("negative normalization leaves standard errors non-negative", {
  dat <- make_lp_data(80)

  het <- hetiv(y = dat$y, O = dat$O, Ind = dat$Ind, P = 1, H = 2, norm = -1)
  prox <- proxyiv(
    y = dat$y, O = dat$O, Z = dat$Z, Ind = dat$Ind, P = 1,
    H = 2, norm = -1
  )

  expect_true(all(het$se >= 0, na.rm = TRUE))
  expect_true(all(prox$se >= 0, na.rm = TRUE))
})

test_that("proxyiv marks missing control-day proxies as contaminated", {
  dat <- make_lp_data(80)
  y <- dat$y
  Ind <- dat$Ind
  Z <- dat$Z
  Z[2:4, 1] <- NA
  Z[5, 1] <- NA

  prox <- proxyiv(y = y, O = y, Z = Z, Ind = Ind, P = 1, H = 2, details = TRUE)

  expect_equal(prox$Obs$To, 3)
  expect_lt(prox$Obs$Tiv, prox$Obs$Tt)
})

test_that("hetiv and proxyiv reject invalid indicators and missing groups", {
  dat <- make_lp_data(40)

  expect_error(
    hetiv(y = dat$y, O = dat$O, Ind = rep(1L, 40), P = 1, H = 2),
    "at least one event day"
  )
  bad_ind <- dat$Ind
  bad_ind[1] <- 9L
  expect_error(
    proxyiv(y = dat$y, O = dat$O, Z = dat$Z, Ind = bad_ind, P = 1, H = 2),
    "Ind must contain only"
  )
})

test_that("transformations validate inputs and transform values", {
  expect_equal(as.numeric(firstdiff(c(1, 3, 6))), c(NA, 2, 3))
  expect_equal(as.numeric(logdiff(c(1, exp(1), exp(3)))), c(NA, 1, 2))
  expect_equal(filllinear(data.frame(x = c(1, NA, 3)), gap = 1)$x, c(1, 2, 3))
  expect_equal(as.numeric(normalize(c(1, 2, 3))), c(-1, 0, 1))
  expect_error(normalize(c(1, 1, 1)), "non-constant")
  expect_error(firstdiff("x"), "TS must be numeric")
})

test_that("plotting functions return one ggplot per variable-shock pair", {
  irf <- array(seq_len(12) / 10, dim = c(3, 2, 2))
  dimnames(irf)[[1]] <- 0:2
  se <- array(0.1, dim = dim(irf), dimnames = dimnames(irf))
  labels <- c("a", "b")

  p1 <- plotirf(irf, se, HTick = 1, Labels = labels)
  p2 <- plot2irf(irf, se, irf * 0.8, se, HTick = 1, Labels = labels)
  p3 <- plotpval(array(0.05, dim = dim(irf), dimnames = dimnames(irf)),
    HTick = 1, Labels = labels
  )

  expect_length(p1, 4)
  expect_length(p2, 4)
  expect_length(p3, 4)
  expect_s3_class(p1[[1]], "ggplot")
  expect_error(plotirf(irf, se, HTick = 1, Labels = "only one"), "Labels")
})

test_that("gweakivtest is deterministic with default seed and preserves RNG state", {
  set.seed(5)
  n <- 100
  Z <- cbind(rnorm(n), rnorm(n), rnorm(n))
  X <- matrix(rnorm(n), ncol = 1)
  Y <- cbind(
    0.8 * Z[, 1] + 0.3 * Z[, 2] + rnorm(n),
    0.7 * Z[, 2] + 0.2 * Z[, 3] + rnorm(n)
  )
  y <- 1.2 * Y[, 1] - 0.5 * Y[, 2] + X[, 1] + rnorm(n)

  set.seed(99)
  seed_before <- .Random.seed
  out1 <- gweakivtest(y, Y, X, Z, points = 2)
  seed_after <- .Random.seed
  out2 <- gweakivtest(y, Y, X, Z, points = 2)

  expect_equal(seed_after, seed_before)
  expect_equal(
    out1$gmin_generalized_critical_value,
    out2$gmin_generalized_critical_value
  )
  expect_true(is.finite(out1$gmin_generalized))
})

test_that("gweakivtest reports singular designs clearly", {
  n <- 30
  z <- rnorm(n)
  expect_error(
    gweakivtest(
      y = rnorm(n),
      Y = matrix(rnorm(n), ncol = 1),
      X = matrix(1, n, 1),
      Z = cbind(z, z)
    ),
    "Z is rank deficient"
  )
})
