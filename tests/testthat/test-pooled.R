make_pool_data <- function(Tobs = 300, seed = 7) {
  Phi <- array(0.3 * diag(2), c(2, 2, 1))
  d <- simulatedata(
    Phi = Phi, SigE = 1, PsiE = matrix(c(1, 0.5), 2), PsiR = diag(2),
    Nobs = Tobs, Nbin = 50, N = 2, R = 2, E = 1, Nevn = 10, P = 1,
    eDist = 0, seed = seed
  )
  list(y = d$y, Ind = as.integer(d$IndE[, 1]), eE = d$eE)
}

test_that(".window_has_event flags control days with an event in the window", {
  Ind <- c(0, 0, 1, 0, 0, 2, 0, 0)
  expect_equal(.window_has_event(Ind, 0), rep(FALSE, 8))
  expect_equal(.window_has_event(Ind, 1), c(FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE))
  expect_equal(.window_has_event(Ind, 2), c(TRUE, TRUE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE))
})

test_that("pooled = TRUE with hpredict = 0 reproduces the single-horizon prediction", {
  d <- make_pool_data()
  for (f in list(hetiv, recols)) {
    a <- suppressWarnings(f(y = d$y, O = d$y, Ind = d$Ind, P = 1, H = 4, details = TRUE))
    b <- suppressWarnings(f(y = d$y, O = d$y, Ind = d$Ind, P = 1, H = 4, pooled = TRUE, details = TRUE))
    expect_equal(a$Shocks, b$Shocks)
    expect_equal(a$Psi, b$Psi)
    expect_equal(colnames(a$et), c("y1.h0", "y2.h0"))
  }
})

test_that("pooled prediction stacks residuals and responses over horizons", {
  d <- make_pool_data()
  Z <- d$eE + 0.5
  Z[d$Ind == 0] <- 0
  fits <- list(
    het = suppressWarnings(hetiv(y = d$y, O = d$y, Ind = d$Ind, P = 1, H = 4,
      hpredict = 2, pooled = TRUE, details = TRUE)),
    prox = suppressWarnings(proxyiv(y = d$y, O = d$y, Z = Z, Ind = d$Ind, P = 1, H = 4,
      hpredict = 2, pooled = TRUE, details = TRUE)),
    rec = suppressWarnings(recols(y = d$y, O = d$y, Ind = d$Ind, P = 1, H = 4,
      hpredict = 2, pooled = TRUE, details = TRUE))
  )
  for (r in fits) {
    expect_equal(dim(r$Psi), c(6, 1))
    expect_equal(colnames(r$et), c("y1.h0", "y2.h0", "y1.h1", "y2.h1", "y1.h2", "y2.h2"))
    expect_equal(dim(r$Sig), c(6, 6))
    expect_equal(dim(r$SigR), c(6, 6))
    # Psi stacks the IRF rows 1..3 horizon by horizon
    expect_equal(unname(r$Psi[, 1]), c(r$irf[1, , 1], r$irf[2, , 1], r$irf[3, , 1]))
    # Shocks reproduce kfpredict() on the stored inputs
    expect_equal(r$Shocks, kfpredict(r$Sig, r$SigR, r$Psi, r$et))
    # ... and on the 3-D array interface
    et3 <- array(as.matrix(r$et), dim = c(nrow(r$et), 2, 3))
    Psi3 <- array(r$irf[1:3, , 1], dim = c(3, 2, 1))
    Psi3 <- aperm(Psi3, c(2, 3, 1)) # N x E x M
    expect_equal(r$Shocks, kfpredict(r$Sig, r$SigR, Psi3, et3))
    # Residuals at horizon 0 equal the single-horizon ones
    expect_true(all(is.na(r$Shocks[d$Ind == 0, ])))
  }
  single <- suppressWarnings(hetiv(y = d$y, O = d$y, Ind = d$Ind, P = 1, H = 4, details = TRUE))
  expect_equal(fits$het$et[, 1:2], single$et)
  expect_equal(names(fits$het$OLSRes)[1:3], c("OLS.n1.h0", "OLS.n1.h1", "OLS.n1.h2"))
  # Pooled shocks are highly, but not perfectly, correlated with the impact ones
  rho <- cor(single$Shocks, fits$het$Shocks, use = "complete.obs")
  expect_gt(rho, 0.5)
  expect_lt(rho, 1)
})

test_that(".shock_residuals drops control days with an event in the window", {
  set.seed(1)
  n <- 40
  Ind <- rep(0, n)
  Ind[c(10, 20, 30)] <- 1
  Ind[25] <- 2
  x <- rnorm(n)
  Dat <- data.frame(depVar.h = x + rnorm(n), x1 = x, Ind = Ind,
    Event = Ind == 1, NoEvent = Ind == 0)
  beg <- 2
  end <- n - 1
  r0 <- .shock_residuals(Dat[beg:end, ], "x1", Ind, beg, end, n, h = 0)
  r2 <- .shock_residuals(Dat[beg:end, ], "x1", Ind, beg, end, n, h = 2)
  expect_length(r0$et, n)
  expect_length(r2$vt, n)
  # Event-day residuals are unaffected by the horizon
  expect_equal(is.na(r0$et), is.na(r2$et))
  expect_equal(which(!is.na(r0$et)), c(10, 20, 30))
  # Control residuals are missing on event, contaminated and padded rows ...
  expect_true(all(is.na(r0$vt[c(1, 10, 20, 25, 30, n)])))
  # ... and at horizon 2 also on the two days before every event
  extra <- setdiff(which(is.na(r2$vt)), which(is.na(r0$vt)))
  expect_equal(extra, c(8, 9, 18, 19, 23, 24, 28, 29))
})

test_that("kfpredict warns when the residual covariance is singular", {
  et <- matrix(rnorm(12), nrow = 3, ncol = 4)
  Sig <- var(et)
  expect_warning(
    kfpredict(Sig = Sig, SigR = Sig, Psi = matrix(1, 4, 1), et = et, scale = FALSE),
    "singular"
  )
  expect_error(
    hetiv(y = make_pool_data()$y, O = make_pool_data()$y, Ind = make_pool_data()$Ind,
      P = 1, H = 4, pooled = "yes"),
    "pooled"
  )
})
