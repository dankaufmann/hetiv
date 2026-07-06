make_known_var_dgp <- function() {
  N <- 2
  P <- 1
  E <- 1
  R <- 2

  Phi <- array(0, dim = c(N, N, P))
  Phi[, , 1] <- matrix(
    c(0.35, 0.10,
      -0.05, 0.25),
    nrow = N,
    byrow = TRUE
  )

  PsiE <- matrix(c(1, 0.65), nrow = N, ncol = E)
  PsiR <- diag(c(1, 0.8), nrow = N)

  sim <- simulatedata(
    Phi = Phi,
    SigE = 6,
    PsiE = PsiE,
    PsiR = PsiR,
    Nobs = 8000,
    Nbin = 400,
    N = N,
    R = R,
    E = E,
    Nevn = 5,
    P = P,
    eDist = 0,
    seed = 42
  )

  list(
    sim = sim,
    Phi = Phi,
    PsiE = PsiE,
    P = P,
    H = 5,
    Ind = as.integer(sim$IndE[, 1])
  )
}

test_that("hetiv recovers known heteroskedastic-IV impulse responses", {
  dgp <- make_known_var_dgp()
  truth <- computeirf(Psi = dgp$PsiE, Phi = dgp$Phi, H = dgp$H, cum = FALSE)

  fit <- hetiv(
    y = dgp$sim$y,
    O = dgp$sim$y,
    Ind = dgp$Ind,
    P = dgp$P,
    H = dgp$H,
    details = TRUE
  )

  expect_equal(dim(fit$irf), dim(truth))
  expect_equal(as.numeric(fit$Psi), as.numeric(dgp$PsiE), tolerance = 0.04)
  expect_equal(as.numeric(fit$irf[, , 1]), as.numeric(truth[, , 1]),
               tolerance = 0.04)
})

test_that("proxyiv recovers known proxy-IV impulse responses", {
  dgp <- make_known_var_dgp()
  truth <- computeirf(Psi = dgp$PsiE, Phi = dgp$Phi, H = dgp$H, cum = FALSE)

  fit <- proxyiv(
    y = dgp$sim$y,
    O = dgp$sim$y,
    Z = dgp$sim$eE[, 1, drop = FALSE],
    Ind = dgp$Ind,
    P = dgp$P,
    H = dgp$H,
    details = TRUE
  )

  expect_equal(dim(fit$irf), dim(truth))
  expect_equal(as.numeric(fit$Psi), as.numeric(dgp$PsiE), tolerance = 0.04)
  expect_equal(as.numeric(fit$irf[, , 1]), as.numeric(truth[, , 1]),
               tolerance = 0.04)
})
