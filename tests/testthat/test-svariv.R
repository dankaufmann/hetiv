make_svar_dgp <- function(Nobs = 4000, seed = 42) {
  N <- 2
  P <- 1
  E <- 1
  R <- 2

  Phi <- array(0, dim = c(N, N, P))
  Phi[, , 1] <- matrix(c(0.8, -0.10, 0.15, 0.5), nrow = N, byrow = TRUE)
  PsiE <- matrix(c(1, -0.65), nrow = N, ncol = E)
  PsiR <- diag(c(1, 0.8), nrow = N)

  sim <- simulatedata(
    Phi = Phi, SigE = 6, PsiE = PsiE, PsiR = PsiR,
    Nobs = Nobs, Nbin = 400, N = N, R = R, E = E, Nevn = 2, P = P,
    eDist = 0, seed = seed
  )
  set.seed(seed)
  Z <- sim$eE + rnorm(length(sim$eE), 0, 5) # noisy proxy

  list(y = sim$y, Z = Z, Phi = Phi, PsiE = PsiE, P = P, H = 10)
}

# Benchmark implementations, written independently of svariv() internals
var_ols <- function(y, P) {
  Nobs <- nrow(y)
  N <- ncol(y)
  W <- 1
  for (p in seq_len(P)) W <- cbind(W, y[(P + 1 - p):(Nobs - p), , drop = FALSE])
  Ydep <- y[(P + 1):Nobs, , drop = FALSE]
  coef <- solve(crossprod(W), crossprod(W, Ydep))
  list(coef = coef, U = Ydep - W %*% coef)
}

test_that("svariv reproduces OLS VAR estimates and the proxy identification", {
  dgp <- make_svar_dgp()
  fit <- svariv(y = dgp$y, Z = dgp$Z, P = 1, H = dgp$H, details = TRUE)

  bm <- var_ols(dgp$y, 1)
  expect_equal(unname(fit$C), unname(bm$coef[1, ]), tolerance = 1e-10)
  expect_equal(fit$Phi[, , 1], t(bm$coef[2:3, ]), tolerance = 1e-10, check.attributes = FALSE)
  expect_equal(fit$U, bm$U, tolerance = 1e-10, check.attributes = FALSE)

  # Single-instrument closed form: b21/b11 = cov(u2, z) / cov(u1, z)
  U <- bm$U
  z <- dgp$Z[2:nrow(dgp$Z), 1]
  ratio <- cov(U[, 2], z) / cov(U[, 1], z)
  expect_equal(fit$Psi[2, 1] / fit$Psi[1, 1], ratio, tolerance = 1e-10)
  expect_equal(fit$Psi[1, 1], 1)

  # First-stage statistics agree with lm()
  fs <- summary(lm(U[, 1] ~ z))
  expect_equal(unname(fit$FStat), unname(fs$fstatistic[1]), tolerance = 1e-8)
  expect_equal(unname(fit$R2), fs$r.squared, tolerance = 1e-10)

  # IRFs follow from computeirf() with the estimated VAR and impact matrix
  irf_bm <- computeirf(Psi = fit$Psi, Phi = fit$Phi, H = dgp$H, cum = FALSE)
  expect_equal(fit$irf, irf_bm)
  expect_equal(dimnames(fit$irf)[[1]], as.character(0:(dgp$H - 1)))
})

test_that("svariv recovers known impulse responses", {
  dgp <- make_svar_dgp()
  truth <- computeirf(Psi = dgp$PsiE, Phi = dgp$Phi, H = dgp$H, cum = FALSE)
  fit <- svariv(y = dgp$y, Z = dgp$Z, P = dgp$P, H = dgp$H, details = TRUE)

  expect_equal(dim(fit$irf), dim(truth))
  expect_equal(as.numeric(fit$Psi), as.numeric(dgp$PsiE), tolerance = 0.05)
  expect_lt(max(abs(fit$irf - truth)), 0.05)
  expect_equal(fit$Method, "SVAR-IV")
  expect_null(fit$se)
  expect_null(fit$lower)
  expect_equal(fit$Obs$Tiv, nrow(dgp$y) - dgp$P)

  # Cumulative responses and normalization
  fitc <- svariv(y = dgp$y, Z = dgp$Z, P = dgp$P, H = dgp$H, cum = TRUE, norm = -2)
  expect_equal(fitc$irf[, , 1], -2 * apply(fit$irf[, , 1], 2, cumsum),
    check.attributes = FALSE
  )
  expect_equal(unname(fitc$irf[1, 1, 1]), -2)
})

test_that("svariv drops missing proxies from identification only", {
  dgp <- make_svar_dgp()
  Zna <- dgp$Z
  Zna[seq(3, nrow(Zna), by = 3), ] <- NA
  fit <- svariv(y = dgp$y, Z = Zna, P = 1, H = 4, details = TRUE)
  full <- svariv(y = dgp$y, Z = dgp$Z, P = 1, H = 4, details = TRUE)

  expect_equal(fit$Obs$Tiv, sum(complete.cases(Zna[-1, , drop = FALSE])))
  expect_equal(fit$Phi, full$Phi) # VAR uses the full sample
  z <- Zna[-1, 1]
  keep <- !is.na(z)
  ratio <- cov(fit$U[keep, 2], z[keep]) / cov(fit$U[keep, 1], z[keep])
  expect_equal(fit$Psi[2, 1], ratio, tolerance = 1e-10)
})

test_that("svariv wild bootstrap matches a manual implementation", {
  dgp <- make_svar_dgp(Nobs = 600)
  H <- 5
  nboot <- 20
  set.seed(5)
  fit <- svariv(
    y = dgp$y, Z = dgp$Z, P = 1, H = H, nboot = nboot, ci = c(0.68, 0.90),
    details = TRUE
  )

  # Manual wild bootstrap (Mertens and Ravn, 2013) with the same RNG stream
  set.seed(5)
  Nobs <- nrow(dgp$y)
  N <- ncol(dgp$y)
  bm <- var_ols(dgp$y, 1)
  C <- bm$coef[1, ]
  B <- t(bm$coef[2:3, ])
  manual <- array(NA_real_, dim = c(H, N, nboot))
  for (j in seq_len(nboot)) {
    rr <- 1 - 2 * (runif(Nobs - 1) > 0.5)
    bootU <- bm$U * rr
    bootZ <- dgp$Z[-1, , drop = FALSE] * rr
    yb <- dgp$y
    for (t in 2:Nobs) yb[t, ] <- C + B %*% yb[t - 1, ] + bootU[t - 1, ]
    bb <- var_ols(yb, 1)
    z <- bootZ[, 1]
    ratio <- cov(bb$U[, 2], z) / cov(bb$U[, 1], z)
    Phib <- array(t(bb$coef[2:3, ]), dim = c(N, N, 1))
    manual[, , j] <- computeirf(matrix(c(1, ratio), N, 1), Phib, H, cum = FALSE)[, , 1]
  }

  expect_equal(unname(fit$boot[, , 1, ]), manual, tolerance = 1e-8)
  expect_equal(fit$nfail, 0L)
  expect_equal(dim(fit$boot), c(H, N, 1, nboot))
  expect_equal(dim(fit$se), c(H, N, 1))
  expect_equal(dim(fit$lower), c(H, N, 1, 2))
  expect_equal(dimnames(fit$lower)[[4]], c("0.68", "0.9"))
  expect_equal(unname(fit$se[, , 1]), apply(manual, c(1, 2), sd))
  expect_equal(
    unname(fit$lower[, , 1, 2]),
    apply(manual, c(1, 2), quantile, probs = 0.05, names = FALSE)
  )
  expect_equal(
    unname(fit$upper[, , 1, 1]),
    apply(manual, c(1, 2), quantile, probs = 0.84, names = FALSE)
  )
  # Normalization holds in every draw: zero uncertainty on impact of variable 1
  expect_equal(unname(fit$se[1, 1, 1]), 0)
  expect_true(all(fit$lower <= fit$upper))

  # Output works with the plotting helpers
  g <- plotirf(fit$irf, fit$se, HTick = 1, Labels = c("Y1", "Y2"))
  expect_length(g, N)
  g2 <- plot2irf(fit$irf, fit$se, fit$irf, fit$se * 0, HTick = 1, Labels = c("Y1", "Y2"))
  expect_length(g2, N)
})

test_that("svariv handles several shocks and deterministic variables", {
  set.seed(7)
  N <- 3
  Nobs <- 2500
  Phi <- array(0, dim = c(N, N, 1))
  Phi[, , 1] <- diag(c(0.5, 0.3, 0.6))
  Psi <- matrix(c(1, 0.5, 0.2, 0, 1, -0.4, 0.3, 0.1, 1), N, N)
  eps <- matrix(rnorm(Nobs * N), ncol = N)
  trend <- (seq_len(Nobs)) / Nobs
  y <- matrix(0, Nobs, N)
  for (t in 2:Nobs) {
    y[t, ] <- c(0.5, 0, -1) * trend[t] + Phi[, , 1] %*% y[t - 1, ] + Psi %*% eps[t, ]
  }
  Z <- cbind(eps[, 1] + rnorm(Nobs), eps[, 2] + rnorm(Nobs))

  fit <- svariv(
    y = y, Z = Z, X = matrix(trend, ncol = 1), P = 1, H = 5, E = 2,
    details = TRUE
  )
  truth <- computeirf(Psi[, 1:2], Phi, 5, cum = FALSE)
  for (e in 1:2) truth[, , e] <- truth[, , e] / truth[1, e, e]

  expect_equal(dim(fit$irf), c(5, N, 2))
  expect_lt(max(abs(fit$irf - truth)), 0.08)
  expect_equal(unname(fit$Psi[1, 2]), 0) # recursive policy block
  expect_equal(diag(fit$Psi[1:2, 1:2]), c(1, 1))
  expect_equal(dim(fit$D), c(N, 1))
  expect_equal(as.numeric(fit$D), c(0.5, 0, -1), tolerance = 0.15)
  expect_length(fit$FStat, 2)
  expect_length(fit$R2, 2)
})

test_that("svariv rejects malformed inputs", {
  dgp <- make_svar_dgp(Nobs = 200)
  expect_error(svariv(y = dgp$y, Z = dgp$Z, P = 0, H = 4), "P must be")
  expect_error(svariv(y = dgp$y, Z = dgp$Z, P = 1, H = 0), "H must be")
  expect_error(svariv(y = dgp$y, Z = dgp$Z, P = 1, H = 4, E = 2), "E must be smaller")
  expect_error(svariv(y = dgp$y, Z = dgp$Z[1:10, , drop = FALSE], P = 1, H = 4), "rows")
  expect_error(svariv(y = dgp$y, Z = dgp$Z, P = 1, H = 4, nboot = -1), "nboot must be")
  expect_error(svariv(y = dgp$y, Z = dgp$Z, P = 1, H = 4, ci = 1.2), "ci must")
  expect_error(svariv(y = dgp$y, Z = dgp$Z, P = 1, H = 4, cum = "yes"), "cum must")
  expect_error(svariv(y = dgp$y, Z = dgp$Z, P = 1, H = 4, details = 1), "details must")
  yna <- dgp$y
  yna[5, 1] <- NA
  expect_error(svariv(y = yna, Z = dgp$Z, P = 1, H = 4), "missing values")
  Zall <- dgp$Z
  Zall[] <- NA
  expect_error(svariv(y = dgp$y, Z = Zall, P = 1, H = 4), "non-missing proxy")
  # Over-identification (more instruments than shocks) is allowed
  set.seed(3)
  Z2 <- cbind(dgp$Z, dgp$Z + rnorm(nrow(dgp$Z)))
  fit <- svariv(y = dgp$y, Z = Z2, P = 1, H = 4, E = 1, details = TRUE)
  expect_equal(dim(fit$irf), c(4, 2, 1))
  expect_length(fit$FStat, 1)
})
