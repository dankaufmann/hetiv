make_fs_data <- function(Tobs = 400, N = 3, E = 2) {
  set.seed(202)
  eps <- matrix(rnorm(Tobs * N), Tobs, N)
  Ind <- rep(0L, Tobs)
  Ind[seq(5, Tobs, by = 5)] <- 1L
  # Event shocks with larger variance on event days and a lower-triangular impact
  Psi <- matrix(c(1, 0.5, 0.3, 0, 1, -0.4, 0, 0, 1), N, N)
  eE <- matrix(rnorm(Tobs * N, sd = 3), Tobs, N) * Ind
  y <- eps + eE %*% t(Psi)
  Z <- (eE[, 1:E] + matrix(rnorm(Tobs * E), Tobs, E)) * Ind
  list(y = y, O = y, Ind = Ind, Z = Z)
}

test_that("hetiv and proxyiv report first-stage F statistics for E = 1", {
  d <- make_fs_data()
  H <- 4
  het <- suppressWarnings(hetiv(y = d$y, O = d$O, Ind = d$Ind, P = 1, H = H, details = TRUE))
  prox <- suppressWarnings(proxyiv(
    y = d$y, O = d$O, Z = d$Z[, 1, drop = FALSE], Ind = d$Ind, P = 1, H = H,
    details = TRUE
  ))

  for (res in list(het, prox)) {
    expect_named(res$FirstStage, c("h", "Tobs", "F_hom", "F_rob"))
    expect_equal(res$FirstStage$h, 0:(H - 1))
    expect_null(res$WeakTest)
    # The sample shrinks with the horizon
    expect_true(all(diff(res$FirstStage$Tobs) <= 0))

    # On impact: F_hom = t^2 and F_rob = HC3 t^2 from the first stage on WeakData
    wd <- res$WeakData
    ctrl <- setdiff(names(wd), c(paste0("y", 1:3), "Z1"))
    fs <- lm(wd$y1 ~ wd$Z1 + as.matrix(wd[, ctrl]))
    t_hom <- coef(fs)[2] / sqrt(vcov(fs)[2, 2])
    t_hc3 <- coef(fs)[2] / sqrt(sandwich::vcovHC(fs, type = "HC3")[2, 2])
    expect_equal(res$FirstStage$F_hom[1], unname(t_hom^2), tolerance = 1e-8)
    expect_equal(res$FirstStage$F_rob[1], unname(t_hc3^2), tolerance = 1e-8)
    expect_equal(res$FirstStage$Tobs[1], stats::nobs(fs))
  }

  # Newey-West version
  het_nw <- suppressWarnings(hetiv(
    y = d$y, O = d$O, Ind = d$Ind, P = 1, H = 2, cov_type = "NW", details = TRUE
  ))
  wd <- het_nw$WeakData
  ctrl <- setdiff(names(wd), c(paste0("y", 1:3), "Z1"))
  fs <- lm(wd$y1 ~ wd$Z1 + as.matrix(wd[, ctrl]))
  t_nw <- coef(fs)[2] /
    sqrt(sandwich::NeweyWest(fs, prewhite = FALSE, adjust = TRUE)[2, 2])
  expect_equal(het_nw$FirstStage$F_rob[1], unname(t_nw^2), tolerance = 1e-8)
  expect_equal(het_nw$FirstStage$F_hom[1], het$FirstStage$F_hom[1])

  # Not computed without details
  expect_null(hetiv(y = d$y, O = d$O, Ind = d$Ind, P = 1, H = 2)$FirstStage)
})

test_that("hetiv and proxyiv run the Lewis-Mertens test by horizon and outcome for E > 1", {
  d <- make_fs_data()
  H <- 3
  N <- 3
  E <- 2
  het <- suppressWarnings(hetiv(y = d$y, O = d$O, Ind = d$Ind, P = 1, H = H, E = E, details = TRUE))
  prox <- suppressWarnings(proxyiv(
    y = d$y, O = d$O, Z = d$Z, Ind = d$Ind, P = 1, H = H, E = E,
    recursive = TRUE, details = TRUE
  ))

  for (res in list(het, prox)) {
    wt <- res$WeakTest
    expect_null(res$FirstStage)
    expect_named(wt, c("h", "i", "Tobs", "gmin", "cv", "cv_simplified", "sy_stat", "sy_cv"))
    # All (h, i) pairs except the policy variables on impact
    expect_equal(nrow(wt), N * H - E)
    expect_false(any(wt$h == 0 & wt$i <= E))
    # The statistic depends on the horizon only, the critical value on (h, i)
    expect_equal(length(unique(round(wt$gmin[wt$h == 1], 8))), 1)

    # Manual call on impact for outcome 3 reproduces the stored row
    wd <- res$WeakData
    ctrl <- setdiff(names(wd), c(paste0("y", 1:N), paste0("Z", 1:E)))
    man <- gweakivtest(
      y = wd$y3, Y = as.matrix(wd[, paste0("y", 1:E)]),
      X = as.matrix(wd[, ctrl]), Z = as.matrix(wd[, paste0("Z", 1:E)]),
      cov_type = "EHW"
    )
    row <- wt[wt$h == 0 & wt$i == 3, ]
    expect_equal(row$gmin, man$gmin_generalized, tolerance = 1e-8)
    expect_equal(row$cv, man$gmin_generalized_critical_value, tolerance = 1e-8)
    expect_equal(row$Tobs, man$nobs)
  }

  expect_null(hetiv(y = d$y, O = d$O, Ind = d$Ind, P = 1, H = 2, E = E)$WeakTest)
})

test_that(".firststage_f handles controls, collinear controls and missing values", {
  set.seed(9)
  Tobs <- 200
  Z <- matrix(rnorm(Tobs * 2), Tobs, 2)
  X <- cbind(rnorm(Tobs), 1) # second column collinear with the constant
  y1 <- Z[, 1] + 0.5 * Z[, 2] + X[, 1] + rnorm(Tobs)
  y1[3] <- NA

  fs <- hetiv:::.firststage_f(y1, Z, X = X, cov_type = "HC1")
  ok <- !is.na(y1)
  reg <- lm(y1[ok] ~ Z[ok, ] + X[ok, 1])
  expect_equal(fs$Tobs, Tobs - 1)
  expect_equal(fs$K, 2)
  expect_equal(fs$F_hom, unname(summary(reg)$fstatistic[1]) * 0 +
    drop(crossprod(coef(reg)[2:3], solve(vcov(reg)[2:3, 2:3], coef(reg)[2:3]))) / 2,
  tolerance = 1e-8
  )
  b <- coef(reg)[2:3]
  V <- sandwich::vcovHC(reg, type = "HC1")[2:3, 2:3]
  Zt <- resid(lm(Z[ok, ] ~ X[ok, 1]))
  Qzz <- crossprod(Zt) / sum(ok)
  expect_equal(fs$F_rob, drop(t(b) %*% Qzz %*% b) / sum(diag(Qzz %*% V)), tolerance = 1e-8)

  # Collinear instruments give NA
  fs_na <- hetiv:::.firststage_f(y1, cbind(Z[, 1], Z[, 1]))
  expect_true(is.na(fs_na$F_hom))
  expect_error(hetiv:::.firststage_f(y1, Z, cov_type = "HC9"), "cov_type")
})
