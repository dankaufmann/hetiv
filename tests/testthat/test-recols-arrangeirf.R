make_recursive_data <- function(Tobs = 2000) {
  set.seed(42)
  e1 <- rnorm(Tobs)
  e2 <- rnorm(Tobs)
  y <- cbind(e1, 0.5 * e1 + e2)
  Ind <- rep(0L, Tobs)
  Ind[seq(5, Tobs, by = 5)] <- 1L
  list(y = y, Ind = Ind)
}

test_that("recols imposes recursive zero restrictions and unit normalization", {
  d <- make_recursive_data(200)
  fit <- recols(y = d$y, O = d$y, Ind = d$Ind, P = 1, H = 3, E = 2)

  expect_equal(dim(fit$irf), c(3, 2, 2))
  expect_equal(dimnames(fit$irf)[[1]], c("0", "1", "2"))
  expect_equal(unname(fit$irf[1, 1, 1]), 1)
  expect_equal(unname(fit$irf[1, 2, 2]), 1)
  expect_equal(unname(fit$irf[1, 1, 2]), 0, tolerance = 1e-10)
  expect_identical(fit$Method, "Recursive-OLS")
})

test_that("recols recovers the impact response of a recursive DGP", {
  d <- make_recursive_data()
  fit <- recols(y = d$y, O = d$y, Ind = d$Ind, P = 1, H = 2, norm = 2)

  expect_equal(unname(fit$irf[1, 1, 1]), 2)
  expect_equal(unname(fit$irf[1, 2, 1]), 1, tolerance = 0.1)
})

test_that("recols returns details, excludes contaminated days, and validates input", {
  d <- make_recursive_data(200)
  Ind <- d$Ind
  Ind[seq(3, 200, by = 10)] <- 2L
  # iid data: kfpredict() warns about weak heteroskedasticity, which is expected here
  fit <- suppressWarnings(
    recols(y = d$y, O = d$y, Ind = Ind, P = 1, H = 3, details = TRUE)
  )

  expect_named(fit, c(
    "irf", "se", "Shocks", "OLSRes", "Obs", "Method",
    "et", "Sig", "SigR", "Psi"
  ))
  expect_equal(fit$Obs$Tp, sum(Ind[-1] == 1))
  expect_equal(fit$Obs$To, sum(Ind[-1] == 2))
  expect_equal(fit$Obs$Tt, fit$Obs$Tp + fit$Obs$Tc)
  expect_equal(length(fit$OLSRes), 3 * 2)
  expect_equal(stats::nobs(fit$OLSRes[["OLS.h1.n1.e1"]]), sum(Ind[-1] < 2))
  expect_true(all(is.na(fit$et[Ind != 1, ])))
  expect_equal(as.numeric(fit$Psi), as.numeric(fit$irf[1, , 1]))

  expect_error(
    recols(y = d$y, O = d$y, Ind = d$Ind, P = 1, H = 5, Hstep = 2,
      hpredict = 2, details = TRUE),
    "hpredict"
  )
  expect_error(recols(y = d$y, O = d$y, Ind = d$Ind, P = 1, H = 3, E = 3), "E cannot")
})

test_that("recols cumulative and level responses agree on impact", {
  d <- make_recursive_data(200)
  lev <- recols(y = d$y, O = d$y, Ind = d$Ind, P = 1, H = 3)
  cum <- recols(y = d$y, O = d$y, Ind = d$Ind, P = 1, H = 3, cum = c(TRUE, FALSE))

  expect_equal(cum$irf[1, , ], lev$irf[1, , ])
  expect_equal(cum$irf[, 2, ], lev$irf[, 2, ])
})

test_that("arrangeirf returns one figure per shock and saves PDFs", {
  irf <- array(c(1, 0.5, 0.2, 0, 0.3, 0.1, 1, 0.4, 0.2, 0, 0.1, 0.1), dim = c(3, 2, 2))
  dimnames(irf)[[1]] <- 0:2
  se <- array(0.1, dim = dim(irf), dimnames = dimnames(irf))
  g <- plotirf(irf, se, HTick = 1, Labels = c("A", "B"))

  stem <- file.path(tempdir(), "arrangeirf_test")
  pdf(NULL)
  figs <- arrangeirf(g, E = 2, N = 2, nRows = 1, file = stem)
  grDevices::dev.off()

  expect_length(figs, 2)
  expect_s3_class(figs[[1]], "gtable")
  expect_true(file.exists(paste0(stem, "_1.pdf")))
  expect_true(file.exists(paste0(stem, "_2.pdf")))

  expect_error(arrangeirf(g, E = 1, N = 2, nRows = 1), "length E \\* N")
  expect_error(arrangeirf(g, E = 2, N = 2, nRows = 3), "nRows")
})
