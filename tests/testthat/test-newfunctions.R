# Tests for plotseries(), gss_pca(), lag_opt() and rescale_to_match()

test_that("plotseries plots any number of series with one legend", {
  set.seed(42)
  d <- data.frame(t = 1:30, a = rnorm(30), b = rnorm(30), c = rnorm(30))

  p3 <- plotseries(d)
  expect_s3_class(p3, "ggplot")
  expect_equal(nlevels(p3$data$series_id), 3)
  expect_equal(nrow(p3$data), 90)
  expect_identical(p3$theme$legend.position, "bottom")

  # Named series give legend labels and select columns; two series allow the
  # correlation label
  p2 <- plotseries(d, x_var = "t", series = c(b = "Series B", a = "Series A"),
    standardize = TRUE, show_cor = TRUE
  )
  expect_equal(levels(p2$data$series_id), c("b", "a"))
  expect_equal(mean(p2$data$value[p2$data$series_id == "a"]), 0)
  expect_equal(stats::sd(p2$data$value[p2$data$series_id == "a"]), 1)
  lab <- p2$layers[[2]]$aes_params$label
  expect_match(lab, sprintf("Correlation: %.2f", stats::cor(d$a, d$b)), fixed = TRUE)
  expect_match(lab, "T = 30", fixed = TRUE)

  # Legend inside the panel and saving to a file
  f <- tempfile(fileext = ".pdf")
  p <- plotseries(d, legend_position = c(0.5, 0.5), xlim = c(0, 10), file = f)
  expect_identical(p$theme$legend.position, "inside")
  expect_true(file.exists(f))

  # Missing values: by default the line connects across gaps (rows with
  # missing values are dropped); connect = FALSE keeps the gaps
  dm <- d
  dm$a[c(5, 6, 20)] <- NA
  p_con <- plotseries(dm, series = c("a", "b"))
  expect_equal(sum(p_con$data$series_id == "a"), 27)
  expect_false(anyNA(p_con$data$value))
  expect_equal(nrow(ggplot2::layer_data(p_con, 1)), 57)
  p_gap <- plotseries(dm, series = c("a", "b"), connect = FALSE)
  expect_equal(sum(is.na(p_gap$data$value)), 3)

  expect_error(plotseries(d, show_cor = TRUE), "exactly two")
  expect_error(plotseries(d, series = "zz"), "distinct columns")
  expect_error(plotseries(d, x_var = "zz"), "x_var")
  expect_error(plotseries(d, colors = "red"), "at least one entry")
  expect_error(plotseries(d["t"]), "at least one series")
})

test_that("plotseries accepts time series objects through tsbox", {
  skip_if_not_installed("tsbox")
  x <- stats::ts(cumsum(rnorm(24)), start = c(2020, 1), frequency = 12)
  y <- stats::ts(cumsum(rnorm(24)), start = c(2020, 1), frequency = 12)
  p <- plotseries(list(First = x, Second = y), show_cor = TRUE)
  expect_s3_class(p, "ggplot")
  expect_equal(levels(p$data$series_id), c("First", "Second"))
  expect_equal(nrow(p$data), 48)
  # A multivariate ts object works as well
  pm <- plotseries(cbind(First = x, Second = y))
  expect_equal(nlevels(pm$data$series_id), 2)
})

test_that("rescale_to_match matches mean and variance on the overlap", {
  set.seed(42)
  x <- rnorm(40, 2, 3)
  y <- 0.2 * x + rnorm(40)
  x[1:3] <- NA
  y[38:40] <- NA
  out <- rescale_to_match(x, y)
  ok <- !is.na(x) & !is.na(y)
  expect_length(out, 40)
  expect_equal(mean(out[ok]), mean(x[ok]))
  expect_equal(stats::sd(out[ok]), stats::sd(x[ok]))
  expect_true(all(is.na(out[38:40])))
  expect_error(rescale_to_match(1:3, 1:2), "same length")
  expect_error(rescale_to_match(c(1, NA, NA), c(1, 2, 3)), "overlapping")
  expect_error(rescale_to_match(c(1, 2, 3), c(1, 1, 1)), "zero variance")
})

make_panel <- function(n = 300, noise = 0.05) {
  set.seed(42)
  f <- matrix(rnorm(2 * n), ncol = 2)
  L <- cbind(c(1, 0.8, 0.6, 0.4, 0.2), c(0, 0.3, 0.6, 0.8, 1))
  X <- f %*% t(L) + matrix(rnorm(5 * n, sd = noise), ncol = 5)
  colnames(X) <- paste0("r", 1:5)
  list(X = X, f = f)
}

test_that("gss_pca recovers target and path factors", {
  d <- make_panel()
  res <- gss_pca(d$X, n_factors = 2)
  expect_s3_class(res, "gss_pca")
  expect_equal(colnames(res$factors), c("target", "path"))
  expect_equal(dim(res$loadings), c(5L, 2L))
  # Path has no loading on the first series; factors match the true ones
  expect_lt(abs(res$restriction), 1e-10)
  expect_gt(abs(stats::cor(res$factors[, "target"], d$f[, 1])), 0.98)
  expect_gt(abs(stats::cor(res$factors[, "path"], d$f[, 2])), 0.98)
  # Normalisation: unit coefficient of the first series on the target factor,
  # equal coefficients of both factors on the second series
  expect_equal(unname(coef(lm(d$X[, 1] ~ res$factors[, "target"]))[2]), 1)
  expect_equal(
    unname(coef(lm(d$X[, 2] ~ res$factors[, "path"]))[2]),
    unname(coef(lm(d$X[, 2] ~ res$factors[, "target"]))[2])
  )
  expect_true(res$converged)
  expect_output(print(res), "gss_pca: 2 factor")

  # One factor: first principal component, unit coefficient on the anchor
  r1 <- gss_pca(d$X, n_factors = 1)
  expect_equal(colnames(r1$factors), "overall")
  expect_equal(unname(coef(lm(d$X[, 1] ~ r1$factors[, 1]))[2]), 1)
  pc1 <- stats::prcomp(d$X, scale. = TRUE)$x[, 1]
  expect_equal(abs(stats::cor(r1$factors[, 1], pc1)), 1)
  expect_null(r1$restriction)
})

test_that("gss_pca handles missing values and validates inputs", {
  d <- make_panel()
  full <- gss_pca(d$X, n_factors = 2)
  X <- d$X
  X[sample(length(X), 60)] <- NA
  res <- gss_pca(X, n_factors = 2)
  expect_equal(res$n_missing, 60)
  expect_true(res$converged)
  expect_false(anyNA(res$factors))
  expect_gt(stats::cor(res$factors[, "target"], full$factors[, "target"]), 0.98)

  expect_error(gss_pca(d$X, n_factors = 3), "n_factors")
  expect_error(gss_pca(d$X, norm_anchor = c(1, 1)), "cannot be 1")
  expect_error(gss_pca(d$X[, 1:2]), "more series")
  expect_error(gss_pca(cbind(d$X, 1)), "Constant")
  X2 <- d$X
  X2[5, ] <- NA
  expect_warning(gss_pca(X2), "no observed series")
})

test_that("lag_opt selects the lag length of a VAR", {
  expect_error(lag_opt(matrix(rnorm(40), ncol = 2), crit = "BIC"), "crit")
  skip_if_not_installed("vars")
  set.seed(42)
  Y <- matrix(rnorm(1000), ncol = 2)
  for (t in 3:500) Y[t, ] <- 0.5 * Y[t - 1, ] + 0.3 * Y[t - 2, ] + Y[t, ]
  expect_equal(lag_opt(Y, crit = "SC", lag.max = 6), 2L)
  X <- matrix(rnorm(500), ncol = 1)
  expect_type(lag_opt(Y, X = X, lag.max = 4), "integer")
})
