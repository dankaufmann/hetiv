make_irf <- function(H = 4, N = 2, E = 1) {
  irf <- array(seq_len(H * N * E) / 10, dim = c(H, N, E))
  dimnames(irf)[[1]] <- seq_len(H) - 1
  irf
}

# Extract ribbon and line layers of a ggplot for inspection
layer_classes <- function(g) {
  vapply(g$layers, function(l) class(l$geom)[1], character(1))
}

test_that("plotirf draws a line only when no uncertainty is supplied", {
  irf <- make_irf()
  g <- plotirf(irf, HTick = 1, Labels = c("a", "b"))

  expect_length(g, 2)
  expect_s3_class(g[[1]], "ggplot")
  expect_false("GeomRibbon" %in% layer_classes(g[[1]]))
  expect_equal(sum(layer_classes(g[[1]]) == "GeomLine"), 1)
})

test_that("plotirf draws normal bands for up to two confidence levels", {
  irf <- make_irf()
  se <- array(0.1, dim = dim(irf), dimnames = dimnames(irf))

  g <- plotirf(irf, se, HTick = 1, Labels = c("a", "b"), ci = c(0.90, 0.95))
  expect_equal(sum(layer_classes(g[[1]]) == "GeomRibbon"), 2)

  # Widest band is drawn first and matches irf +/- z * se
  ribbon <- g[[1]]$layers[[which(layer_classes(g[[1]]) == "GeomRibbon")[1]]]$data
  z <- qnorm((1 + 0.95) / 2)
  expect_equal(ribbon$lower, unname(irf[, 1, 1]) - z * 0.1)
  expect_equal(ribbon$upper, unname(irf[, 1, 1]) + z * 0.1)

  expect_error(
    plotirf(irf, se, HTick = 1, Labels = c("a", "b"), ci = c(0.9, 0.95, 0.99)),
    "one or two confidence levels"
  )
  # Zero standard errors suppress the bands
  g0 <- plotirf(irf, se * 0, HTick = 1, Labels = c("a", "b"))
  expect_false("GeomRibbon" %in% layer_classes(g0[[1]]))
})

test_that("plotirf computes recentred percentile bands from bootstrap draws", {
  irf <- make_irf()
  set.seed(3)
  B <- 300
  # Skewed draws that are not centred on irf: band must still be centred on irf
  boot <- array(rexp(length(irf) * B), dim = c(dim(irf), B)) + 1
  res <- list(irf = irf, se = NULL, boot = boot, Method = "Test")

  g <- plotirf(res, HTick = 1, Labels = c("a", "b"), ci = 0.90)
  ribbon <- g[[1]]$layers[[which(layer_classes(g[[1]]) == "GeomRibbon")[1]]]$data

  dev <- sweep(boot[, 1, 1, ], 1, rowMeans(boot[, 1, 1, ]), "-")
  expect_equal(ribbon$lower, unname(irf[, 1, 1]) + apply(dev, 1, quantile, 0.05, names = FALSE))
  expect_equal(ribbon$upper, unname(irf[, 1, 1]) + apply(dev, 1, quantile, 0.95, names = FALSE))
  # Recentring: the band is not the raw percentile interval
  expect_false(isTRUE(all.equal(
    ribbon$lower, apply(boot[, 1, 1, ], 1, quantile, 0.05, names = FALSE)
  )))

  # boot_type = "percentile" gives the raw percentile interval instead
  gp <- plotirf(res, HTick = 1, Labels = c("a", "b"), ci = 0.90, boot_type = "percentile")
  ribbonp <- gp[[1]]$layers[[which(layer_classes(gp[[1]]) == "GeomRibbon")[1]]]$data
  expect_equal(ribbonp$lower, apply(boot[, 1, 1, ], 1, quantile, 0.05, names = FALSE))
  expect_equal(ribbonp$upper, apply(boot[, 1, 1, ], 1, quantile, 0.95, names = FALSE))
  expect_error(plotirf(res, HTick = 1, Labels = c("a", "b"), boot_type = "foo"))

  # Draws take precedence over standard errors; failed draws (NA) are ignored
  boot_na <- boot
  boot_na[, , , 1:5] <- NA
  res2 <- list(irf = irf, se = irf * 0 + 1, boot = boot_na)
  g2 <- plotirf(res2, HTick = 1, Labels = c("a", "b"), ci = 0.90)
  ribbon2 <- g2[[1]]$layers[[which(layer_classes(g2[[1]]) == "GeomRibbon")[1]]]$data
  expect_false(anyNA(ribbon2$lower))
  expect_false(isTRUE(all.equal(ribbon2$lower, irf[, 1, 1] - qnorm(0.95))))
})

test_that("plotirf overlays up to three IRFs and plot2irf is a wrapper", {
  irf <- make_irf(E = 2)
  se <- array(0.1, dim = dim(irf), dimnames = dimnames(irf))
  labels <- c("a", "b")

  g3 <- plotirf(
    list(list(irf = irf, se = se), irf * 0.8, list(irf = irf * 0.5, Method = "M")),
    HTick = 1, Labels = labels, ci = 0.90, legend = TRUE
  )
  expect_length(g3, 4)
  expect_equal(sum(layer_classes(g3[[1]]) == "GeomLine"), 3)
  expect_equal(sum(layer_classes(g3[[1]]) == "GeomRibbon"), 1)
  expect_identical(g3[[1]]$theme$legend.position, "bottom")

  expect_error(
    plotirf(list(irf, irf, irf, irf), HTick = 1, Labels = labels),
    "up to three"
  )
  expect_error(
    plotirf(list(irf, se), IRFse = se, HTick = 1, Labels = labels),
    "IRFse cannot"
  )
  expect_error(
    plotirf(list(irf, irf[1:3, , ]), HTick = 1, Labels = labels),
    "same dimensions"
  )

  p_old <- plot2irf(irf, se, irf * 0.8, se * 0, HTick = 1, Labels = labels)
  p_new <- plotirf(list(list(irf = irf, se = se), irf * 0.8),
    HTick = 1, Labels = labels, ci = 0.90
  )
  expect_length(p_old, 4)
  expect_equal(layer_classes(p_old[[1]]), layer_classes(p_new[[1]]))
  expect_equal(p_old[[1]]$layers[[1]]$data$lower, p_new[[1]]$layers[[1]]$data$lower)
  expect_identical(p_old[[1]]$theme$legend.position, "none")
})

test_that("plotirf accepts two-dimensional input with standard errors", {
  irf <- make_irf()[, , 1]
  se <- irf * 0 + 0.1
  g <- plotirf(irf, se, HTick = 1, Labels = c("a", "b"))
  expect_length(g, 2)
  expect_true("GeomRibbon" %in% layer_classes(g[[1]]))

  boot <- array(irf, dim = c(dim(irf), 20))
  g2 <- plotirf(list(irf = irf, boot = boot), HTick = 1, Labels = c("a", "b"))
  expect_length(g2, 2)
})

test_that("plotirf validates names, colours and inputs", {
  irf <- make_irf()
  expect_error(plotirf(irf, HTick = 1, Labels = "only one"), "Labels")
  expect_error(plotirf(list(a = 1), HTick = 1, Labels = c("a", "b")), "element 'irf'")
  expect_error(
    plotirf(list(irf, irf), HTick = 1, Labels = c("a", "b"), names = "one"),
    "names must"
  )
  expect_error(
    plotirf(list(irf, irf), HTick = 1, Labels = c("a", "b"), colors = "red"),
    "colors must"
  )
  bad_boot <- array(1, dim = c(3, 2, 1, 5))
  expect_error(
    plotirf(list(irf = irf, boot = bad_boot), HTick = 1, Labels = c("a", "b")),
    "Bootstrap draws"
  )
})
