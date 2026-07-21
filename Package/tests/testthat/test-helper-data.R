test_that("README-style simulated data have compatible dimensions", {
  dat <- make_readme_data(
    N = 60L,
    M = 100L,
    P = 6L,
    Q1 = 3L,
    Q2 = 2L
  )

  expect_length(dat$y_train, 60L)
  expect_length(dat$pi_train, 60L)
  expect_equal(dim(dat$x_train), c(60L, 6L))
  expect_equal(dim(dat$cov_train), c(60L, 3L))

  expect_equal(dim(dat$x_new), c(100L, 6L))
  expect_equal(dim(dat$cov_new), c(100L, 2L))
  expect_length(dat$outcome, 100L)

  expect_equal(colnames(dat$x_train), colnames(dat$x_new))
  expect_equal(nrow(dat$x_annotation), ncol(dat$x_train))

  expect_true(all(dat$pi_train > 0 & dat$pi_train < 1))
  expect_true(all(is.finite(dat$y_train)))
  expect_true(all(is.finite(dat$outcome)))
})


test_that("SNP annotation contains fields used by MiXcan2_ensemble", {
  dat <- make_readme_data(P = 5L)

  expect_named(
    dat$x_annotation,
    c("varID", "position", "rsid", "ref", "eff"),
    ignore.order = FALSE
  )

  expect_equal(length(unique(dat$x_annotation$varID)), 5L)
  expect_true(all(is.finite(dat$x_annotation$position)))
})
