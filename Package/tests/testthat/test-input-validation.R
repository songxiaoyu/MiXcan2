test_that("test data generator returns aligned inputs", {
  dat <- make_mixcan_data(n = 60L, p = 5L, q = 2L)

  expect_length(dat$y, 60L)
  expect_equal(dim(dat$x), c(60L, 5L))
  expect_equal(dim(dat$cov), c(60L, 2L))
  expect_length(dat$pi, 60L)

  expect_true(all(is.finite(dat$y)))
  expect_true(all(is.finite(dat$x)))
  expect_true(all(is.finite(dat$cov)))
  expect_true(all(dat$pi > 0 & dat$pi < 1))
})
