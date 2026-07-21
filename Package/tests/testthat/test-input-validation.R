test_that("MiXcan2_model rejects incompatible dimensions", {
  dat <- make_mixcan_data()

  expect_error(
    MiXcan2_model(
      y = dat$y[-1], x = dat$x, cov = dat$cov, pi = dat$pi,
      foldid = dat$foldid
    )
  )

  expect_error(
    MiXcan2_model(
      y = dat$y, x = dat$x, cov = dat$cov, pi = dat$pi[-1],
      foldid = dat$foldid
    )
  )
})
