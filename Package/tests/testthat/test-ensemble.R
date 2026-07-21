test_that("MiXcan2_ensemble completes a small smoke test", {
  skip_on_cran()
  skip_if_not_installed("foreach")

  expect_true(
    exists(
      "MiXcan2_ensemble",
      envir = asNamespace("MiXcan2"),
      inherits = FALSE
    )
  )

  # Use a deterministic sequential backend during testing.
  foreach::registerDoSEQ()

  dat <- make_mixcan_data(n = 100L, p = 6L, q = 0L)

  set.seed(123)
  invisible(capture.output(
    out <- MiXcan2_ensemble(
      y = dat$y,
      x = dat$x,
      cov = NULL,
      pi = dat$pi,
      yName = dat$y_name,
      B = 2,
      seed = 123
    )
  ))

  # Check the current broad return contract without assuming
  # undocumented component names.
  expect_type(out, "list")
  expect_gt(length(out), 0L)
})

test_that("MiXcan2_ensemble accepts covariates", {
  skip_on_cran()
  skip_if_not_installed("foreach")

  expect_true(
    exists(
      "MiXcan2_ensemble",
      envir = asNamespace("MiXcan2"),
      inherits = FALSE
    )
  )

  foreach::registerDoSEQ()
  dat <- make_mixcan_data(n = 100L, p = 6L, q = 2L)

  set.seed(456)
  invisible(capture.output(
    out <- MiXcan2_ensemble(
      y = dat$y,
      x = dat$x,
      cov = dat$cov,
      pi = dat$pi,
      yName = dat$y_name,
      B = 2,
      seed = 456
    )
  ))

  expect_type(out, "list")
  expect_gt(length(out), 0L)
})
