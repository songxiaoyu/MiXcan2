test_that("MiXcan2_ensemble completes a small smoke test", {
  foreach::registerDoSEQ()

  out <- MiXcan2_ensemble(
    # existing arguments
  )

  # assertions below
})

test_that("MiXcan2_ensemble completes a small smoke test", {
  skip_on_cran()
  dat <- make_mixcan_data(n = 100L, p = 6L, q = 0L)

  set.seed(123)
  out <- MiXcan2_ensemble(
    y = dat$y,
    x = dat$x,
    cov = NULL,
    pi = dat$pi,
    yName = dat$y_name,
    B = 2,
    seed = 123
  )

  expect_type(out, "list")
  expect_true(all(c("summary", "ensemble_weight", "all_weights") %in% names(out)))
  expect_s3_class(out$summary, "data.frame")
  expect_s3_class(out$ensemble_weight, "data.frame")
  expect_s3_class(out$all_weights, "data.frame")
})
