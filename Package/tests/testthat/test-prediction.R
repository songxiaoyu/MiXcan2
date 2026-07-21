test_that("MiXcan2_prediction returns two cell-type predictions", {
  dat <- make_prediction_inputs(n = 80L, p = 4L)

  pred <- MiXcan2_prediction(
    weight = dat$weight,
    new_x = dat$new_x
  )

  expect_true(is.matrix(pred) || is.data.frame(pred))
  expect_equal(dim(pred), c(80L, 2L))
  expect_true(all(is.finite(as.matrix(pred))))
})


test_that("MiXcan2_prediction agrees with matrix multiplication", {
  dat <- make_prediction_inputs(n = 50L, p = 3L)

  pred <- MiXcan2_prediction(
    weight = dat$weight,
    new_x = dat$new_x
  )

  expected <- dat$new_x %*% dat$weight

  expect_equal(
    unname(as.matrix(pred)),
    unname(expected),
    tolerance = 1e-10
  )
})


test_that("MiXcan2_prediction rejects incompatible dimensions", {
  dat <- make_prediction_inputs(n = 30L, p = 4L)

  bad_weight <- dat$weight[-1, , drop = FALSE]

  expect_error(
    MiXcan2_prediction(
      weight = bad_weight,
      new_x = dat$new_x
    )
  )
})
