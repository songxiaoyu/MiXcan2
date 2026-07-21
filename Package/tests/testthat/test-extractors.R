test_that("MiXcan2_extract_weight returns a usable weight table", {
  fit <- fit_small_mixcan_model()
  out <- MiXcan2_extract_weight(model = fit$model)

  expect_s3_class(out, "data.frame")
  expect_true(all(c(
    "yName", "xNameMatrix", "weight_cell_1", "weight_cell_2", "type"
  ) %in% names(out)))
  expect_true(all(out$yName == fit$data$y_name))
  expect_true(all(out$type %in% c("CellTypeSpecific", "NonSpecific", "NoPredictor")))
  expect_true(all(is.finite(out$weight_cell_1)))
  expect_true(all(is.finite(out$weight_cell_2)))
})

test_that("MiXcan2_extract_summary returns finite model diagnostics", {
  fit <- fit_small_mixcan_model()
  dat <- fit$data

  out <- MiXcan2_extract_summary(
    x = dat$x,
    y = dat$y,
    cov = dat$cov,
    pi = dat$pi,
    model = fit$model
  )

  expect_s3_class(out, "data.frame")
  expect_equal(nrow(out), 1L)
  expect_true(all(c(
    "yName", "n_snp_input", "n_snp_model", "model_type",
    "cv_r2", "in_sample_r2"
  ) %in% names(out)))
  expect_equal(out$yName, dat$y_name)
  expect_equal(out$n_snp_input, ncol(dat$x))
  expect_true(out$n_snp_model >= 0)
  expect_true(out$model_type %in% c("CellTypeSpecific", "NonSpecific", "NoPredictor"))
  expect_true(is.finite(out$cv_r2))
  expect_true(is.finite(out$in_sample_r2))
})
