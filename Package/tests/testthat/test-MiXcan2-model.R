test_that("MiXcan2_model returns the documented model components", {
  fit <- fit_small_mixcan_model()
  model <- fit$model
  dat <- fit$data

  expect_type(model, "list")
  expect_true(model$type %in% c("CellTypeSpecific", "NonSpecific", "NoPredictor"))

  expect_s3_class(model$beta.SNP.cell1, "data.frame")
  expect_s3_class(model$beta.SNP.cell2, "data.frame")
  expect_equal(nrow(model$beta.SNP.cell1), ncol(dat$x))
  expect_equal(nrow(model$beta.SNP.cell2), ncol(dat$x))
  expect_named(model$beta.SNP.cell1, c("xNameMatrix", "weight"))
  expect_named(model$beta.SNP.cell2, c("xNameMatrix", "weight"))
  expect_equal(model$beta.SNP.cell1$xNameMatrix, dat$x_names)
  expect_equal(model$beta.SNP.cell2$xNameMatrix, dat$x_names)
  expect_true(all(is.finite(model$beta.SNP.cell1$weight)))
  expect_true(all(is.finite(model$beta.SNP.cell2$weight)))

  expect_true(is.matrix(model$beta.all.models))
  expect_equal(colnames(model$beta.all.models), c("Tissue", "Cell1", "Cell2"))
  expect_length(model$intercept, 2L)
  expect_equal(model$yName, dat$y_name)
})

test_that("MiXcan2_model accepts covariates", {
  fit <- fit_small_mixcan_model(with_covariates = TRUE)
  expect_type(fit$model, "list")
  expect_true(fit$model$type %in% c("CellTypeSpecific", "NonSpecific", "NoPredictor"))
  expect_equal(nrow(fit$model$beta.SNP.cell1), ncol(fit$data$x))
})

test_that("MiXcan2_model is reproducible with fixed folds and seed", {
  fit1 <- fit_small_mixcan_model(seed = 999L)
  fit2 <- fit_small_mixcan_model(seed = 999L)

  expect_equal(fit1$model$type, fit2$model$type)
  expect_equal(fit1$model$beta.SNP.cell1, fit2$model$beta.SNP.cell1, tolerance = 1e-10)
  expect_equal(fit1$model$beta.SNP.cell2, fit2$model$beta.SNP.cell2, tolerance = 1e-10)
})
