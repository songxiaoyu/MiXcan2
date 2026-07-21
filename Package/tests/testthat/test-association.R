test_that("MiXcan2_association returns documented association summaries", {
  skip_if_not_installed("ACAT")

  dat <- make_association_inputs(n = 300L)

  result <- MiXcan2_association(
    new_y = dat$new_y,
    new_cov = dat$new_cov,
    new_outcome = dat$outcome,
    family = gaussian
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)

  documented_names <- c(
    "cell1_est",
    "cell1_se",
    "cell1_p",
    "cell2_est",
    "cell2_se",
    "cell2_p",
    "p_combined"
  )

  expect_true(all(documented_names %in% names(result)))

  numeric_results <- as.matrix(
    result[, documented_names, drop = FALSE]
  )
  storage.mode(numeric_results) <- "double"

  expect_true(all(is.finite(numeric_results)))
  expect_true(result$cell1_p >= 0 && result$cell1_p <= 1)
  expect_true(result$cell2_p >= 0 && result$cell2_p <= 1)
  expect_true(result$p_combined >= 0 && result$p_combined <= 1)
})


test_that("MiXcan2_association detects a strong simulated cell-1 signal", {
  skip_if_not_installed("ACAT")

  dat <- make_association_inputs(n = 500L, seed = 11L)

  result <- MiXcan2_association(
    new_y = dat$new_y,
    new_cov = dat$new_cov,
    new_outcome = dat$outcome,
    family = gaussian
  )

  expect_gt(result$cell1_est, 0)
  expect_lt(result$cell1_p, 0.01)
  expect_lt(result$p_combined, 0.05)
})


test_that("MiXcan2_association works without covariates", {
  skip_if_not_installed("ACAT")

  dat <- make_association_inputs(n = 250L, seed = 22L)

  result <- MiXcan2_association(
    new_y = dat$new_y,
    new_cov = NULL,
    new_outcome = dat$outcome,
    family = gaussian
  )

  expect_s3_class(result, "data.frame")
  expect_equal(nrow(result), 1L)

  numeric_result <- as.matrix(result)
  storage.mode(numeric_result) <- "double"

  expect_true(all(is.finite(numeric_result)))
})


test_that("MiXcan2_association rejects mismatched sample sizes", {
  skip_if_not_installed("ACAT")

  dat <- make_association_inputs(n = 100L)

  expect_error(
    MiXcan2_association(
      new_y = dat$new_y,
      new_cov = dat$new_cov[-1, , drop = FALSE],
      new_outcome = dat$outcome,
      family = gaussian
    )
  )
})
