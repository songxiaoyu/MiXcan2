test_that("MiXcan2_ensemble follows the README training interface", {
  skip_on_cran()
  skip_if_not_installed("foreach")

  foreach::registerDoSEQ()

  dat <- make_readme_data(
    N = 100L,
    M = 100L,
    P = 8L,
    Q1 = 3L,
    Q2 = 2L,
    seed = 123L
  )

  set.seed(123L)

  invisible(capture.output(
    ensbl <- MiXcan2_ensemble(
      y = dat$y_train,
      x = dat$x_train,
      cov = dat$cov_train,
      pi = dat$pi_train,
      xNameMatrix = dat$x_annotation,
      yName = dat$gene_name,
      B = 2L,
      seed = 123L
    )
  ))

  expect_type(ensbl, "list")

  expected_components <- c(
    "ensemble_intecept",
    "ensemble_weight",
    "all_weights",
    "ensemble_summary_by_type",
    "all_summary"
  )

  expect_true(all(expected_components %in% names(ensbl)))
})


test_that("MiXcan2_ensemble returns documented weight columns", {
  skip_on_cran()
  skip_if_not_installed("foreach")

  foreach::registerDoSEQ()

  dat <- make_readme_data(
    N = 100L,
    M = 100L,
    P = 8L,
    seed = 321L
  )

  invisible(capture.output(
    ensbl <- MiXcan2_ensemble(
      y = dat$y_train,
      x = dat$x_train,
      cov = dat$cov_train,
      pi = dat$pi_train,
      xNameMatrix = dat$x_annotation,
      yName = dat$gene_name,
      B = 2L,
      seed = 321L
    )
  ))

  expect_s3_class(ensbl$ensemble_weight, "data.frame")

  expected_weight_columns <- c(
    "yName",
    "varID",
    "position",
    "rsid",
    "ref",
    "eff",
    "type",
    "weight_cell_1",
    "weight_cell_2"
  )

  expect_true(
    all(expected_weight_columns %in% names(ensbl$ensemble_weight))
  )

  expect_true(
    all(
      is.finite(
        ensbl$ensemble_weight$weight_cell_1
      )
    )
  )

  expect_true(
    all(
      is.finite(
        ensbl$ensemble_weight$weight_cell_2
      )
    )
  )
})


test_that("ensemble weights can be passed to MiXcan2_prediction", {
  skip_on_cran()
  skip_if_not_installed("foreach")

  foreach::registerDoSEQ()

  dat <- make_readme_data(
    N = 100L,
    M = 120L,
    P = 8L,
    seed = 456L
  )

  invisible(capture.output(
    ensbl <- MiXcan2_ensemble(
      y = dat$y_train,
      x = dat$x_train,
      cov = dat$cov_train,
      pi = dat$pi_train,
      xNameMatrix = dat$x_annotation,
      yName = dat$gene_name,
      B = 2L,
      seed = 456L
    )
  ))

  selected_idx <- match(
    ensbl$ensemble_weight$varID,
    dat$x_annotation$varID
  )

  expect_false(anyNA(selected_idx))

  pred <- MiXcan2_prediction(
    weight = as.matrix(
      ensbl$ensemble_weight[
        ,
        c("weight_cell_1", "weight_cell_2"),
        drop = FALSE
      ]
    ),
    new_x = dat$x_new[, selected_idx, drop = FALSE]
  )

  expect_equal(nrow(pred), nrow(dat$x_new))
  expect_equal(ncol(pred), 2L)
  expect_true(all(is.finite(as.matrix(pred))))
})
