test_that("the README-style end-to-end pipeline completes", {
  skip_on_cran()
  skip_if_not_installed("foreach")
  skip_if_not_installed("ACAT")

  foreach::registerDoSEQ()

  dat <- make_readme_data(
    N = 100L,
    M = 200L,
    P = 8L,
    Q1 = 3L,
    Q2 = 2L,
    seed = 789L
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
      seed = 789L
    )
  ))

  selected_idx <- match(
    ensbl$ensemble_weight$varID,
    dat$x_annotation$varID
  )

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

  association <- MiXcan2_association(
    new_y = pred,
    new_cov = dat$cov_new,
    new_outcome = dat$outcome,
    family = gaussian
  )

  expect_equal(dim(pred), c(200L, 2L))
  expect_s3_class(association, "data.frame")
  expect_equal(nrow(association), 1L)
  expect_true(
    all(
      c(
        "cell1_est",
        "cell1_se",
        "cell1_p",
        "cell2_est",
        "cell2_se",
        "cell2_p",
        "p_combined"
      ) %in% names(association)
    )
  )
})
