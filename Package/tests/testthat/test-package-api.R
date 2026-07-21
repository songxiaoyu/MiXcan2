test_that("documented public functions are exported", {
  expected <- c(
    "MiXcan2_association",
    "MiXcan2_ensemble",
    "MiXcan2_extract_summary",
    "MiXcan2_extract_weight",
    "MiXcan2_model",
    "MiXcan2_prediction",
    "deNet_purity",
    "pi_estimation"
  )

  exports <- getNamespaceExports("MiXcan2")
  expect_setequal(intersect(exports, expected), expected)
})
