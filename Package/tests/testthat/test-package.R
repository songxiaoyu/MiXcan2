test_that("MiXcan2 loads and documents the public pipeline functions", {
  expect_true(requireNamespace("MiXcan2", quietly = TRUE))

  exports <- getNamespaceExports("MiXcan2")

  expect_true(
    all(
      c(
        "MiXcan2_ensemble",
        "MiXcan2_prediction",
        "MiXcan2_association"
      ) %in% exports
    )
  )
})


test_that("MiXcan2 has a valid package version", {
  version <- utils::packageVersion("MiXcan2")

  expect_s3_class(version, "package_version")
  expect_true(version >= package_version("0.1.0"))
})
