test_that("MiXcan2 package loads and exposes functions", {
  expect_true(requireNamespace("MiXcan2", quietly = TRUE))

  exports <- getNamespaceExports("MiXcan2")
  expect_type(exports, "character")
  expect_gt(length(exports), 0L)

  exported_objects <- mget(
    exports,
    envir = asNamespace("MiXcan2"),
    inherits = FALSE
  )
  exported_functions <- vapply(exported_objects, is.function, logical(1))

  expect_true(any(exported_functions))
})

test_that("MiXcan2 version is available", {
  version <- utils::packageVersion("MiXcan2")
  expect_s3_class(version, "package_version")
})
