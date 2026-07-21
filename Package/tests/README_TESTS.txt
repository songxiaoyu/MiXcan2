MiXcan2 testthat suite aligned with the current GitHub README
==============================================================

Files
-----
tests/testthat.R
tests/testthat/helper-readme-data.R
tests/testthat/test-package.R
tests/testthat/test-helper-data.R
tests/testthat/test-prediction.R
tests/testthat/test-association.R
tests/testthat/test-ensemble.R
tests/testthat/test-readme-pipeline.R

What is tested
--------------
1. Package loading and exported public functions.
2. README-style simulated training and application data.
3. MiXcan2_prediction dimensions and exact matrix multiplication.
4. MiXcan2_association output columns, p-value ranges, signal direction,
   no-covariate support, and mismatched-input handling.
5. MiXcan2_ensemble README arguments and documented output components.
6. Compatibility between ensemble weights and prediction.
7. A complete ensemble -> prediction -> association pipeline.

Run
---
devtools::test()

Notes
-----
- Ensemble tests use foreach::registerDoSEQ() to avoid parallel-worker
  instability in package tests and GitHub Actions.
- Ensemble tests use B = 2 to keep runtime manageable.
- Computational ensemble tests use skip_on_cran().
- Association tests use skip_if_not_installed("ACAT").
- The tests intentionally use the README spelling "ensemble_intercept",
  because that is the currently documented component name.
