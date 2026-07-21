make_mixcan_data <- function(n = 100L, p = 6L, q = 0L, seed = 123L) {
  stopifnot(n > 10L, p > 1L, q >= 0L)

  set.seed(seed)

  x <- matrix(
    sample(0:2, n * p, replace = TRUE),
    nrow = n,
    ncol = p
  )
  colnames(x) <- paste0("V", seq_len(p))

  pi <- stats::runif(n, min = 0.15, max = 0.85)

  cov <- if (q > 0L) {
    out <- matrix(stats::rnorm(n * q), nrow = n, ncol = q)
    colnames(out) <- paste0("cov", seq_len(q))
    out
  } else {
    NULL
  }

  beta <- seq(0.20, 0.05, length.out = p)
  y <- as.numeric(x %*% beta + 0.5 * pi + stats::rnorm(n, sd = 0.5))

  list(
    y = y,
    x = x,
    cov = cov,
    pi = pi,
    y_name = "test_gene"
  )
}
