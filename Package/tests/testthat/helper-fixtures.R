make_mixcan_data <- function(n = 120L, p = 8L, q = 2L, seed = 123L) {
  stopifnot(n >= 30L, p >= 2L, q >= 0L)
  set.seed(seed)

  x <- matrix(rbinom(n * p, size = 2, prob = 0.30), nrow = n, ncol = p)
  colnames(x) <- paste0("SNP", seq_len(p))

  pi <- stats::rbeta(n, shape1 = 2, shape2 = 3)
  cov <- if (q > 0L) {
    out <- matrix(stats::rnorm(n * q), nrow = n, ncol = q)
    colnames(out) <- paste0("Cov", seq_len(q))
    out
  } else {
    NULL
  }

  # A reproducible signal that varies by cell type.
  y <- 0.8 * x[, 1] * pi + 0.2 * x[, 1] * (1 - pi) +
    0.15 * x[, 2] + stats::rnorm(n, sd = 0.35)

  foldid <- rep(seq_len(10L), length.out = n)

  list(
    y = y,
    x = x,
    cov = cov,
    pi = pi,
    foldid = foldid,
    x_names = colnames(x),
    y_name = "Gene1"
  )
}

fit_small_mixcan_model <- function(seed = 123L, with_covariates = FALSE) {
  dat <- make_mixcan_data(seed = seed, q = if (with_covariates) 2L else 0L)
  set.seed(seed)
  model <- MiXcan2_model(
    y = dat$y,
    x = dat$x,
    cov = dat$cov,
    pi = dat$pi,
    xNameMatrix = dat$x_names,
    yName = dat$y_name,
    foldid = dat$foldid
  )
  list(data = dat, model = model)
}
