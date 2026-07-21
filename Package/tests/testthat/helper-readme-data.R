make_readme_data <- function(
    N = 100L,
    M = 200L,
    P = 10L,
    Q1 = 3L,
    Q2 = 2L,
    seed = 123L
) {
  stopifnot(
    N >= 30L,
    M >= 30L,
    P >= 3L,
    Q1 >= 0L,
    Q2 >= 0L
  )

  set.seed(seed)

  gene_name <- "BRCA1"

  # Training data
  pi_train <- stats::runif(N, min = 0.1, max = 0.9)

  x_train <- matrix(
    sample(0:2, N * P, replace = TRUE),
    nrow = N,
    ncol = P
  )

  snp_names <- paste0("V", seq_len(P))
  colnames(x_train) <- snp_names

  x_annotation <- data.frame(
    varID = paste0(
      "chr21_",
      seq.int(1000001L, by = 1000L, length.out = P),
      "_A_G_b38"
    ),
    position = seq.int(1000001L, by = 1000L, length.out = P),
    rsid = paste0("rs", seq_len(P) + 1000L),
    ref = rep("A", P),
    eff = rep("G", P),
    stringsAsFactors = FALSE
  )

  cov_train <- if (Q1 > 0L) {
    out <- matrix(stats::rnorm(N * Q1), nrow = N, ncol = Q1)
    colnames(out) <- paste0("train_cov", seq_len(Q1))
    out
  } else {
    NULL
  }

  expression_cell_1 <- stats::rnorm(
    N,
    mean = 2 * x_train[, 1],
    sd = 1
  )
  expression_cell_2 <- stats::rnorm(N, mean = 0, sd = 1)

  y_train <- (
    expression_cell_1 * pi_train +
      expression_cell_2 * (1 - pi_train)
  )

  # Application data
  x_new <- matrix(
    sample(0:2, M * P, replace = TRUE),
    nrow = M,
    ncol = P
  )
  colnames(x_new) <- snp_names

  true_cell_1 <- stats::rnorm(
    M,
    mean = 2 * x_new[, 1],
    sd = 1
  )
  true_cell_2 <- stats::rnorm(M, mean = 0, sd = 1)

  cov_new <- if (Q2 > 0L) {
    out <- matrix(stats::rnorm(M * Q2), nrow = M, ncol = Q2)
    colnames(out) <- paste0("new_cov", seq_len(Q2))
    out
  } else {
    NULL
  }

  outcome <- stats::rnorm(
    M,
    mean = 2 * true_cell_1 + 0.5 * true_cell_2,
    sd = 1
  )

  list(
    gene_name = gene_name,
    pi_train = pi_train,
    x_train = x_train,
    x_annotation = x_annotation,
    cov_train = cov_train,
    y_train = y_train,
    x_new = x_new,
    cov_new = cov_new,
    outcome = outcome
  )
}


make_prediction_inputs <- function(n = 80L, p = 4L, seed = 2026L) {
  stopifnot(n >= 10L, p >= 2L)

  set.seed(seed)

  new_x <- matrix(
    sample(0:2, n * p, replace = TRUE),
    nrow = n,
    ncol = p
  )

  weight <- cbind(
    weight_cell_1 = seq(0.1, 0.1 * p, by = 0.1),
    weight_cell_2 = seq(-0.05, -0.05 * p, by = -0.05)
  )

  list(weight = weight, new_x = new_x)
}


make_association_inputs <- function(
    n = 300L,
    seed = 2026L
) {
  set.seed(seed)

  new_y <- cbind(
    cell_1 = stats::rnorm(n),
    cell_2 = stats::rnorm(n)
  )

  new_cov <- cbind(
    age = stats::rnorm(n),
    score = stats::rnorm(n)
  )

  outcome <- (
    1.5 * new_y[, 1] -
      0.75 * new_y[, 2] +
      0.25 * new_cov[, 1] +
      stats::rnorm(n, sd = 0.75)
  )

  list(
    new_y = new_y,
    new_cov = new_cov,
    outcome = outcome
  )
}
