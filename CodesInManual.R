rm(list=ls())
library(MiXcan2)
library(doParallel)
library(tidyverse)
nCores=detectCores()-1; registerDoParallel(nCores) # use parallel computing for speed,
# Pseudo data creation
gene_name <- "BRCA1"

set.seed(123)
N <- 100  # number of training samples
M <- 1000 # number of GWAS samples
P <- 10   # number of genetic predictors
Q <- 3    # number of covariates

# Covariate matrix (e.g. age, principal components)
Cov_test <- matrix(rnorm(N * Q), nrow = N, ncol = Q)
colnames(Cov_test) <- c("age", "PC1", "PC2")

# Cell-type fraction estimates (bounded between 0 and 1)
Pi_test <- runif(N, min = 0.1, max = 0.9)

# Genotype matrix with V-prefixed column names
X_test <- matrix(sample(0:2, N * P, replace = TRUE), nrow = N, ncol = P)
colnames(X_test) <- paste0("V", sample(40:100, P))

# SNP annotation matrix
X_rows_test <- data.frame(
  varID    = paste0("chr21_", sample(1e6:5e6, P), "_",
                    sample(c("A","T","C","G"), P, replace = TRUE), "_",
                    sample(c("A","T","C","G"), P, replace = TRUE), "_b38"),
  position = sample(1e6:5e6, P),
  rsid     = paste0("rs", sample(1000:9999, P)),
  ref      = sample(c("A","T","C","G"), P, replace = TRUE),
  eff      = sample(c("A","T","C","G"), P, replace = TRUE)
)

# Expression levels
Y_test <- rnorm(N, mean = 0, sd = 1)


X_test[1:6,]
X_rows_test[1:6,]
Cov_test[1:6,]
Pi_test[1:6]
Y_test[1:6]


ensbl <- MiXcan2_ensemble(y = Y_test, x = X_test, cov = Cov_test,
                          pi = Pi_test, xNameMatrix = X_rows_test,
                          yName = gene_name, B = 5, seed = 123)
ensbl$ensemble_intecept

all_weights_df <- ensbl$all_weights
head(all_weights_df)

ensbl_summary_by_type =ensbl$ensemble_summary_by_type
ensbl_summary_by_type

ensbl_all_summary =ensbl$all_summary 
ensbl_all_summary




