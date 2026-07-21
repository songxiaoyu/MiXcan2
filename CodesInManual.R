rm(list=ls())
library(MiXcan2)
library(doParallel)
library(tidyverse)
nCores=detectCores()-1; registerDoParallel(nCores) # use parallel computing for speed,
# Pseudo data creation
gene_name <- "BRCA1"

set.seed(123)

P <- 10   # number of genetic predictors in both training and application data
N <- 100  # number of training data
Q1 <- 3    # number of covariates in training data
M <- 1000 # number of application data
Q2 <- 2    # number of covariates in application data

# --- Simulation training data 

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

# Covariate matrix (e.g. age, principal components)
Cov_test <- matrix(rnorm(N * Q1), nrow = N, ncol = Q1)
colnames(Cov_test) <- c("age", "PC1", "PC2")

# Expression levels
Y_test1 <- rnorm(N, mean = 2*X_test[,1], sd = 1)
Y_test2 <- rnorm(N, mean = 0, sd = 1)
Y_test <- Y_test1*Pi_test+Y_test2*(1-Pi_test)

# ------ Simulate application data

# Genotype matrix; column names are the same as training data
XM_test <- matrix(sample(0:2, M * P, replace = TRUE), nrow = M, ncol = P)
colnames(XM_test) <- paste0("V", sample(40:100, P))

# Expression levels (unobservable in real data)
YM1_test <- rnorm(M, mean = 2*XM_test[,1], sd = 1)
YM2_test <- rnorm(M, mean = 0, sd = 1)

# Covariate matrix (can be different from training data)
CovM_test <- matrix(rnorm(M * Q2), nrow = M, ncol = Q2)
colnames(CovM_test) <- c("age", "score")

# phenotype associated with expression in both cell types
DM_test <- rnorm(M, mean = 2*YM1_test+0.5*YM2_test, sd = 1)

# We can observe the following variables in the training data:
Pi_test[1:6]
X_test[1:6,]
X_rows_test[1:6,]
Y_test[1:6]
Cov_test[1:6,]

# We can observe the following variables in the application data:
XM_test[1:6,]
DM_test[1:6]
CovM_test[1:6,]


# train 
ensbl <- MiXcan2_ensemble(y = Y_test, x = X_test, cov = Cov_test,
                          pi = Pi_test, xNameMatrix = X_rows_test,
                          yName = gene_name, B = 9, seed = 123)
# show results 
ensbl$ensemble_intecept
ensbl$ensemble_weight

all_weights_df <- ensbl$all_weights
head(all_weights_df)

ensbl_summary_by_type =ensbl$ensemble_summary_by_type
ensbl_summary_by_type

ensbl_all_summary =ensbl$all_summary 
ensbl_all_summary

# application to the new GWAS data
selected_X_idx=match(ensbl$ensemble_weight$varID, X_rows_test$varID)
selected_X_idx
pred=MiXcan2_prediction(weight=ensbl$ensemble_weight[selected_X_idx,c("weight_cell_1", "weight_cell_2")], 
                        new_x=as.matrix(XM_test[,selected_X_idx]))
pred[1:6,]
MiXcan2_association(new_y=pred, new_cov=CovM_test, 
                    new_outcome=DM_test, family = gaussian) 
  


