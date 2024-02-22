rm(list=ls())
library(MiXcan2)
library(doParallel)
library(tidyverse)
nCores=detectCores()-1; registerDoParallel(nCores) # use parallel computing for speed,

# make up a dataset------
# set.seed(123)
# n=200
# pi=rbeta(n, 2, 3)
# x=matrix(rbinom(n*10, 2, 0.3), ncol=10)
# y=x[,1]*pi+rnorm(n, sd=0.2)
# cov=NULL
# 
# var(2*x[,1]*pi)/var(y) # real 2
# 
# a=lm(y~x[,1]+pi+x[,1]*pi)
# var(predict(a))/var(y) # r2
# a=lm(y~x[,1])
# var(predict(a))/var(y) # r2
# a=lm(y~pi)
# var(predict(a))/var(y) # r2
# 
# library(MiXcan)
data(example_data)
set.seed(123)
pi_estimation_result <- pi_estimation(expression_matrix = GTEx_epithelial_genes,
                                      prior = GTEx_prior, n_iteration = 5)
pi=pi_estimation_result$mean_trim_0.05
x=x_example
y=y_example
cov=cov_example

# MiXcan2  -------
set.seed(222)
foldid <- sample(1:10, length(y), replace=T)

model <- MiXcan2_model(y=y, x=x, cov = cov,
                        pi= pi,
                        foldid = foldid, yName="Gene1",
                       xNameMatrix = paste0("X", 1:ncol(x)))
model$beta.SNP.cell1
model$beta.SNP.cell2

# Extract Info
MiXcan_weight_result <- MiXcan2_extract_weight(model = model)
MiXcan_weight_result

MiXcan_summary_result <- MiXcan2_extract_summary(model=model)
MiXcan_summary_result

# Refit
MiXcan_refit <- MiXcan2_refit(model = model)
MiXcan_refit$weight
MiXcan_refit$summary
# MiXcan2 Ensemble -------
ensemble=MiXcan2_ensemble(y=y, x=x, cov=cov, pi=pi, 
                 yName="Gene1", B=10, seed=123) 
ensemble$ensemble_summary_by_type
ensemble$ensemble_summary
ensemble$ensemble_weight
ensemble$all_weights
