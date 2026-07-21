
<!-- README.md is NOT generated from README.Rmd. Please edit README.Rmd file -->

# `Cell-type-aware transcriptome-wide association study of mammographic density phenotypes`

A full description of the manuscript can be found 
[here](https://www.medrxiv.org/content/10.1101/2025.11.25.25341027v1).



## 1.Introduction

Analyses method **MiXcan2:**

- Constructs cell-type-level prediction models for genetically regulated
  expression (GReX) via an ensemble strategy;
- Predicts cell-type-level GReX in new genotype data; and
- Performs cell-type-aware TWAS.
- A full description of the method can be found 
[here](https://www.nature.com/articles/s41467-023-35888-4).

Analyses method  **S-MiXcan**

- With cell-type-level prediction models for GReX, infers cell-type-level GReX associations directly from GWAS summary statistics without requiring individual-level genotype data;
- A full description of the method can be found 
[here](https://www.medrxiv.org/content/10.64898/2026.03.22.26349035v1)

**Advantages over tissue-level TWAS:**

- Improves GReX prediction accuracy;
- Boosts the study power, especially for genes that function in minor
  cell types or have different association directions in different cell
  types;
- Sheds light on the responsible cell type(s) of associations.

**Disadvantages over tissue-level TWAS:**

- Requires prior knowledge on disease-critical cell types and their
  proportions in tissue;
- Has more model parameters;
- May be less powerful than tissue-level TWAS for genes that have
  similar disease associations in different cell types or function in
  major cell types.

**Input:**

- Prediction model construction: genotype, covariates, and gene
  expression data (same as in PrediXcan) + cell-type composition
  estimates (e.g. from existing methods, such as BayesDebulk).

- Association Analysis: individual level data (e.g. genotype, covariates and phenotype data; same as
  in PrediXcan) or GWAS summary statistics.

**Output:**

- Prediction model construction: Cell-type-specific or nonspecific
  prediction weights for different genes.
- Association Analysis: Tissue-level association p-values and
  cell-type-level association summaries including estimates, standard
  error and p-values.


## 2.Installation

#### Hardware Requirements

The MiXcan package requires only a standard computer with a reasonable
RAM/CPU to support the operations. The minimal RAM recommended by
authors is 2 GB. Authors used a computer with the following
specification:

RAM: 32 GB
CPU: 2.3 GHz 8-core

#### Software Requirements

The github package is developed on macOS operating systems. The MiXcan
pacakge is developed under an open source software R (version 4.1.2).
Different versions of R can be downloaded at
<https://cran.r-project.org/>.

#### Package Installation

With R, users can install the MiXcan package directly from GitHub with
[devtools](https://github.com/hadley/devtools):

``` r
knitr::opts_chunk$set(echo = TRUE)
library(MiXcan2) 
```
``` r
# installation of MiXcan2
install.packages("devtools")
devtools::install_github("songxiaoyu/MiXcan2/Package")
```

For install for S-MiXcan, please refer its own [GitHub](https://github.com/songxiaoyu/SMiXcan) page.

The typical install time of the package is less than 5 minutes.


## 3. Reproduce the Publication Results

The entire reproducing pipeline to generate results from raw data is available [here](https://github.com/songxiaoyu/MiXcan2/blob/main/REPRODUCING_RESULTS.md).
Please note that the genomics data are not accessible due to participant privacy and data governance 
requirements, and for those parts we only provide codes and data access instructions (in manuscript).
For validation studies and breast cancer enrichment analyses, where raw data is not needed, we provided
end-to-end pipeline for analysis. 

## 4. Toy Example to Illustrate the Usage

As the input data used in our study is not publicly available, the
example below uses simulated pseudo data to demonstrate the MiXcan
analysis pipeline on a single gene. In reality, multiple genes can be
analyzed in parallel, and a holdout set can be pre-excluded to allow
model training on the remaining samples.

#### Simuate Data

``` r
## Pseudo data creation
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
```
#### Simulated data that can be obeserved

We can observe the following variables in the training data:
``` r
> Pi_test[1:6]
[1] 0.3300620 0.7306441 0.4271815 0.8064139 0.8523738 0.1364452
> X_test[1:6,]
     V94 V78 V76 V47 V86 V53 V87 V40 V62 V65
[1,]   0   0   1   0   2   1   0   0   1   2
[2,]   1   1   0   1   2   1   1   2   1   1
[3,]   0   0   0   0   2   1   2   0   1   1
[4,]   1   2   0   0   1   1   2   2   1   1
[5,]   0   0   1   0   0   2   1   2   1   2
[6,]   2   0   0   1   0   0   0   2   0   0
> X_rows_test[1:6,]
                  varID position   rsid ref eff
1 chr21_3958870_T_T_b38  4397319 rs4711   A   C
2 chr21_2968830_C_A_b38  2637837 rs6551   T   A
3 chr21_2359184_C_G_b38  3015021 rs7347   G   T
4 chr21_1426028_T_A_b38  2611015 rs7188   G   G
5 chr21_4053899_A_T_b38  3753398 rs2362   A   G
6 chr21_3078381_G_C_b38  4407639 rs9549   G   G
> Y_test[1:6]
[1]  1.2214868  0.9029712 -1.0829126  2.2228998 -0.2397309  0.8757725
> Cov_test[1:6,]
            age          PC1         PC2
[1,] -1.3501461  0.929726653  0.71756885
[2,] -0.2645207  0.312910943 -0.61896534
[3,] -0.4718988 -0.000204297 -0.06492498
[4,] -0.5691862 -0.250647707  0.75658230
[5,]  2.3327194  1.520275176 -1.01422951
[6,] -0.4890443 -2.093256519 -0.69207082
```
We can observe the following variables in the application data:
``` r
> XM_test[1:6,]
     V96 V87 V97 V68 V74 V77 V45 V40 V73 V94
[1,]   0   1   0   2   2   2   1   2   2   1
[2,]   1   2   2   2   1   0   0   1   2   2
[3,]   0   0   1   0   2   0   1   2   0   1
[4,]   1   1   2   1   2   1   2   1   0   1
[5,]   1   2   2   1   1   2   2   0   0   0
[6,]   1   2   0   1   1   2   0   0   2   1
> DM_test[1:6]
[1] 0.4074296 4.5332502 1.1829998 4.0772315 3.1780357 1.5153982
> CovM_test[1:6,]
             age      score
[1,] -0.36602729 -0.5751098
[2,] -0.12025695 -1.7366401
[3,]  1.58980932 -0.2018158
[4,] -1.58950575 -0.3643292
[5,]  0.25009968  0.6180302
[6,] -0.09430046  0.8100630

```


#### Traing MiXcan2 prediction model
Estimating cell-type-specific (and nonspecific) GReX prediction
weights of a gene using the MiXcan function.

``` r
library(doParallel)
library(tidyverse)
library(rlist)
nCores=detectCores()-1; registerDoParallel(nCores) # use parallel computing for speed, but leave 1 core out for other activities. 
```




``` r
ensbl <- MiXcan2_ensemble(y = Y_test, x = X_test, cov = Cov_test,
+                           pi = Pi_test, xNameMatrix = X_rows_test,
+                           yName = gene_name, B = 9, seed = 123)

# show results - assembled intercepts
ensbl$ensemble_intercept
                 intercept_cell_1 intercept_cell_2
CellTypeSpecific         1.096866       -0.0784024
# show results - assembled weights (will be used in the downstream analyses)
ensbl$ensemble_weight
# A tibble: 8 × 9
# Groups:   yName, varID, position, rsid, ref, eff [8]
  yName varID    position rsid  ref   eff   type  weight_cell_1 weight_cell_2
  <chr> <chr>       <int> <chr> <chr> <chr> <chr>         <dbl>         <dbl>
1 BRCA1 chr21_1…  2611015 rs71… G     G     Cell…        0.0157      -0.00525
2 BRCA1 chr21_2…  3162112 rs19… T     A     Cell…        0.0537       0.0537 
3 BRCA1 chr21_2…  1282765 rs18… T     C     Cell…       -0.0616      -0.0616 
4 BRCA1 chr21_2…  3015021 rs73… G     T     Cell…       -0.0136      -0.0136 
5 BRCA1 chr21_2…  2637837 rs65… T     A     Cell…       -0.0839      -0.0839 
6 BRCA1 chr21_3…  4407639 rs95… G     G     Cell…       -0.0582       0.0516 
7 BRCA1 chr21_3…  4397319 rs47… A     C     Cell…        1.23         0.231  
8 BRCA1 chr21_4…  3753398 rs23… A     G     Cell…       -0.103       -0.103  

# Weights for all the B boostrap samples 
all_weights_df <- ensbl$all_weights
head(all_weights_df)
  ID yName                 varID position   rsid ref eff weight_cell_1
1  1 BRCA1 chr21_3958870_T_T_b38  4397319 rs4711   A   C    1.33605586
2  1 BRCA1 chr21_2968830_C_A_b38  2637837 rs6551   T   A   -0.07936197
3  1 BRCA1 chr21_2359184_C_G_b38  3015021 rs7347   G   T   -0.05338579
4  1 BRCA1 chr21_1426028_T_A_b38  2611015 rs7188   G   G    0.00000000
5  1 BRCA1 chr21_4053899_A_T_b38  3753398 rs2362   A   G   -0.19830951
6  1 BRCA1 chr21_3078381_G_C_b38  4407639 rs9549   G   G   -0.16260200
  weight_cell_2             type
1    0.12890964 CellTypeSpecific
2   -0.07936197 CellTypeSpecific
3   -0.05338579 CellTypeSpecific
4    0.00000000 CellTypeSpecific
5   -0.19830951 CellTypeSpecific
6    0.10326876 CellTypeSpecific

# Model summary for the selected model type 
ensbl_summary_by_type =ensbl$ensemble_summary_by_type
> ensbl_summary_by_type
# A tibble: 1 × 8
  Gene  model_type  n_snp_input n_snp_model in.sample.cor in.sample.R2 cv.cor
  <chr> <chr>             <dbl>       <dbl>         <dbl>        <dbl>  <dbl>
1 BRCA1 CellTypeSp…          10        6.33         0.804        0.626  0.732
# ℹ 1 more variable: cv.R2 <dbl>

# Model summary for each of the B boostrap samples
ensbl_all_summary =ensbl$all_summary 
ensbl_all_summary
  yName n_snp_input n_snp_model       model_type in.sample.cor in.sample.R2
1 BRCA1          10           7 CellTypeSpecific     0.8271396    0.6683251
2 BRCA1          10           7 CellTypeSpecific     0.8239286    0.6597136
3 BRCA1          10           6 CellTypeSpecific     0.7969319    0.6131068
4 BRCA1          10           5 CellTypeSpecific     0.7790375    0.5812645
5 BRCA1          10           5 CellTypeSpecific     0.7942746    0.6074591
6 BRCA1          10           6 CellTypeSpecific     0.7918577    0.6092398
7 BRCA1          10           6 CellTypeSpecific     0.8136668    0.6433772
8 BRCA1          10           8 CellTypeSpecific     0.7965697    0.6181685
9 BRCA1          10           7 CellTypeSpecific     0.8091436    0.6376223
     cv.cor     cv.R2
1 0.7517858 0.5616249
2 0.7527023 0.5597877
3 0.7391163 0.5346809
4 0.6885775 0.4686230
5 0.7374831 0.5328995
6 0.7084881 0.4988803
7 0.7639915 0.5731816
8 0.7229446 0.5189709
9 0.7245171 0.5178526
```

#### Apply prediction weights to the new genomics data
``` r
# select SNP with nonzero weights
selected_X_idx=match(ensbl$ensemble_weight$varID, X_rows_test$varID)
> selected_X_idx
[1] 4 8 7 3 2 6 1 5

# apply the weights to new genotype data
pred=MiXcan2_prediction(weight=ensbl$ensemble_weight[selected_X_idx,c("weight_cell_1", "weight_cell_2")], 
                       new_x=as.matrix(XM_test[,selected_X_idx]))
pred[1:6,]
          cell_1      cell_2
[1,]  0.76815421 -0.01360668
[2,] -0.21432724 -0.23523150
[3,]  0.79654834 -0.20497674
[4,]  2.06816667  0.15409435
[5,]  2.25082166  0.44663144
[6,] -0.09084614  0.10801379
```
#### Associate predicted cell-type-level expression with phenotype
``` r
MiXcan2_association(new_y=pred, new_cov=CovM_test, 
                    new_outcome=DM_test, family = gaussian) 
                    
  cell1_est  cell1_se    cell1_p cell2_est cell2_se   cell2_p p_combined
1 0.5258613 0.2478619 0.03411882 -1.489066 1.091636 0.1728543 0.05750228
```