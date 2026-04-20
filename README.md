
<!-- README.md is generated from README.Rmd. Please edit README.Rmd file -->

# `MiXcan: Statistical Framework for Cell-type-Aware Transcriptome-Wide Association Studies with Bulk Tissue Data`

## Introduction to **MiXcan**

**Goal:**

- Constructs cell-type-level prediction models for genetically regulated
  expression (GReX);

- Predicts cell-type-level GReX in new genotype data; and

- Performs cell-type-aware TWAS.

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
  estimates (e.g. from existing methods, such as ESTIMATE, CIBERSORT,
  xCell).

- Association Analysis: genotype, covariates and phenotype data (same as
  in PrediXcan).

**Output:**

- Prediction model construction: Cell-type-specific or nonspecific
  prediction weights for different genes.

- Association Analysis: Tissue-level association p-values and
  cell-type-level association summaries including estimates, standard
  error and p-values.

A full description of the method can be found in our
[paper](https://www.biorxiv.org/content/10.1101/2022.03.15.484509v1.abstract).

``` r
knitr::opts_chunk$set(echo = TRUE)
library(MiXcan2) 
```

## Installation

### Hardware Requirements

The MiXcan package requires only a standard computer with a reasonable
RAM/CPU to support the operations. The minimal RAM recommended by
authors is 2 GB. Authors used a computer with the following
specification:

RAM: 32 GB

CPU: 2.3 GHz 8-core

### Software Requirements

The github package is developed on macOS operating systems. The MiXcan
pacakge is developed under an open source software R (version 4.1.2).
Different versions of R can be downloaded at
<https://cran.r-project.org/>.

### Package Installation

With R, users can install the MiXcan package directly from GitHub with
[devtools](https://github.com/hadley/devtools):

``` r
install.packages("devtools")
devtools::install_github("songxiaoyu/MiXcan2/Package")
```

The typical install time of the package is less than 5 minutes.

## Example of use

Below demonstrates the MiXcan analysis pipeline on a single peusdo gene.
In reality, multiple genes can be analyzed in parallel. Holdout set can
be pre-excluded to allow model training on the remaining samples.

### Data

Create peusdo data for demonstration:

```r
# Pseudo data creation
gene_names <- "BRCA1"

set.seed(123)
N <- 100  # number of samples
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
```

### MiXcan analysis pipeline

``` r
library(doParallel)
```

    ## Loading required package: foreach

    ## Loading required package: iterators

    ## Loading required package: parallel

``` r
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.4
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.4.4     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.3     ✔ tidyr     1.3.0
    ## ✔ purrr     1.0.2

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ purrr::accumulate() masks foreach::accumulate()
    ## ✖ dplyr::filter()     masks stats::filter()
    ## ✖ dplyr::lag()        masks stats::lag()
    ## ✖ purrr::when()       masks foreach::when()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(rlist)
nCores=detectCores()-1; registerDoParallel(nCores) # use parallel computing for speed, but leave 1 core out for other activities. 
```

Step 2: Estimating cell-type-specific (and nonspecific) GReX prediction weights of a gene using the MiXcan function
```
# X_test
#      V70 V53 V59 V79 V63 V76 V80 V68 V49 V88
# [1,]   1   1   0   2   0   1   2   0   2   1
# [2,]   0   2   2   0   2   2   2   2   1   0
# [3,]   0   0   2   2   0   2   1   0   2   2
# [4,]   2   1   2   1   2   0   1   0   2   2
# [5,]   0   0   1   0   2   1   1   1   2   0
# [6,]   2   1   2   0   1   2   2   2   1   0

# X_rows_test
#             xNameMatrix position   rsid ref eff
# 1 chr21_4380166_A_A_b38  4598800 rs2457   G   C
# 2 chr21_2424043_A_C_b38  4171961 rs5424   A   T
# 3 chr21_1669969_T_C_b38  4261569 rs9509   C   G
# 4 chr21_1455554_T_G_b38  2729790 rs7777   C   A
# 5 chr21_4849901_G_A_b38  2055678 rs9293   C   C
# 6 chr21_2125566_T_T_b38  3991202 rs6735   C   A

# Cov_test
#             age         PC1        PC2
# [1,] -0.5604757 -0.7104066  2.1988103
# [2,] -0.2301775  0.2568837  1.3124130
# [3,]  1.5587083 -0.2466919 -0.2651451
# [4,]  0.0705084 -0.3475426  0.5431941
# [5,]  0.1292877 -0.9516186 -0.4143399
# [6,]  1.7150650 -0.0450277 -0.4762469

# Pi_test
# [1] 0.2897838 0.6491923 0.2806547 0.3547957 0.2391871 0.7411437

# Y_test (numeric vector of length N)
# [1]  1.7791027  0.5351380 -0.3719449 -1.0255422 -0.5824017  0.3428884
```


``` r
set.seed(123)
ensbl <- MiXcan2_ensemble(y = Y_test, x = X_test, cov = Cov_test,
                              pi = Pi_test, xNameMatrix = X_rows_test,
                              yName = gene_name, B = 3, seed = 123)
ensbl$ensemble_intecept
```

    #              intercept_cell_1 intercept_cell_2
    # NonSpecific   -0.0246929       -0.0246929
    # NoPredictor   -0.0869106       -0.0869106

Step 3: Extracting the weights and model summaries from the MiXcan output.

Raw SNP weights from every individual ensemble model (across all B iterations).
No averaging applied across B model interations and retaining all weights including 0.
``` r
all_weights_df <- ensbl$all_weights
head(all_weights_df)

```
    # all_weights_df
    # ID yName xNameMatrix         position rsid ref eff weight_cell_1 weight_cell_2 type
    # 1 BRCA1 chr21_4380166_A_A_b38 4598800 rs2457 G C           0           0       NoPredictor
    # 1 BRCA1 chr21_2424043_A_C_b38 4171961 rs5424 A T           0           0       NoPredictor
    # 1 BRCA1 chr21_1669969_T_C_b38 4261569 rs9509 C G           0           0       NoPredictor


Overall averaged model performance metrics (e.g. cross-validated R², number of SNPs)
Results are summarized across all B models into a single row per gene
``` r
## ensemble_summary
ensbl_summary = ensbl$ensemble_summary
ensbl_summary
```
    ##Gene n_snp_input n_snp_model in.sample.cor in.sample.R2     cv.cor    cv.R2        CTS        NS        NP
    ##1 BRCA1      10      0.6666667     0.1401488    0.1401488    -0.2263911   -0.2263911   0       0.6666667 0.3333333

Overall averaged model performance metrics (e.g. cross-validated R², number of SNPs)
Summarized separately for each model type (Cell-Type-Specific, Non-Specific, NoPredictor)
``` r
ensbl_summary_by_type =ensbl$ensemble_summary_by_type
ensbl_summary_by_type
```

    ##   Gene  model_type  n_snp_input n_snp_model in.sample.cor in.sample.R2 cv.cor
    ## 1 BRCA1 NoPredictor          10           0         0            0      0    
    ## 2 BRCA1 NonSpecific          10           1         0.210        0.210 -0.340

Model performance metrics from each ensemble model (across B interations). 
Reults are prior to averaging, one row per model iteration.
``` r
ensbl_all_summary =ensbl$all_summary 
ensbl_all_summary
```

    # yName n_snp_input n_snp_model   model_type   in.sample.cor in.sample.R2   cv.cor     cv.R2
    # BRCA1   10           0          NoPredictor   0.0000000     0.0000000   0.0000000   0.0000000
    # BRCA1   10           1          NonSpecific   0.1986117     0.1986117   -0.3679533   -0.3679533
    # BRCA1   10           1          NonSpecific   0.2218347     0.2218347   -0.3112200 -  0.3112200
