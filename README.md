
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

``` r
set.seed(123)
n=200
pi=rbeta(n, 2, 3)
x=matrix(rbinom(n*10, 2, 0.3), ncol=10)
y=x[,1]*pi+rnorm(n, sd=0.2)
cov=NULL
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

Step 2: Estimating cell-type-specific (and nonspecific) GReX prediction
weights of a gene using the MiXcan function

``` r
set.seed(111)
foldid <- sample(1:10, length(y), replace=T)

model <- MiXcan2_model(y=y, x=x, cov = cov,
                        pi= pi,
                        foldid = foldid, yName="Gene1",
                       xNameMatrix = paste0("X", 1:10))
model$beta.SNP.cell1
```

    ##    xNameMatrix    weight
    ## 1           X1 0.4985004
    ## 2           X2 0.0000000
    ## 3           X3 0.0000000
    ## 4           X4 0.0000000
    ## 5           X5 0.0000000
    ## 6           X6 0.0000000
    ## 7           X7 0.0000000
    ## 8           X8 0.0000000
    ## 9           X9 0.0000000
    ## 10         X10 0.0000000

``` r
model$beta.SNP.cell2
```

    ##    xNameMatrix   weight
    ## 1           X1 0.181293
    ## 2           X2 0.000000
    ## 3           X3 0.000000
    ## 4           X4 0.000000
    ## 5           X5 0.000000
    ## 6           X6 0.000000
    ## 7           X7 0.000000
    ## 8           X8 0.000000
    ## 9           X9 0.000000
    ## 10         X10 0.000000

Step 3: Extracting the weights and model summaries from the MiXcan
output.

``` r
MiXcan_weight_result <- MiXcan2_extract_weight(model = model)
```

    ## Joining with `by = join_by(xNameMatrix)`

``` r
MiXcan_weight_result
```

    ##   yName xNameMatrix weight_cell_1 weight_cell_2             type
    ## 1 Gene1          X1     0.4985004      0.181293 CellTypeSpecific

``` r
MiXcan_summary_result <- MiXcan2_extract_summary(x=x, y=y, cov=cov,
                                                pi=pi, model=model)
```

    ## Joining with `by = join_by(xNameMatrix)`

``` r
MiXcan_summary_result
```

    ##   yName n_snp_input n_snp_model       model_type     cv_r2 in_sample_r2
    ## 1 Gene1          10           1 CellTypeSpecific 0.5500884    0.5734584

Note, the MiXcan estimated models are from penalized regression
(elastic-net), which shrinks the effect size towards zero. If users are
interested to use less penalized weights for the MiXcan selected SNPs,
they can employ the following function:

``` r
MiXcan_refit <- MiXcan2_refit(model = model,
                                           y=y,
                                           x=x, cov = cov,
                                           pi= pi)
```

    ## Joining with `by = join_by(xNameMatrix)`
    ## Joining with `by = join_by(xNameMatrix)`
    ## Joining with `by = join_by(xNameMatrix)`

``` r
MiXcan_refit$weight
```

    ##   yName xNameMatrix weight_cell_1 weight_cell_2             type
    ## 1 Gene1          X1     0.9149063    0.03824874 CellTypeSpecific

``` r
MiXcan_refit$summary
```

    ##   yName n_snp_input n_snp_model       model_type     cv_r2 in_sample_r2
    ## 1 Gene1          10           1 CellTypeSpecific 0.5500884    0.5734584
    ##   cv_r2_refit in_sample_r2_refit
    ## 1   0.6145838           0.636116

One can run MiXcan for multiple times to ensemble models for downstream
anlaysis. Here are the codes:

``` r
ensemble=MiXcan2_ensemble(y=y, x=x, cov=cov, pi=pi, 
                 yName="Gene1", B=10, seed=123) 
  
ensemble$summary
```

    ##    Gene n_snp_input n_snp_model no_of_models     cv_r2 in_sample_r2 cv_r2_refit
    ## . Gene1          10         2.2           10 0.5630061    0.5846373   0.6228923
    ##   in_sample_r2_refit NP CTS NS
    ## .          0.6411247  0   1  0

``` r
ensemble$ensemble_weight
```

    ## # A tibble: 3 × 5
    ## # Groups:   yName, xNameMatrix [3]
    ##   yName xNameMatrix type             weight_cell_1 weight_cell_2
    ##   <chr> <chr>       <chr>                    <dbl>         <dbl>
    ## 1 Gene1 SNP1        CellTypeSpecific        0.909         0.0442
    ## 2 Gene1 SNP4        CellTypeSpecific        0.0267        0.0267
    ## 3 Gene1 SNP6        CellTypeSpecific        0.0129        0.0129

``` r
ensemble$all_weights
```

    ##     ID yName xNameMatrix weight_cell_1 weight_cell_2             type
    ## 1    1 Gene1        SNP1    0.91042963    0.04152770 CellTypeSpecific
    ## 2    1 Gene1        SNP2    0.00000000    0.00000000 CellTypeSpecific
    ## 3    1 Gene1        SNP3    0.00000000    0.00000000 CellTypeSpecific
    ## 4    1 Gene1        SNP4    0.03301747    0.03301747 CellTypeSpecific
    ## 5    1 Gene1        SNP5    0.00000000    0.00000000 CellTypeSpecific
    ## 6    1 Gene1        SNP6    0.00000000    0.00000000 CellTypeSpecific
    ## 7    1 Gene1        SNP7    0.00000000    0.00000000 CellTypeSpecific
    ## 8    1 Gene1        SNP8    0.00000000    0.00000000 CellTypeSpecific
    ## 9    1 Gene1        SNP9    0.00000000    0.00000000 CellTypeSpecific
    ## 10   1 Gene1       SNP10    0.00000000    0.00000000 CellTypeSpecific
    ## 11   2 Gene1        SNP1    0.91042963    0.04152770 CellTypeSpecific
    ## 12   2 Gene1        SNP2    0.00000000    0.00000000 CellTypeSpecific
    ## 13   2 Gene1        SNP3    0.00000000    0.00000000 CellTypeSpecific
    ## 14   2 Gene1        SNP4    0.03301747    0.03301747 CellTypeSpecific
    ## 15   2 Gene1        SNP5    0.00000000    0.00000000 CellTypeSpecific
    ## 16   2 Gene1        SNP6    0.00000000    0.00000000 CellTypeSpecific
    ## 17   2 Gene1        SNP7    0.00000000    0.00000000 CellTypeSpecific
    ## 18   2 Gene1        SNP8    0.00000000    0.00000000 CellTypeSpecific
    ## 19   2 Gene1        SNP9    0.00000000    0.00000000 CellTypeSpecific
    ## 20   2 Gene1       SNP10    0.00000000    0.00000000 CellTypeSpecific
    ## 21   3 Gene1        SNP1    0.91490630    0.03824874 CellTypeSpecific
    ## 22   3 Gene1        SNP2    0.00000000    0.00000000 CellTypeSpecific
    ## 23   3 Gene1        SNP3    0.00000000    0.00000000 CellTypeSpecific
    ## 24   3 Gene1        SNP4    0.00000000    0.00000000 CellTypeSpecific
    ## 25   3 Gene1        SNP5    0.00000000    0.00000000 CellTypeSpecific
    ## 26   3 Gene1        SNP6    0.00000000    0.00000000 CellTypeSpecific
    ## 27   3 Gene1        SNP7    0.00000000    0.00000000 CellTypeSpecific
    ## 28   3 Gene1        SNP8    0.00000000    0.00000000 CellTypeSpecific
    ## 29   3 Gene1        SNP9    0.00000000    0.00000000 CellTypeSpecific
    ## 30   3 Gene1       SNP10    0.00000000    0.00000000 CellTypeSpecific
    ## 31   4 Gene1        SNP1    0.91042963    0.04152770 CellTypeSpecific
    ## 32   4 Gene1        SNP2    0.00000000    0.00000000 CellTypeSpecific
    ## 33   4 Gene1        SNP3    0.00000000    0.00000000 CellTypeSpecific
    ## 34   4 Gene1        SNP4    0.03301747    0.03301747 CellTypeSpecific
    ## 35   4 Gene1        SNP5    0.00000000    0.00000000 CellTypeSpecific
    ## 36   4 Gene1        SNP6    0.00000000    0.00000000 CellTypeSpecific
    ## 37   4 Gene1        SNP7    0.00000000    0.00000000 CellTypeSpecific
    ## 38   4 Gene1        SNP8    0.00000000    0.00000000 CellTypeSpecific
    ## 39   4 Gene1        SNP9    0.00000000    0.00000000 CellTypeSpecific
    ## 40   4 Gene1       SNP10    0.00000000    0.00000000 CellTypeSpecific
    ## 41   5 Gene1        SNP1    0.90482455    0.04996658 CellTypeSpecific
    ## 42   5 Gene1        SNP2    0.00000000    0.00000000 CellTypeSpecific
    ## 43   5 Gene1        SNP3    0.00000000    0.00000000 CellTypeSpecific
    ## 44   5 Gene1        SNP4    0.03382514    0.03382514 CellTypeSpecific
    ## 45   5 Gene1        SNP5    0.00000000    0.00000000 CellTypeSpecific
    ## 46   5 Gene1        SNP6    0.03230958    0.03230958 CellTypeSpecific
    ## 47   5 Gene1        SNP7    0.00000000    0.00000000 CellTypeSpecific
    ## 48   5 Gene1        SNP8    0.00000000    0.00000000 CellTypeSpecific
    ## 49   5 Gene1        SNP9    0.00000000    0.00000000 CellTypeSpecific
    ## 50   5 Gene1       SNP10    0.00000000    0.00000000 CellTypeSpecific
    ## 51   6 Gene1        SNP1    0.90482455    0.04996658 CellTypeSpecific
    ## 52   6 Gene1        SNP2    0.00000000    0.00000000 CellTypeSpecific
    ## 53   6 Gene1        SNP3    0.00000000    0.00000000 CellTypeSpecific
    ## 54   6 Gene1        SNP4    0.03382514    0.03382514 CellTypeSpecific
    ## 55   6 Gene1        SNP5    0.00000000    0.00000000 CellTypeSpecific
    ## 56   6 Gene1        SNP6    0.03230958    0.03230958 CellTypeSpecific
    ## 57   6 Gene1        SNP7    0.00000000    0.00000000 CellTypeSpecific
    ## 58   6 Gene1        SNP8    0.00000000    0.00000000 CellTypeSpecific
    ## 59   6 Gene1        SNP9    0.00000000    0.00000000 CellTypeSpecific
    ## 60   6 Gene1       SNP10    0.00000000    0.00000000 CellTypeSpecific
    ## 61   7 Gene1        SNP1    0.90482455    0.04996658 CellTypeSpecific
    ## 62   7 Gene1        SNP2    0.00000000    0.00000000 CellTypeSpecific
    ## 63   7 Gene1        SNP3    0.00000000    0.00000000 CellTypeSpecific
    ## 64   7 Gene1        SNP4    0.03382514    0.03382514 CellTypeSpecific
    ## 65   7 Gene1        SNP5    0.00000000    0.00000000 CellTypeSpecific
    ## 66   7 Gene1        SNP6    0.03230958    0.03230958 CellTypeSpecific
    ## 67   7 Gene1        SNP7    0.00000000    0.00000000 CellTypeSpecific
    ## 68   7 Gene1        SNP8    0.00000000    0.00000000 CellTypeSpecific
    ## 69   7 Gene1        SNP9    0.00000000    0.00000000 CellTypeSpecific
    ## 70   7 Gene1       SNP10    0.00000000    0.00000000 CellTypeSpecific
    ## 71   8 Gene1        SNP1    0.90482455    0.04996658 CellTypeSpecific
    ## 72   8 Gene1        SNP2    0.00000000    0.00000000 CellTypeSpecific
    ## 73   8 Gene1        SNP3    0.00000000    0.00000000 CellTypeSpecific
    ## 74   8 Gene1        SNP4    0.03382514    0.03382514 CellTypeSpecific
    ## 75   8 Gene1        SNP5    0.00000000    0.00000000 CellTypeSpecific
    ## 76   8 Gene1        SNP6    0.03230958    0.03230958 CellTypeSpecific
    ## 77   8 Gene1        SNP7    0.00000000    0.00000000 CellTypeSpecific
    ## 78   8 Gene1        SNP8    0.00000000    0.00000000 CellTypeSpecific
    ## 79   8 Gene1        SNP9    0.00000000    0.00000000 CellTypeSpecific
    ## 80   8 Gene1       SNP10    0.00000000    0.00000000 CellTypeSpecific
    ## 81   9 Gene1        SNP1    0.91490630    0.03824874 CellTypeSpecific
    ## 82   9 Gene1        SNP2    0.00000000    0.00000000 CellTypeSpecific
    ## 83   9 Gene1        SNP3    0.00000000    0.00000000 CellTypeSpecific
    ## 84   9 Gene1        SNP4    0.00000000    0.00000000 CellTypeSpecific
    ## 85   9 Gene1        SNP5    0.00000000    0.00000000 CellTypeSpecific
    ## 86   9 Gene1        SNP6    0.00000000    0.00000000 CellTypeSpecific
    ## 87   9 Gene1        SNP7    0.00000000    0.00000000 CellTypeSpecific
    ## 88   9 Gene1        SNP8    0.00000000    0.00000000 CellTypeSpecific
    ## 89   9 Gene1        SNP9    0.00000000    0.00000000 CellTypeSpecific
    ## 90   9 Gene1       SNP10    0.00000000    0.00000000 CellTypeSpecific
    ## 91  10 Gene1        SNP1    0.91042963    0.04152770 CellTypeSpecific
    ## 92  10 Gene1        SNP2    0.00000000    0.00000000 CellTypeSpecific
    ## 93  10 Gene1        SNP3    0.00000000    0.00000000 CellTypeSpecific
    ## 94  10 Gene1        SNP4    0.03301747    0.03301747 CellTypeSpecific
    ## 95  10 Gene1        SNP5    0.00000000    0.00000000 CellTypeSpecific
    ## 96  10 Gene1        SNP6    0.00000000    0.00000000 CellTypeSpecific
    ## 97  10 Gene1        SNP7    0.00000000    0.00000000 CellTypeSpecific
    ## 98  10 Gene1        SNP8    0.00000000    0.00000000 CellTypeSpecific
    ## 99  10 Gene1        SNP9    0.00000000    0.00000000 CellTypeSpecific
    ## 100 10 Gene1       SNP10    0.00000000    0.00000000 CellTypeSpecific
