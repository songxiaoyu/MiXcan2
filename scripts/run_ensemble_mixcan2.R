## =============================================================================
## OVERVIEW
## Run MiXcan2 ensemble models across a list of genes
## The output files of this script can be found here: https://zenodo.org/records/19075631
## This link contains the MiXcan2 models trained using GTEx v8 mammary tissue samples 
## from female subjects with European ancestry

## For each gene, fits an ensemble of MiXcan2 models (bootstrap-aggregated) to
## estimate cell-type-specific vs. non-specific gene expression models, using
## per-gene genotype/expression data plus BayesDeBulk cell-type proportion
## estimates.
##
## NOTE: The input data referenced by this script (per-gene genotype/expression
## files, covariates, and BayesDeBulk cell-type proportions) are NOT publicly
## available and are not included in this repository. Paths below are
## placeholders — supply your own data in the same structure to run this script.
##
## Required MiXcan2 functions (MiXcan2_ensemble, MiXcan2_extract_summary,
## MiXcan2_extract_weight, MiXcan2_model) are provided by the MiXcan2 package:
##   devtools::install_github("songxiaoyu/MiXcan2/Package")
## =============================================================================

## ---------------------------------------------------------------------------
## Input file structure (per gene)
## ---------------------------------------------------------------------------
##
## Each gene has its own directory containing:
##
##   Cov.csv      - Sample-level covariates (e.g. age, sex, ancestry PCs, batch effects)
##   X.csv        - Genotype matrix (N samples x P SNPs) for the gene
##   X_rows.csv   - SNP identifiers corresponding to the columns/features in X.csv
##   Y.csv        - Gene expression vector for the gene across N samples
##
## Cell-type proportions are NOT read from a per-gene file; instead they are
## loaded once from a BayesDeBulk output file (3 cell types: Epi, Fibr, Adipose)
## and must be aligned to the same samples/order as X.csv and Y.csv.
##
## MiXcan2_ensemble expects:
##   y   : expression vector for a single gene (N samples)
##   x   : genotype matrix for that gene (N x P SNPs)
##   cov : covariate matrix (N x Q covariates)
##   pi  : cell-type proportion estimates for the selected cell type (N samples)

## ---------------------------------------------------------------------------
## Setup
## ---------------------------------------------------------------------------

options(warn = -1)

library(doParallel)
library(dplyr)
library(readr)
library(tibble)
library(MiXcan2)   # devtools::install_github("songxiaoyu/MiXcan2/Package")

nCores <- detectCores() - 1
registerDoParallel(nCores)
print(nCores)

## ---------------------------------------------------------------------------
## Config -- EDIT THESE
## ---------------------------------------------------------------------------

input_dir       <- "path/to/per_gene_inputs"          # one subfolder per gene (see structure above)
output_dir      <- "path/to/results"
pi_file         <- "path/to/BayesDeBulk_pi_3ct.tsv"   # BayesDeBulk cell-type proportions (Epi, Fibr, Adipose)
gene_list_file  <- "ensembl_ids_intersecting.txt"      # newline-separated list of gene IDs to process

cell_type <- "Epi"     # which cell type column to pull from pi_file
n_boot    <- 101        # number of bootstrap iterations for the ensemble
seed      <- 123

## ---------------------------------------------------------------------------
## Run MiXcan2 ensemble for each gene
## ---------------------------------------------------------------------------

gene_names <- readLines(gene_list_file)

for (gene_name in gene_names) {
  
  print(paste("Processing:", gene_name))
  
  # Skip genes that have already been processed.
  out_marker <- paste0("ensbl_all_summary_", gene_name, ".RData")
  if (file.exists(out_marker)) {
    print(paste("Already done, skipping:", gene_name))
    next
  }
  
  gene_dir <- file.path(input_dir, gene_name)
  
  # Genotype matrix (N samples x P SNPs).
  X <- read.csv(file.path(gene_dir, "X.csv"), header = TRUE)
  names(X) <- NULL
  X <- as.matrix(X)
  
  # Gene expression vector (N samples).
  Y <- read.csv(file.path(gene_dir, "Y.csv"), header = TRUE)
  names(Y) <- NULL
  Y <- as.double(unlist(Y))
  
  # SNP identifiers for the columns of X.
  X_rows <- read.csv(file.path(gene_dir, "X_rows.csv"), header = TRUE)
  names(X_rows) <- NULL
  X_rows <- unlist(X_rows)
  
  # Covariates (N samples x Q covariates).
  Cov <- read.csv(file.path(gene_dir, "Cov.csv"))
  Cov <- mutate_all(Cov, as.numeric)
  Cov <- as_tibble(Cov)
  
  # Cell-type proportion estimates for the selected cell type.
  Pi_table <- read.delim(pi_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  if (!(cell_type %in% colnames(Pi_table))) {
    stop("Error: cell_type is not a valid column in Pi_table.")
  }
  Pi <- as.double(unlist(Pi_table %>% pull(!!cell_type)))
  
  # Fit the MiXcan2 ensemble model for this gene.
  ensbl <- MiXcan2_ensemble(
    y = Y, x = X, cov = Cov, pi = Pi,
    xName = X_rows, yName = gene_name,
    B = n_boot, seed = seed
  )
  
  # Save results.
  cell_dir <- file.path(output_dir, cell_type)
  if (!dir.exists(cell_dir)) dir.create(cell_dir, recursive = TRUE)
  
  write.csv(
    ensbl$ensemble_intecept,
    file = file.path(cell_dir, paste0("intercept_", gene_name, ".csv")),
    row.names = TRUE
  )
  saveRDS(ensbl$ensemble_weight,          file.path(cell_dir, paste0("ensemble_weights_", gene_name, ".rds")))
  saveRDS(ensbl$all_weights,              file.path(cell_dir, paste0("all_weights_df_", gene_name, ".rds")))
  saveRDS(ensbl$ensemble_summary,         file.path(cell_dir, paste0("ensbl_summary_", gene_name, ".rds")))
  saveRDS(ensbl$ensemble_summary_by_type, file.path(cell_dir, paste0("ensbl_summary_by_type_", gene_name, ".rds")))
  saveRDS(ensbl$all_summary,              file.path(cell_dir, paste0("ensbl_all_summary_", gene_name, ".rds")))
}
