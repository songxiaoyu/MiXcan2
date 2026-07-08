## =============================================================================
## Step 1 of association analysis: predict cell-type-specific gene expression
## from genotype dosages using MiXcan2 ensemble weights.
##
## NOTE: The genotype dosage data used here is NOT publicly available and is
## not included in this repository. Paths below are placeholders.
##
## Input:
##   - MiXcan2 ensemble weights CSV (output of run_mixcan2_ensemble.R), with
##     columns including yName (gene), xNameMatrix (SNP ID, "chrN_..."),
##     weight_cell_1, weight_cell_2.
##   - Per-chromosome genotype dosage files (samples x SNPs), one per
##     autosome, with a `variant_gtex_id` column identifying each SNP.
##
## Output:
##   - One prediction CSV per gene: predicted expression for cell_1 and
##     cell_2 in every sample, written to predictions/.
## =============================================================================

library(MiXcan2)
library(dplyr)
library(foreach)
library(doParallel)

## ---------------------------------------------------------------------------
## Config -- EDIT THESE
## ---------------------------------------------------------------------------

weights_file      <- "path/to/nonrefit_weights_wProportion_greater0.5.csv"  # from ensemble step
genotype_dir      <- "path/to/genotype_dosages"                            # per-chromosome dosage files
genotype_pattern  <- "dosages-chr(0[1-9]|1[0-9]|2[0-2])"                   # regex matching chr01-chr22 files
output_dir        <- "path/to/association_analysis"                        # predictions/ is created under here

n_cores <- max(1, detectCores() - 1)

## ---------------------------------------------------------------------------
## Helper functions
## ---------------------------------------------------------------------------

# Align a gene's SNP weights and genotype matrix to the same SNP set/order.
order_snps <- function(weight, new_x) {
  new_x  <- subset(new_x, rownames(new_x) %in% rownames(weight))
  weight <- subset(weight, rownames(weight) %in% rownames(new_x))
  
  common_row_names <- intersect(rownames(weight), rownames(new_x))
  weight_subset_ordered <- weight[match(common_row_names, rownames(weight)), ]
  new_x_ordered         <- new_x[match(common_row_names, rownames(new_x)), ]
  
  list(weight_subset_ordered = weight_subset_ordered, new_x_ordered = new_x_ordered)
}

# Predict cell-type-specific expression from genotypes and MiXcan2 weights.
MiXcan_prediction <- function(weight, new_x) {
  yhat_cell_1 <- t(new_x) %*% as.matrix(weight[, "weight_cell_1"])
  yhat_cell_2 <- t(new_x) %*% as.matrix(weight[, "weight_cell_2"])
  prediction  <- cbind(yhat_cell_1, yhat_cell_2)
  colnames(prediction) <- c("cell_1", "cell_2")
  prediction
}

## ---------------------------------------------------------------------------
## Run
## ---------------------------------------------------------------------------

setwd(output_dir)
dir.create("predictions", showWarnings = FALSE)

weight_result <- read.csv(weights_file)

# Pre-load one genotype dosage file per chromosome (keyed by zero-padded chromosome number).
registerDoParallel(cores = n_cores)

genotype_files <- list.files(path = genotype_dir, pattern = genotype_pattern, full.names = TRUE)
unique_chromosomes <- unique(sub("chr([0-9]+)_.*", "\\1", weight_result$xNameMatrix))

chromosome_files <- list()
for (chr in unique_chromosomes) {
  chr_padded <- ifelse(as.numeric(chr) < 10, sprintf("0%d", as.numeric(chr)), chr)
  genotype_file <- file.path(genotype_dir, paste0("dosages-chr", chr_padded))
  chromosome_files[[chr_padded]] <- read.table(genotype_file, header = TRUE, sep = "\t")
  rownames(chromosome_files[[chr_padded]]) <- chromosome_files[[chr_padded]]$variant_gtex_id
}

# Predict expression for each gene, in parallel.
foreach(gene = unique(weight_result$yName), .packages = c("foreach", "doParallel")) %dopar% {
  
  weight_test <- subset(weight_result, weight_result$yName == gene)
  rownames(weight_test) <- weight_test$xNameMatrix
  
  chromosome <- sub("chr([0-9]+)_.*", "\\1", weight_test$xNameMatrix[1])
  chromosome <- ifelse(as.numeric(chromosome) < 10, sprintf("0%d", as.numeric(chromosome)), chromosome)
  
  out_file <- paste0("predictions/prediction_", gene, "_chr", chromosome, ".csv")
  if (file.exists(out_file)) {
    cat("Prediction already exists for gene:", gene, "on chromosome", chromosome, "\n")
    return(NULL)
  }
  
  new_x <- chromosome_files[[chromosome]]
  
  aligned <- order_snps(weight_test, new_x)
  weight_subset_ordered <- aligned$weight_subset_ordered
  new_x_ordered <- aligned$new_x_ordered[, 4:ncol(aligned$new_x_ordered)]  # drop leading ID/annotation columns
  
  prediction <- MiXcan_prediction(weight_subset_ordered, new_x_ordered)
  rownames(prediction) <- gsub("^X", "", rownames(prediction))
  
  write.csv(prediction, file = out_file, row.names = TRUE)
  cat("Saved:", out_file, "\n")
}
