## =============================================================================
## Step 2 of association analysis: test predicted cell-type-specific gene
## expression (from predict_expression.R) for association with a phenotype
## (e.g. Dense Area, Non-Dense Area, Percent Density), stratified by
## measurement device/manufacturer.
##
## NOTE: The covariate, phenotype, and manufacturer-group data used here are
## NOT publicly available and are not included in this repository. Paths and
## group definitions below are placeholders.
##
## Input:
##   - Per-gene prediction CSVs from predict_expression.R
##     (predictions/prediction_<gene>_<chromosome>.csv)
##   - Covariates file (samples x covariates)
##   - Phenotype/outcome file (samples x 1)
##   - Manufacturer/group assignment file (sample_id, group)
##
## Output:
##   - One association result CSV per gene, per group
##   - One combined CSV of association results across all genes, per group
## The link to the genome-wide associaiton analysis results can be found here: https://zenodo.org/uploads/21268781
## =============================================================================

library(MiXcan2)
library(dplyr)
library(broom)
library(ACAT)

## ---------------------------------------------------------------------------
## Config -- EDIT THESE
## ---------------------------------------------------------------------------

weights_file       <- "path/to/weights.csv"  # from ensemble step, weights subset for cv.R2>0 and with a CTS or NS proportion > 0.5 (for MiXcan2 models)
predictions_dir    <- "path/to/association_analysis/predictions"
covariates_file     <- "path/to/covariates"
outcome_file        <- "path/to/outcome"     # e.g. dense area (DA), non-dense area (NDA), percent density (PD)
outcome_label       <- "DA"                  # used in output filenames

manufacturer_file  <- "path/to/manufacturer"  # sample_id, group (e.g. "GE" / "Hologic")
groups             <- c("Hologic", "GE")      # run the association analysis separately within each group

output_dir <- "path/to/association_analysis"

## ---------------------------------------------------------------------------
## Helper functions
## ---------------------------------------------------------------------------

# Align samples across prediction, covariate, and outcome tables.
order_individuals <- function(new_y, new_cov, new_outcome) {
  common_row_names <- intersect(intersect(rownames(new_y), rownames(new_cov)), rownames(new_outcome))
  list(
    new_y       = new_y[match(common_row_names, rownames(new_y)), ],
    new_cov     = new_cov[match(common_row_names, rownames(new_cov)), ],
    new_outcome = new_outcome[match(common_row_names, rownames(new_outcome)), ]
  )
}

# Fit a per-gene association model of the outcome on predicted cell_1/cell_2
# expression + covariates, combining the two cell-type p-values with ACAT.
MiXcan2_association <- function(new_y, new_cov, new_outcome, family = gaussian) {
  dat <- data.frame(new_y, new_cov, y = new_outcome)
  ft  <- glm(y ~ ., data = dat, family = family)
  res <- broom::tidy(ft)
  
  res1 <- res %>% filter(term == "cell_1")
  res2 <- res %>% filter(term == "cell_2")
  
  # If the two cell-type predictions are (anti-)collinear, cell_2's model is
  # degenerate -- reuse cell_1's result (flipping sign for anti-correlation).
  correlation <- tryCatch({
    cor_val <- cor(new_y[, 1], new_y[, 2], use = "complete.obs")
    if (is.na(cor_val)) stop("Correlation calculation resulted in NA")
    cor_val
  }, error = function(e) {
    cat("Error calculating correlation:", e$message, "\n")
    NA
  })
  
  if (!is.na(correlation)) {
    if (correlation > 0.99)  res2 <- res1
    if (correlation < -0.99) res2 <- res1 %>% mutate(estimate = -estimate)
  }
  
  p_values1 <- na.omit(res1$p.value)
  p_values2 <- na.omit(res2$p.value)
  if (length(p_values1) == 0) p_values1 <- 1
  if (length(p_values2) == 0) p_values2 <- 1
  
  p_combined <- tryCatch(
    ACAT::ACAT(c(p_values1, p_values2)),
    error = function(e) {
      cat("Error in ACAT:", e$message, "\n")
      NA
    }
  )
  
  res1 %>%
    dplyr::select(cell1_est = estimate, cell1_se = std.error, cell1_p = p.value) %>%
    bind_cols(res2 %>% dplyr::select(cell2_est = estimate, cell2_se = std.error, cell2_p = p.value)) %>%
    mutate(p_combined = p_combined) %>%
    as.data.frame()
}

## ---------------------------------------------------------------------------
## Run association analysis, separately per manufacturer/group
## ---------------------------------------------------------------------------

setwd(output_dir)

weight_result <- read.csv(weights_file)
manufacturer  <- read.table(manufacturer_file, header = FALSE)

for (group in groups) {
  
  dir.create(paste0("association_", group), showWarnings = FALSE)
  group_ids <- subset(manufacturer, V2 == group)$V1
  association_allgenes <- list()
  
  for (gene in unique(weight_result$yName)) {
    
    weight_test <- subset(weight_result, weight_result$yName == gene)
    rownames(weight_test) <- weight_test$xNameMatrix
    
    chromosome <- gsub("^chr([0-9]+)_.*", "\\1", weight_test$xNameMatrix[1])
    if (!is.na(as.numeric(chromosome))) chromosome <- sprintf("chr%02d", as.numeric(chromosome))
    
    filename <- file.path(predictions_dir, paste0("prediction_", gene, "_", chromosome, ".csv"))
    if (!file.exists(filename) || file.info(filename)$size == 0) {
      cat("Missing/empty prediction for gene:", gene, "chromosome:", chromosome, "\n")
      next
    }
    
    cov     <- read.table(covariates_file, header = FALSE, sep = "\t", row.names = 1)
    outcome <- read.table(outcome_file, header = FALSE, sep = "\t", row.names = 1)
    
    prediction <- read.csv(filename, row.names = 1)
    rownames(prediction) <- gsub("^X", "", rownames(prediction))
    
    # Restrict to samples in this manufacturer/group.
    cov        <- subset(cov, rownames(cov) %in% group_ids)
    outcome    <- subset(outcome, rownames(outcome) %in% group_ids)
    prediction <- subset(prediction, rownames(prediction) %in% group_ids)
    
    aligned <- order_individuals(prediction, cov, outcome)
    
    association_result <- MiXcan2_association(aligned$new_y, aligned$new_cov, aligned$new_outcome)
    association_result$Chromosome <- chromosome
    association_allgenes[[gene]] <- association_result
    
    out_file <- paste0("association_", group, "/association_", gene, "_", chromosome, ".csv")
    write.csv(association_result, file = out_file)
    cat("Saved:", out_file, "\n")
  }
  
  association_combined <- do.call(rbind, association_allgenes)
  combined_file <- paste0("association_allgenes_", outcome_label, "_", group, ".csv")
  write.csv(association_combined, file = combined_file, row.names = TRUE)
  cat("Saved combined results:", combined_file, "\n")
}
