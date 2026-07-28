#' Predict the cell-type specific or non-specific expression levels of a gene with MiXcan model in new genetic data.
#'
#' @param weight Weight matrix produced by MiXcan_extract_weight()
#' @param new_x The genetic data used for prediction in a new data set.
#'
#' @return A N by 2 matrix indicating the predicted gene expression levels in two
#' cell types. If a non-specific model is used for prediction, the predicted values should be the same in two cell types.
#' @export
#'

# wrapper function for MiXcan_predicion (from MiXcan function) that subsets the weight matrix by SNPs (rownames) in new_x and orders
# the dataframes 

MiXcan_prediction_wrapped <- function(weight, new_x) {
  # Extract the SNP IDs from new_x
  snp_ids <- rownames(new_x)
  
  # Subset weights by SNP IDs
  weight_subset <- weight[weight$SNP %in% snp_ids, ]
  
  # Order both weight_subset and new_x by SNP IDs
  weight_subset <- weight_subset[match(snp_ids, weight_subset$SNP), ]
  new_x_ordered <- new_x[match(rownames(weight_subset), rownames(new_x)), ]
  
  # Perform prediction
  yhat_MiXcan_cell_1 <- new_x_ordered %*% as.matrix(weight_subset[,"weight_cell_1"])
  yhat_MiXcan_cell_2 <- new_x_ordered %*% as.matrix(weight_subset[,"weight_cell_2"])
  yhat_MiXcan_prediction <- cbind(yhat_MiXcan_cell_1, yhat_MiXcan_cell_2)
  colnames(yhat_MiXcan_prediction) <- c("cell_1", "cell_2")
  
  return(yhat_MiXcan_prediction)
}
