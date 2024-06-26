#' Extract the gene expression prediction weights matrix
#'
#' Extract the gene expression prediction weights matrix (elastic net regression coefficients) from MiXcan() function.
#' The output can be directly applied to GWAS data for cell-type specific TWAS.
#' @import dplyr
#' @param model A direct output from MiXcan2() function, which includes the
#' prediction weights as well as other summaries of the prediction models.
#' @param keepZeroWeight Default=F. Whether to keep SNPs that have zero weights in all cell types.
#'
#' @return A data frame with weight for cell 1 and 2, including the potential meta data for the SNP/gene.
#' @export
#'
MiXcan2_extract_weight <- function(model, keepZeroWeight=F) {
  res <-
    model$beta.SNP.cell1 %>% dplyr::rename(weight_cell_1 = weight) %>%
    inner_join(model$beta.SNP.cell2 %>% dplyr::rename(weight_cell_2 = weight)) %>%
    mutate(type = model$type) %>%
    tibble::add_column(yName=model$yName, .before = 1)

  if (keepZeroWeight==F) {
    res = res %>% dplyr::filter(!(weight_cell_1 == 0 & weight_cell_2 == 0))
  }

  return(res)
}

