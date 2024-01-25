#' Extract the gene expression prediction summary
#'
#' Extract the gene expression prediction summary (e.g. No of SNPs, R^2, correlation) function.
#' The output can be directly applied to GWAS data for cell-type specific TWAS.
#'
#' @param model A direct output from MiXcan() function, which includes the
#' prediction weights as well as other summaries of the prediction models.
#' @param y: The pre-cleaned expression level data for a single gene in N samples.
#' @param x: A N by P matrix for all the genetic predictors used to predict the genetically regulated expression  of the gene.
#' @param cov: A N by Q matrix for the covariates adjusted in the model (e.g. age, population stratification).
#' @param pi: An estimation of cell-type faction of the cell type of interest.
#'
#' @return A data frame with weight for cell 1 and 2, including the potential meta data for the SNP
#' @export
#'
MiXcan2_extract_summary <- function(model, x, y, cov=NULL, pi) {

  weights=model$beta.all.models

  # clean predictors
  x=as.matrix(x);
  y=as.matrix(y)
  if (is.null(cov)==F) {
    cov=as.matrix(cov)
    res=lm(y~cov)$residuals
    res=as.matrix(res)
  } else {res=y}

  # r calculation includes intercepts

  if (model$type=="NonSpecific") {
    cv_r2=model$ns.cv.r2
    in_sample_r2=model$ns.in.r2
  }
  if (model$type=="CellTypeSpecific") {
    cv_r2=model$cts.cv.r2
    in_sample_r2=model$cts.in.r2
  }
  if (model$type=="NoPredictor") {
    cv_r2=in_sample_r2=0
  }

  w2=MiXcan2_extract_weight(model=model, keepZeroWeight=F)


  summary=data.frame(yName=model$yName, n_snp_input=ncol(x),
                n_snp_model=nrow(w2),
                model_type=model$type,
                cv_r2=cv_r2,
                in_sample_r2=in_sample_r2)

  return(summary)
}

