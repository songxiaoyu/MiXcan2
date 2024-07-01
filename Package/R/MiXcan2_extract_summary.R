#' Extract the gene expression prediction summary
#'
#' Extract the gene expression prediction summaries, such as  No of SNPs, correlation, R^2.
#'
#' @param model A direct output from MiXcan2() function, which includes the
#' prediction weights as well as other summaries of the prediction models.
#'
#' @return  A row including summary statistics of the prediction model of the gene.
#' @export
#'
MiXcan2_extract_summary <- function(model) {

  w1=MiXcan2_extract_weight(model=model, keepZeroWeight=T)
  w2=MiXcan2_extract_weight(model=model, keepZeroWeight=F)

  summary=data.frame(yName=model$yName, n_snp_input=nrow(w1),
                n_snp_model=nrow(w2),
                model_type=model$type,
                in.sample.cor=model$in.sample.metrics[1],
                in.sample.R2=model$in.sample.metrics[2],
                cv.cor=model$cv.metrics[1],
                cv.R2=model$cv.metrics[2])

  return(summary)
}

