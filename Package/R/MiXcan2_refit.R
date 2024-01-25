#'  Un-penalized weights for MiXcan selected SNPs
#'
#' Refit MiXcan selected SNPs to  ordinary least square for un-penalized weights.
#'
#' @param model A direct output from MiXcan() function, which includes the
#' prediction weights as well as other summaries of the prediction models.
#' @param y: The pre-cleaned expression level data for a single gene in N samples.
#' @param x: A N by P matrix for all the genetic predictors used to predict the genetically regulated expression  of the gene.
#' @param cov: A N by Q matrix for the covariates adjusted in the model (e.g. age, population stratification).
#' @param pi: An estimation of cell-type faction of the cell type of interest (e.g.
#' epithelial). It can be estimated using existing methods
#' in the literature or from the output of pi_estimation function.
#' @param keepZeroWeight: Whether to keep predictors with zero weights.
#'
#' @return A data frame with weight for cell 1 and 2, including the potential meta data for the SNP/gene.
#' @export
#'
MiXcan2_refit <- function(model, x, y, cov, pi, foldid=NULL, keepZeroWeight=F) {
  MiXcan_weight_result <- MiXcan2_extract_weight(model = model, keepZeroWeight = F)
  MiXcan_weight_result2 <- MiXcan2_extract_weight(model = model, keepZeroWeight = T)
  summary <- MiXcan2_extract_summary(model = model, x=x, y=y, cov=cov, pi=pi)

  # clean predictors
  x=as.matrix(x); p=ncol(x); n=nrow(x)
  y=as.matrix(y)
  if (is.null(cov)==F) {
    cov=as.matrix(cov)
    res=lm(y~cov)$residuals
    res=as.matrix(res)
  } else {res=y}
  if (is.null(foldid)) {foldid= sample(1:10, n, replace=T)}


  # NoPredictor - no performance
  if (model$type=="NoPredictor") {
    cv.r2=in.r2=0
  }
  # NonSpecific
  if (model$type=="NonSpecific") {

    xreduced=as.matrix(x[, Matrix::which(as.numeric(model$glmnet.tissue$beta) !=0)])
    snpidx=Matrix::which(model$glmnet.tissue$beta[1:p]!=0)

    if (ncol(xreduced)>1) {
      ft=glmnet::glmnet(x=xreduced, y=res, family = "gaussian", alpha=0, lambda = 1e-3)
      ft0=glmnet::cv.glmnet(x=xreduced, y=res, family = "gaussian", alpha=0, lambda = c(1, 1e-3))
      beta=ft$beta #[1:length(snpidx)]
    } else { #  ncol(xreduced)==1
      xreduced2=cbind(1, xreduced)
      ft=glmnet::glmnet(x=xreduced2, y=res, family = "gaussian", alpha=0, lambda = 1e-3)
      ft0=glmnet::cv.glmnet(x=xreduced2, y=res, family = "gaussian", alpha=0, lambda = c(1, 1e-3))
      beta=ft$beta[-1]
    }

    cv.r2=max(0, 1-ft0$cvm[2]/var(res))
    in.r2=ft$dev.ratio

    MiXcan_weight_result$weight_cell_1=
      MiXcan_weight_result$weight_cell_2=as.numeric(beta)
    MiXcan_weight_result2$weight_cell_1[MiXcan_weight_result2$weight_cell_1!=0]=
      MiXcan_weight_result2$weight_cell_2[MiXcan_weight_result2$weight_cell_1!=0]=
      as.numeric(beta)
  }

  # CellTypeSpecific
  if (model$type=="CellTypeSpecific") {


    ci=pi-0.5; z=ci*x; z=as.matrix(z); xx=as.matrix(cbind(ci, x, z))
    idx=Matrix::which(model$glmnet.cell$beta!=0)
    xreduced=xx[,idx]
    ft=glmnet::glmnet(x=xreduced, y=res, family = "gaussian", alpha=0, lambda = 0.001)
    ft0=glmnet::cv.glmnet(x=xreduced, y=res, family = "gaussian", alpha=0, lambda = c(1, 1e-3))

    est0=as.numeric(ft$beta)
    beta=rep(0, length(model$glmnet.cell$beta))
    beta[idx]=est0
    est=c(0, beta)
    beta11=est[3: (p+2)] + est[(p+3): (2*p+2)]/2
    beta21=est[3: (p+2)] - est[(p+3): (2*p+2)]/2

    MiXcan_weight_result2$weight_cell_1=beta11
    MiXcan_weight_result2$weight_cell_2=beta21

    MiXcan_weight_result = MiXcan_weight_result2 %>%
      dplyr::filter(!(weight_cell_1 == 0 & weight_cell_1 == 0))

    cv.r2=max(0, 1-ft0$cvm[2]/var(res))
    in.r2=ft$dev.ratio

  }

  if (keepZeroWeight==F) {
    res=MiXcan_weight_result
  } else {res=MiXcan_weight_result2}

  summary2= summary %>% data.frame()  %>%
    mutate(cv_r2_refit =cv.r2)%>%
    mutate(in_sample_r2_refit=in.r2)


  return(list(weight=res, summary=summary2))
}

