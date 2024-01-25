#' Estimate the cell-type level prediction weights of a gene.
#'
#' The core function of MiXcan package for estimating the
#' cell-type level prediction weights of a gene.
#' @param y The pre-cleaned expression level data for a single gene in N samples.
#' @param x A N by P matrix for all the genetic predictors used to predict the genetically regulated expression  of the gene.
#' @param cov A N by Q matrix for the covariates adjusted in the model (e.g. age, population stratification).
#' @param pi An estimation of cell-type fraction for the cell type of interest (e.g.
#' epithelial). It can be obtained using existing methods
#' in the literature or from the output of pi_estimation function.
#' @param xNameMatrix Default is NULL. A matrix to save the X matrix information,
#' such as variable ID, position, rsid, ref_allele, eff_allele.
#' @param yName Default is NULL. A row vector to save the expression information, such as gene ID, gene name.
#' @param foldid Default is NULL. 10-fold cross-validation (CV) is used in our pipeline. A random split
#' is considered if foldid is NULL. Otherwise foldid is used to split the data for CV.
#' @return list with 9 elements. It contains
#'
#' \item{type:}{Whether the prediction model is "CellTypeSpecific" or "NonSpecific".}
#' \item{beta.SNP.cell1:}{The prediction weights of the genetic predictors in cell type 1 (the cell type of interest).}
#' \item{beta.SNP.cell2:}{The prediction weights of the genetic predictors in cell type 2 (other cell types).}
#' \item{beta.all.models:}{All regression coefficients are saved in beta.all.models, including intercepts,
#' coefficients of genetic and non-genetic predictors in cell-type specific and non-specific models.}

#' \item{glmnet.cell:}{The cell-type-level prediction model selected using elastic-net.
#' This model may not be the final model as elastic-net selected parameters in the two
#' cell types may not be robustly different.  }
#' \item{glmnet.tissue:}{The tissue-level prediction model, which does not consider
#' cell type composition (same as PrediXcan).}
#' @export
#'
#'
#'
MiXcan2_model=function(y, x, cov=NULL, pi,
                 xNameMatrix=NULL, yName=NULL,
                foldid=NULL) {
  # clean predictors
  x=as.matrix(x);
  n=nrow(x); p=ncol(x)
  ci=pi-0.5; z=ci*x; z=as.matrix(z);
  xx=as.matrix(cbind(ci, x, z))

  # clean outcome -- get rid of the covariate effects
  y=as.matrix(y)
  if (is.null(cov)==F) {
    cov=as.matrix(cov)
    res=lm(y~cov)$residuals
    res=as.matrix(res)
  } else {res=y}

  # clean name
  if(is.null(yName)) {yName="Gene"}
  if(is.null(xNameMatrix)) {xNameMatrix=paste0("SNP", 1:p)}
  if (is.null(foldid)) {foldid= sample(1:10, n, replace=T)}


  # tissue model
  ft00=glmnet::cv.glmnet(x=x, y=res,family="gaussian",  foldid=foldid, alpha=0.5)
  ft0=glmnet::glmnet(x=x, y=res,  family="gaussian", lambda = ft00$lambda.1se, alpha=0.5)
  if (ft00$glmnet.fit$df[ft00$index[2]]==0) {
    est.tissue=c(ft0$a0, rep(0, p))
  } else {est.tissue=c(ft0$a0,as.numeric(ft0$beta))}

  # cell type specific model

  ft11=glmnet::cv.glmnet(x=xx, y=res,
                         penalty.factor=c(0, rep(1, ncol(xx)-1)),
                         family="gaussian", foldid=foldid, alpha=0.5)

  ft=glmnet::glmnet(x=xx, y=res, penalty.factor=c(0, rep(1, ncol(xx)-1)), family="gaussian",
                    lambda = ft11$lambda.1se, alpha=0.5)
  if (ft11$glmnet.fit$df[ft11$index[2]]==0 |
      (ft11$glmnet.fit$df[ft11$index[2]]==1 & ft11$glmnet.fit$beta[1]!=0) ) {
    est=c(ft$a0, ft$beta[1], rep(0, 2*p))
  } else {est=c(ft$a0, as.numeric(ft$beta))}

  beta10=est[1]+est[2]/2
  beta20=est[1]-est[2]/2
  beta11=est[3: (p+2)] + est[(p+3): (2*p+2)]/2
  beta21=est[3: (p+2)] - est[(p+3): (2*p+2)]/2


  ## add inference for difference > 0
  Type ="NonSpecific"
  idx.diff=((p+3): (2*p+2)) -1
  idx.nonzero=which(est[-1]!=0)
  idx.nonzero.diff=intersect(idx.diff, idx.nonzero)

  if (length(idx.nonzero.diff)!=0) {
    xx.select=xx[,idx.nonzero]
    beta.ols.boot=NULL; # beta.en.boot=NULL
    for(boot in 1:200) {
      id=sample(n, n, replace =T)
      gfit = lm(res[id,]~xx.select[id,])
      beta.ols.boot =rbind(beta.ols.boot, coef(gfit)[-1])
    }
    beta.range=apply(beta.ols.boot, 2, function(f) quantile(f, prob=c(0.025, 0.975), na.rm=T))
    beta.diff.range=beta.range[, match(idx.nonzero.diff, idx.nonzero)]
    # print(beta.diff.range)
    if (is.null(dim(beta.diff.range))) {
      any.nonzero= beta.diff.range[1] * beta.diff.range[2]>0} else {
        print(apply(beta.diff.range, 2, function(f) f[1] * f[2]>0))
        any.nonzero= any(apply(beta.diff.range, 2, function(f) f[1] * f[2]>0), na.rm=T)
      }

    if (is.na(any.nonzero)==F & any.nonzero==T) {Type ="CellTypeSpecific"}
  }

  # Given model type get weights


  if (Type!="CellTypeSpecific") {
    beta1=beta2=est.tissue
  } else {
      beta1=c(beta10, beta11)
      beta2=c(beta20, beta21)
  }
  beta.all.models=data.frame(est.tissue, beta1, beta2)
  colnames(beta.all.models)=c("Tissue", "Cell1", "Cell2")
  beta.SNP.cell1=data.frame(xNameMatrix, weight=beta1[2:(p+1)])
  beta.SNP.cell2=data.frame(xNameMatrix, weight=beta2[2:(p+1)])

  if (suppressWarnings(
    all(c(beta.SNP.cell1$weight, beta.SNP.cell2$weight)==0) )) {
    Type ="NoPredictor"}

  # Given model type get CV R2

  cts.cv.r2=1-ft11$cvm[ft11$index[2]]/var(res)
  cts.in.r2=ft$dev.ratio
  ns.cv.r2=1-ft00$cvm[ft00$index[2]]/var(res)
  ns.in.r2=ft0$dev.ratio



  return(list(type=Type,
              beta.SNP.cell1=beta.SNP.cell1,
              beta.SNP.cell2=beta.SNP.cell2,
              beta.all.models=beta.all.models,
              glmnet.cell=ft,
              glmnet.tissue=ft0,
              cts.cv.r2=cts.cv.r2,
              cts.in.r2=cts.in.r2,
              ns.cv.r2=ns.cv.r2,
              ns.in.r2=ns.in.r2,
              yName=yName,
              xNameMatrix=xNameMatrix))

}



































