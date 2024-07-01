#' Estimate the cell-type level prediction weights of a gene.
#'
#' The core function of MiXcan package for estimating the
#' cell-type level prediction weights of a gene.
#' @param y The pre-cleaned expression level data for a single gene in N samples.
#' @param x A N by P matrix for all the genetic predictors used to predict the genetically regulated expression  of the gene.
#' @param cov A N by Q matrix for the covariates adjusted in the model (e.g. age, population stratification).
#' @param pi An estimation of cell-type fraction for the cell type of interest (e.g.
#' epithelial). It can be obtained using existing methods
#' in the literature or from the output of `pi_estimation` function.
#' @param xNameMatrix Default is NULL. A matrix to save the X matrix information,
#' such as variable ID, position, rsid, ref_allele, eff_allele.
#' @param yName Default is NULL. A row vector to save the expression information, such as gene ID, gene name.
#' @param foldid Default is NULL. 10-fold cross-validation (CV) is used in our pipeline. A random split
#' is considered if foldid is NULL. Otherwise foldid is used to split the data for CV.
#' @return list with 9 elements. It contains
#'
#' \item{type}{Whether the prediction model is "CellTypeSpecific" or "NonSpecific" or "NoPredictor".}
#' \item{beta.SNP.cell1}{The prediction weights of the genetic predictors in cell type 1 (the cell type of interest).}
#' \item{beta.SNP.cell2}{The prediction weights of the genetic predictors in cell type 2 (other cell types).}
#' \item{beta.all.models}{All regression coefficients are saved in beta.all.models, including intercepts,
#' coefficients of genetic and non-genetic predictors in cell-type specific and non-specific models.}
#' \item{glmnet.cell}{The cell-type-level prediction model selected using elastic-net.
#' This model may not be the final model as elastic-net selected parameters in the two
#' cell types may not be robustly different.  }
#' \item{glmnet.tissue}{The tissue-level prediction model, which does not consider
#' cell type composition (same as PrediXcan).}
#' \item{in.sample.metrics}{In-sample correlation and R2.}
#' \item{cv.metrics}{CV correlation and R2.}
#' \item{intercept}{Intercepts from the two cell types.}
#' \item{xNameMatrix}{Output xNameMatrix for summary information.}
#' \item{yName}{Output yName for summary information.}
#' @export
#'
#'
#'


MiXcan2_model=function(y, x, cov=NULL, pi,
                 xNameMatrix=NULL, yName=NULL,
                foldid=NULL) {
  # format input
  x=as.matrix(x); y=as.matrix(y);n=nrow(x); p=ncol(x)
  # clean name
  if(is.null(yName)) {yName="Gene"}
  if(is.null(xNameMatrix)) {xNameMatrix=paste0("SNP", 1:p)}
  if (is.null(foldid)) {foldid= sample(1:10, n, replace=T)}

  if (is.null(cov)) {
    pcov=0; xcov=x;
    ci=pi-0.5; z=ci*x; xx=as.matrix(cbind(ci, x, z))
  } else {
    cov=as.matrix(cov)
    pcov=ncol(cov); xcov=as.matrix(cbind(x, cov))
    ci=pi-0.5; z=ci*x;xx=as.matrix(cbind(ci, x, z, cov))
  }

  # tissue model
  ft00=glmnet::cv.glmnet(x=xcov, y=y,family="gaussian",  foldid=foldid, alpha=0.5)
  ft0=glmnet::glmnet(x=xcov, y=y,  family="gaussian", lambda = ft00$lambda.1se, alpha=0.5)
  est.tissue=c(ft0$a0,as.numeric(ft0$beta))

  # cell type specific model
  ft11=glmnet::cv.glmnet(x=xx, y=y,penalty.factor=c(0, rep(1, ncol(xx)-1)),
                         family="gaussian", foldid=foldid, alpha=0.5)
  ft=glmnet::glmnet(x=xx, y=y, penalty.factor=c(0, rep(1, ncol(xx)-1)),
                    family="gaussian", lambda = ft11$lambda.1se, alpha=0.5)
  est=c(ft$a0,as.numeric(ft$beta))
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
      gfit = lm(y[id,]~xx.select[id,])
      beta.ols.boot =rbind(beta.ols.boot, coef(gfit)[-1])
    }
    beta.range=apply(beta.ols.boot, 2, function(f)
      quantile(f, prob=c(0.025, 0.975), na.rm=T))
    beta.diff.range=beta.range[, match(idx.nonzero.diff, idx.nonzero)]
    # print(beta.diff.range)
    if (is.null(dim(beta.diff.range))) {
      any.nonzero= beta.diff.range[1] * beta.diff.range[2]>0} else {
        print(apply(beta.diff.range, 2, function(f) f[1] * f[2]>0))
        any.nonzero= any(apply(beta.diff.range, 2, function(f) f[1] * f[2]>0), na.rm=T)
      }

    if (is.na(any.nonzero)==F & any.nonzero==T) {Type ="CellTypeSpecific"}
  }


  if (Type!="CellTypeSpecific") { # NonSpecific
    beta1=beta2=est.tissue
  } else { # CellTypeSpecific
    if (is.null(cov)) {
      beta1=c(beta10, beta11)
      beta2=c(beta20, beta21)
    } else {
      beta1=c(beta10, beta11, est[ (2*p+3): (2*p+2+pcov)])
      beta2=c(beta20, beta21, est[ (2*p+3): (2*p+2+pcov)])
    }
  }
  beta.all.models=cbind(est.tissue, beta1, beta2)
  colnames(beta.all.models)=c("Tissue", "Cell1", "Cell2")
  beta.SNP.cell1=data.frame(xNameMatrix, weight=beta1[2:(p+1)])
  beta.SNP.cell2=data.frame(xNameMatrix, weight=beta2[2:(p+1)])
  intercept=beta.all.models[1, 2:3]

  if (suppressWarnings(
    all(c(beta.SNP.cell1$weight, beta.SNP.cell2$weight)==0) )) {
    Type ="NoPredictor"}

  # ---- get in sample metrics ------
  if (Type=="NonSpecific" | Type=="CellTypeSpecific") {
    design=cbind(1,x)

    y_hat=pi*(design %*% beta.all.models[1:(1+p), "Cell1"])+
      (1-pi)*(design %*% beta.all.models[1:(1+p), "Cell2"])

    if (is.null(cov)==F) {
      beta_cov = beta1[ (p+2): length(beta1)]
      y_tilde=y-cov %*%beta_cov
    } else {y_tilde=y}

    in.sample=metrics(y_hat=y_hat, y_tilde=y_tilde)

  } else {in.sample=rep(0, 2)}


  # ---- get CV metrics ------
  if (Type=="NonSpecific") {
    y_hat=y_tilde=rep(NA, n)
    for (i in 1:10) {
      temp=glmnet::glmnet(x=as.matrix(xcov[foldid!=i,]), y=y[foldid!=i],
                          family="gaussian", lambda = ft00$lambda.1se, alpha=0.5)

      y_hat[foldid==i]=cbind(1, x[foldid==i,]) %*% c(temp$a0, temp$beta[1:p])

      if (is.null(cov)==F) {
        beta_cov = temp$beta[(p+1): length(temp$beta)]
        y_tilde[foldid==i]=y[foldid==i]-cov[foldid==i,] %*%beta_cov
      } else {y_tilde[foldid==i]=y[foldid==i]}
    }
    cv=metrics(y_hat=y_hat, y_tilde=y_tilde)

  }

  if (Type=="CellTypeSpecific") {
    y_hat=y_tilde=rep(NA, n)
    for (i in 1:10) {
      temp=glmnet::glmnet(x=xx[foldid!=i, ], y=y[foldid!=i],
                          penalty.factor=c(0, rep(1, ncol(xx)-1)),
                          family="gaussian", lambda = ft11$lambda.1se, alpha=0.5)

      test=c(temp$a0,as.numeric(temp$beta))
      tbeta10=test[1]+test[2]/2
      tbeta20=test[1]-test[2]/2
      tbeta11=test[3: (p+2)] + test[(p+3): (2*p+2)]/2
      tbeta21=test[3: (p+2)] - test[(p+3): (2*p+2)]/2
      tbeta1=c(tbeta10, tbeta11)
      tbeta2=c(tbeta20, tbeta21)

      tdesign=cbind(1, x[foldid==i, ] )

      y_hat[foldid==i]= pi[foldid==i] * tdesign %*% tbeta1 +
        (1-pi[foldid==i]) * tdesign %*% tbeta2

      if (is.null(cov)==F) {
        beta_cov = test[ (2*p+3): (2*p+2+pcov)]
        y_tilde[foldid==i]=y[foldid==i]-cov[foldid==i,] %*%beta_cov
      } else {y_tilde[foldid==i]=y[foldid==i]}

    }
    cv=metrics(y_hat=y_hat, y_tilde=y_tilde)
  }

  if (Type =="NoPredictor") {cv=rep(0, 2)}

  return(list(type=Type,
              beta.SNP.cell1=beta.SNP.cell1,
              beta.SNP.cell2=beta.SNP.cell2,
              beta.all.models=beta.all.models,
              glmnet.cell=ft,
              glmnet.tissue=ft0,
              in.sample.metrics=in.sample,
              cv.metrics=cv,
              intercept=intercept,
              xNameMatrix=xNameMatrix,
              yName=yName
              ))

}


metrics=function(y_hat,  # should include intercept
                 y_tilde # y-z*gamma
                 ) {
  adj.cor= cor(y_hat, y_tilde)
  adj.R2= 1-sum((y_tilde-y_hat)^2)/sum((y_tilde-mean(y_tilde))^2)
  return(c(adj.cor, adj.R2))
}
































