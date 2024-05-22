#'  Un-penalized weights for MiXcan selected SNPs
#'
#' Refit MiXcan selected SNPs to  ordinary least square for un-penalized weights.
#'
#' @param model A direct output from MiXcan() function, which includes the
#' prediction weights as well as other summaries of the prediction models.
#' @param keepZeroWeight: Whether to keep predictors with zero weights.
#'
#' @return A data frame with weight for cell 1 and 2, including the potential meta data for the SNP/gene.
#' @export
#'
MiXcan2_refit <- function(model, keepZeroWeight=F) {

  x=model$x;y=model$y; cov=model$cov; n=nrow(x)
  type=model$type; pi=model$pi; foldid=model$foldid
  w <- MiXcan2_extract_weight(model = model, keepZeroWeight = T)
  summary <- MiXcan2_extract_summary(model = model)


  snpidx=Matrix::which(w$weight_cell_1 !=0 | w$weight_cell_2!=0)
  pr=length(snpidx)
  xreduced=as.matrix(x[, snpidx])
  ci=pi-0.5; zreduced=ci*xreduced;
  if (is.null(cov)==F) {
    xr_cov=cbind(xreduced, cov)
    xxreduced=as.matrix(cbind(ci, xreduced, zreduced, cov))
  } else {
    xr_cov=xreduced # nonspecific to use
    xxreduced=as.matrix(cbind(ci, xreduced, zreduced)) # CTS to use
  }

  # NoPredictor - no performance -----
  if (type=="NoPredictor") {in.sample.refit=cv.refit=rep(0, 4); intercept=rep(0, 2)}
  # NonSpecific ----------
  if (type=="NonSpecific") {
    ft=glmnet::glmnet(x=xr_cov, y=y, family = "gaussian", alpha=0, lambda = 1e-3)
    intercept=c(ft$a0, ft$a0)
    beta=ft$beta[1:pr]
    # weight
    w$weight_cell_1[w$weight_cell_1!=0]=
      w$weight_cell_2[w$weight_cell_2!=0]=as.numeric(beta)

    # in sample metrics
    y_hat=cbind(1, xreduced)%*%c(ft$a0, beta)
    if (is.null(cov)==F) {y_tilde=y-cov %*% ft$beta[(pr+1): (length(ft$beta))]} else {y_tilde=y}
    in.sample.refit=metrics(y_hat=y_hat, y_tilde=y_tilde, y=y)

    # cv metrics
    y_hat=y_tilde= rep(NA, n)
    for (i in 1:10) {
      temp=glmnet::glmnet(x=xr_cov[foldid!=i,], y=y[foldid!=i],family="gaussian",lambda = 1e-3, alpha=0)
      y_hat[foldid==i]=cbind(1, xreduced[foldid==i,]) %*% c(temp$a0, temp$beta[1:pr])
      if (is.null(cov)==F) {
        y_tilde[foldid==i]=y[foldid==i]-cov[foldid==i,] %*% temp$beta[(pr+1): (length(temp$beta))]
      } else {y_tilde[foldid==i]=y[foldid==i]}
    }
    cv.refit=metrics(y_hat=y_hat, y_tilde=y_tilde, y=y)
  }

  # CellTypeSpecific --------
  if (type=="CellTypeSpecific") {

    temp=glmnet::glmnet(x=xxreduced, y=y, family ="gaussian", alpha=0, lambda = 0.001)
    test=c(temp$a0,as.numeric(temp$beta))
    tbeta10=test[1]+test[2]/2
    tbeta20=test[1]-test[2]/2
    tbeta11=test[3: (pr+2)] + test[(pr+3): (2*pr+2)]/2
    tbeta21=test[3: (pr+2)] - test[(pr+3): (2*pr+2)]/2
    tbeta1=c(tbeta10, tbeta11)
    tbeta2=c(tbeta20, tbeta21)
    # weight
    w$weight_cell_1[w$weight_cell_1!=0]=tbeta11
    w$weight_cell_2[w$weight_cell_2!=0]=tbeta21
    intercept=c(tbeta10, tbeta20)

    # in sample metrics
    y_hat=pi* cbind(1, xreduced) %*% tbeta1 + (1-pi)* cbind(1, xreduced) %*% tbeta2
    if (is.null(cov)==F) {y_tilde=y-cov %*% test[(2*pr+3): (length(test))]} else {y_tilde=y}
    in.sample.refit=metrics(y_hat=y_hat, y_tilde=y_tilde, y=y)

    # cv metrics
    y_hat=y_tilde= rep(NA, n)
    for (i in 1:10) {
      temp=glmnet::glmnet(x=xxreduced[foldid!=i,], y=y[foldid!=i],family="gaussian",lambda = 1e-3, alpha=0)
      test=c(temp$a0,as.numeric(temp$beta))
      tbeta10=test[1]+test[2]/2
      tbeta20=test[1]-test[2]/2
      tbeta11=test[3: (pr+2)] + test[(pr+3): (2*pr+2)]/2
      tbeta21=test[3: (pr+2)] - test[(pr+3): (2*pr+2)]/2
      tbeta1=c(tbeta10, tbeta11)
      tbeta2=c(tbeta20, tbeta21)

      tdesign=cbind(1, xreduced[foldid==i,] )
      y_hat[foldid==i]= pi[foldid==i] * tdesign %*% tbeta1 +
        (1-pi[foldid==i]) * tdesign %*% tbeta2

      if (is.null(cov)==F) {
        y_tilde[foldid==i]=y[foldid==i]-as.matrix(cov[foldid==i,]) %*% test[(2*pr+3): (length(test))]
      } else (y_tilde[foldid==i]=y[foldid==i])

    }
    cv.refit=metrics(y_hat=y_hat, y_tilde=y_tilde, y=y)
  }
  # Clean output --------
  if (keepZeroWeight==F) {
    w = w %>%
      dplyr::filter(!(weight_cell_1 == 0 & weight_cell_2 == 0))
  }

  s2=data.frame(in.sample.unadj.cor.refit=in.sample.refit[1],
                in.sample.adj.cor.refit=in.sample.refit[2],
                 in.sample.unadj.R2.refit=in.sample.refit[3],
                 in.sample.adj.R2.refit=in.sample.refit[4],
                cv.unadj.cor.refit=cv.refit[1],
                cv.adj.cor.refit=cv.refit[2],
                cv.unadj.R2.refit=cv.refit[3],
                cv.adj.R2.refit=cv.refit[4])
  summary2= cbind(summary, s2)

  return(list(weight=w, summary=summary2, intercept=intercept))
}

