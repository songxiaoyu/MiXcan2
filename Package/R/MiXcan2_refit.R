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
  x=model$x;y=model$y; cov=model$cov;
  type=model$type; pi=model$pi; foldid=model$foldid
  w <- MiXcan2_extract_weight(model = model, keepZeroWeight = T)
  summary <- MiXcan2_extract_summary(model = model)

  # NoPredictor - no performance
  if (type=="NoPredictor") {cv.r2.refit=in.sample.r2.refit=0}

  snpidx=Matrix::which(w$weight_cell_1 !=0 | w$weight_cell_2!=0)
  pr=length(snpidx)
  xreduced=as.matrix(x[, snpidx])
  # NonSpecific
  if (type=="NonSpecific") {

    xr_cov=cbind(xreduced, cov)
    ft=glmnet::glmnet(x=xr_cov, y=y, family = "gaussian",
                      alpha=0, lambda = 1e-3)
    beta=ft$beta[1:pr]
    # weight
    w$weight_cell_1[w$weight_cell_1!=0]=
      w$weight_cell_2[w$weight_cell_2!=0]=as.numeric(beta)

    # in sample c2
    y_hat=xreduced%*%beta
    in.sample.r2.refit=1-sum( (y-y_hat)^2)/sum(y^2)

    # cv r2
    all_r2=NULL
    for (i in 1:10) {
      temp=glmnet::glmnet(x=xr_cov[foldid!=i,], y=y[foldid!=i],
                          family="gaussian",
                          lambda = 1e-3, alpha=0)
      y_hat=as.matrix(xreduced[foldid==i,]) %*% temp$beta[1:pr]
      r2_temp=1-sum( (y[foldid==i]-y_hat)^2)/sum(y[foldid==i]^2)

      all_r2=c(all_r2, r2_temp)
    }
    all_r2[is.na(all_r2)]=0
    cv.r2.refit=mean(all_r2)
  }

  # CellTypeSpecific
  if (type=="CellTypeSpecific") {
    ci=pi-0.5; zreduced=ci*xreduced;
    xxreduced=as.matrix(cbind(ci, xreduced, zreduced))
    temp=glmnet::glmnet(x=xxreduced, y=y, family =
                        "gaussian", alpha=0, lambda = 0.001)
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

    # in sample r2
    y_hat=pi* cbind(1, xreduced) %*% tbeta1 +
      (1-pi)* cbind(1, xreduced) %*% tbeta2

    in.sample.r2.refit=1-sum( (y-y_hat)^2)/sum(y^2)


    # cv r2

    all_r2=NULL
    for (i in 1:10) {
      temp=glmnet::glmnet(x=xxreduced[foldid!=i,], y=y[foldid!=i],
                          family="gaussian",
                          lambda = 1e-3, alpha=0)
      test=c(temp$a0,as.numeric(temp$beta))
      tbeta10=test[1]+test[2]/2
      tbeta20=test[1]-test[2]/2
      tbeta11=test[3: (pr+2)] + test[(pr+3): (2*pr+2)]/2
      tbeta21=test[3: (pr+2)] - test[(pr+3): (2*pr+2)]/2
      tbeta1=c(tbeta10, tbeta11)
      tbeta2=c(tbeta20, tbeta21)

      tdesign=cbind(1, xreduced[foldid==i,] )
      y_hat= pi[foldid==i] * tdesign %*% tbeta1 +
        (1-pi[foldid==i]) * tdesign %*% tbeta2
      r2_temp=1-sum( (y[foldid==i]-y_hat)^2)/sum(y[foldid==i]^2)
      all_r2=c(all_r2, r2_temp)

    }
    all_r2[is.na(all_r2)]=0
    cv.r2.refit=mean(all_r2)
  }

  if (keepZeroWeight==F) {
    w = w %>%
      dplyr::filter(!(weight_cell_1 == 0 & weight_cell_2 == 0))
  }

  summary2= summary %>% data.frame()  %>%
    mutate(in_sample_r2_refit =in.sample.r2.refit)%>%
    mutate(cv_r2_refit =cv.r2.refit)


  return(list(weight=w, summary=summary2))
}

