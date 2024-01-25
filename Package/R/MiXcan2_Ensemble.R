#' MiXcan2 Ensemble
#'
#' MiXcan2 Ensemble
#'
#' @param y The pre-cleaned expression level data for a single gene in N samples.
#' @param x A N by P matrix for all the genetic predictors used to predict the genetically regulated expression  of the gene.
#' @param cov A N by Q matrix for the covariates adjusted in the model (e.g. age, population stratification).
#' @param pi An estimation of cell-type faction of the cell type of interest.
#' @param xNameMatrix Default is NULL. A matrix to save the X matrix information,
#' such as variable ID, position, rsid, ref_allele, eff_allele.
#' @param yName Default is NULL. A row vector to save the expression information, such as gene ID, gene name.
#' @param B No. of ensemble models
#' @param seed seed
#'
#' @return A data frame with weight for cell 1 and 2, including the potential meta data for the SNP
#' @export
#'
#'
MiXcan2_ensemble=function(y, x, cov, pi, yName=NULL, xNameMatrix=NULL,
                          B=10, seed=NULL) {

  # start foreach
  res=foreach(i = 1:B, .combine='c') %dopar% {
    print(i)
    # set seed
    if (is.null(seed)==F) {seed1=seed+i*13; set.seed(seed1)}

    foldid1 <- sample(1:10, length(y), replace=T)


    # run MiXcan
    model <- MiXcan2_model(y=y, x=x, cov = cov,
                            pi= pi,
                            foldid = foldid1, yName=yName,
                           xNameMatrix=xNameMatrix)

    # Refit
    MiXcan_refit <- MiXcan2_refit(model = model,
                                               y=y,
                                               x=x, cov = cov,
                                               pi= pi, keepZeroWeight=T)

    return(list(MiXcan_refit))
  }

  # summarize B results -- summary
  a=lapply(1:B, function(f) res[[f]]$summary)
  b=list.rbind(a)%>%
    mutate(NP= 1*(model_type=="NoPredictor"))%>%
    mutate(CTS= 1*(model_type=="CellTypeSpecific"))%>%
    mutate(NS= 1*(model_type=="NonSpecific"))%>%
    dplyr::select(!"model_type") %>%
    mutate(no_of_models=B, .before="cv_r2")

  sum=apply(b[,-1], 2, mean) %>% as.data.frame() %>% t()%>% as.data.frame() %>%
    mutate(Gene=b[1,1],.before = 1)

  # summarize B results -- weights

  ww=NULL
  for (i in 1:B) {
    w=res[[i]]$weight %>%
      mutate(ID=i, .before = 1)
    ww=rbind(ww, w)
  }

  name_temp=setdiff(colnames(ww), c("ID", "weight_cell_1", "weight_cell_2"))

  ensemble_weight= ww%>%
     group_by_at(name_temp) %>%
    summarise_at(vars("weight_cell_1", "weight_cell_2"), mean) %>%
    filter(weight_cell_1!=0 |weight_cell_2 !=0)


  return(list(summary=sum, ensemble_weight=ensemble_weight, all_weights=ww))

}
