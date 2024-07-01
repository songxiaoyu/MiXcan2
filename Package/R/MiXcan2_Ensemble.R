#' MiXcan2 Ensemble
#'
#' MiXcan2 Ensemble
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
#' @return list with 9 elements. It contains
#'
#' \item{ensemble_summary}{Summarized model descriptions from all B ensemble models.}
#' \item{ensemble_summary_by_type}{Summarized model descriptions for Cell-Type-Specific and Non-Specific
#' models separately.}
#' \item{ensemble_weight}{Prediction weights from B ensemble models summarized separately for
#' Cell-Type-Specific and Non-Specific models.}
#' \item{ensemble_intecept}{Prediction intercepts from B ensemble models summarized separately for
#' Cell-Type-Specific and Non-Specific models.}
#' \item{CTS_weight}{Prediction weights summarized from all
#' Cell-Type-Specific models.}
#' \item{NS_weight}{Prediction weights summarized from all
#' Non-Specific models.}
#' \item{all_summary}{All B model descriptions w/o summarization (average).}
#' \item{all_weights}{All B weights w/o summarization (average).}
#' @export
#'
MiXcan2_ensemble=function(y, x, cov, pi, yName=NULL, xNameMatrix=NULL,
                          B=10, seed=NULL) {

  # start foreach
  `%dopar%` <- foreach::`%dopar%`
  res=foreach(i = 1:B, .combine='c') %dopar% {
    print(i)
    # set seed
    if (is.null(seed)==F) {seed1=seed+i*13; set.seed(seed1)}
    n=nrow(x)
    boot_idx=sample(n, round(n*0.9), replace=F)
    yboot=y[boot_idx]
    xboot=as.matrix(x[boot_idx,])
    piboot=pi[boot_idx]
    if (is.null(cov)==F) {covboot=as.matrix(cov)[boot_idx,]} else(covboot=NULL)

    foldid1 <- sample(1:10, length(yboot), replace=T)

    # run MiXcan
    model <- MiXcan2_model(y=yboot, x=xboot, cov = covboot,
                            pi= piboot, foldid = foldid1, yName=yName,
                           xNameMatrix=xNameMatrix)
    w <- MiXcan2_extract_weight(model = model, keepZeroWeight = T)
    summary <- MiXcan2_extract_summary(model=model)
    intercept=model$intercept

    return(list(list(w=w, intercept=intercept, summary=summary)))
  }


  # summarize B results - summary
  a=lapply(1:B, function(f) res[[f]]$summary) %>% rlist::list.rbind()
  sum1=a %>%  group_by(model_type) %>%
    summarise(across(n_snp_input:cv.R2, mean)) %>%
    mutate(Gene=a[1,1],.before = 1)
  sum2=a %>% mutate(CTS= 1*(model_type=="CellTypeSpecific")) %>%
    mutate(NS=1*(model_type=="NonSpecific")) %>%
    mutate(NP=1*(model_type=="NoPredictor")) %>%
    dplyr::select(!"model_type")  %>%
    summarise(across(n_snp_input:NP, mean)) %>%
    mutate(Gene=a[1,1],.before = 1)


  # summarize B results - weights
  wo=NULL; for (i in 1:B) {w=res[[i]]$w %>% mutate(ID=i, .before = 1);wo=rbind(wo, w)}
  name_temp=setdiff(colnames(wo), c("ID", "weight_cell_1", "weight_cell_2"))

  ensemble_weight= wo%>%
    group_by_at(name_temp) %>%
    summarise_at(vars("weight_cell_1", "weight_cell_2"), mean) %>%
    filter(weight_cell_1!=0 |weight_cell_2 !=0)

  CTS_weight=ensemble_weight %>% filter(type=="CellTypeSpecific")
  NS_weight=ensemble_weight %>% filter(type=="NonSpecific")

  # summarize B results - intercept
  all_intercept=lapply(1:B, function(f) res[[f]]$intercept)%>% rlist::list.rbind()
  mean_intercept=sapply(1:2, function(f) tapply(all_intercept[,f], a$model_type, mean)) %>%
    matrix(., ncol=2)
  colnames(mean_intercept)=c("intercept_cell_1", "intercept_cell_2")
  nt=names(table(a$model_type))
  rownames(mean_intercept)= nt

  return(list(ensemble_summary= sum2,
              ensemble_summary_by_type=sum1,
              ensemble_weight=ensemble_weight,
              ensemble_intecept=mean_intercept,
              CTS_weight=CTS_weight,
              NS_weight=NS_weight,
              all_summary=a,
              all_weights=wo))

}
