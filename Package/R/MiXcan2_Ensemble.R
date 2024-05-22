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
#' @return A data frame with weight for cell 1 and 2, including the potential meta data for the SNP
#' @export
#'
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
                            pi= piboot,
                            foldid = foldid1, yName=yName,
                           xNameMatrix=xNameMatrix)
    w <- MiXcan2_extract_weight(model = model, keepZeroWeight = T)
    intercept=model$intercept

    # Refit
    MiXcan_refit <- MiXcan2_refit(model = model, keepZeroWeight=T)

    return(list(list(original_w=w, original_intercept=intercept, refit=MiXcan_refit)))
  }

  # summarize B results -- summary
  a=lapply(1:B, function(f) res[[f]]$refit$summary) %>% rlist::list.rbind()
  sum1=a %>%  group_by(model_type) %>%
     summarise(across(n_snp_input:cv.adj.R2.refit, mean)) %>%
     mutate(Gene=a[1,1],.before = 1)
  sum2=a %>% mutate(CTS= 1*(model_type=="CellTypeSpecific")) %>%
    mutate(NS=1*(model_type=="NonSpecific")) %>%
    mutate(NP=1*(model_type=="NoPredictor")) %>%
    dplyr::select(!"model_type")  %>%
    summarise(across(n_snp_input:NP, mean)) %>%
    mutate(Gene=a[1,1],.before = 1)
  # summarize B results -- original weights
  wo=NULL; for (i in 1:B) {w=res[[i]]$original_w %>% mutate(ID=i, .before = 1);wo=rbind(wo, w)}
  name_temp=setdiff(colnames(wo), c("ID", "weight_cell_1", "weight_cell_2"))

  original_ensemble_weight= wo%>%
    group_by_at(name_temp) %>%
    summarise_at(vars("weight_cell_1", "weight_cell_2"), mean) %>%
    filter(weight_cell_1!=0 |weight_cell_2 !=0)

  original_CTS_weight=original_ensemble_weight %>%
    filter(type=="CellTypeSpecific")
  original_NS_weight=original_ensemble_weight %>%
    filter(type=="NonSpecific")

  # summarize B results -- original intercept
  o_all_intercept=lapply(1:B, function(f) res[[f]]$original_intercept)%>% rlist::list.rbind()
  original_mean_intercept=sapply(1:2, function(f) tapply(o_all_intercept[,f], a$model_type, mean)) %>%
    matrix(., ncol=2)
  colnames(original_mean_intercept)=c("intercept_cell_1", "intercept_cell_2")
  if (nrow(original_mean_intercept)==2) { rownames(original_mean_intercept)= C("CTS", "NS")}

  # summarize B results -- refit weights
  ww=NULL;
  for (i in 1:B) {w=res[[i]]$refit$weight %>% mutate(ID=i, .before = 1);ww=rbind(ww, w)}
  name_temp=setdiff(colnames(ww), c("ID", "weight_cell_1", "weight_cell_2"))

  refit_ensemble_weight= ww%>%
     group_by_at(name_temp) %>%
    summarise_at(vars("weight_cell_1", "weight_cell_2"), mean) %>%
    filter(weight_cell_1!=0 |weight_cell_2 !=0)

  refit_CTS_weight=refit_ensemble_weight %>%
    filter(type=="CellTypeSpecific")
  refit_NS_weight=refit_ensemble_weight %>%
    filter(type=="NonSpecific")

  # summarize B results -- refit intercept
  refit_all_intercept=lapply(1:B, function(f) res[[f]]$refit$intercept)%>% rlist::list.rbind()
  refit_mean_intercept=sapply(1:2, function(f) tapply(refit_all_intercept[,f], a$model_type, mean))%>%
    matrix(., ncol=2)
  colnames(refit_mean_intercept)=c("intercept_cell_1", "intercept_cell_2")
  if (nrow(refit_mean_intercept)==2) { rownames(refit_mean_intercept)= C("CTS", "NS")}

  return(list(ensemble_summary= sum2,
              ensemble_summary_by_type=sum1,
              refit_ensemble_weight=refit_ensemble_weight,
              refit_CTS_weight=refit_CTS_weight,
              refit_NS_weight=refit_NS_weight,
              refit_intecept=refit_mean_intercept,
              original_ensemble_weight=original_ensemble_weight,
              original_CTS_weight=original_CTS_weight,
              original_NS_weight=original_NS_weight,
              original_intecept=original_mean_intercept,
              all_summary=a,
              refit_all_weights=ww))

}
