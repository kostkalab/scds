#' Find doublets/multiples in UMI scRNA-seq data;
#'
#' @param sce single cell experiment to analyze; needs "counts" in assays slot.
#' @param ntop number of top variance genes to consider. Default: 500
#' @param srat ratio between orginal "cells" and simulated doubltes; Default:1
#' @param verb progress messages. Default: FALSE
#' @param retRes shuld the trained classifier be returned? Default: FALSE
#' @param nmax maximum number of training rounds; interger or "tune". Default: "tune"
#' @param varImp shold variable (i.e., gene) importance be returned? Default: FALSE
#' @return sce Input sce with doublet scores added to colData as "bcds_score" column, and possibly more (details)
#' @importFrom Matrix Matrix rowSums rowMeans
#' @importFrom stats var predict
#' @importFrom xgboost xgboost xgb.cv xgb.importance xgb.DMatrix
#' @import  SingleCellExperiment
#' @importFrom SummarizedExperiment assay assay<- assays assays<- assayNames rowData rowData<- colData colData<-
#' @importFrom S4Vectors SimpleList
#' @importFrom S4Vectors metadata 'metadata<-'
#' @export
bcds <- function(sce,ntop=500,srat=1,verb=FALSE, retRes=FALSE, nmax="tune", varImp=FALSE){
#=========================================================================================

  if(!is(counts(sce),"sparseMatrix")) counts(sce) = Matrix::Matrix(counts(sce),sparse=TRUE)

  #- select variable genes
  if(verb) cat("-> selecting genes\n")
  ind1            = Matrix::rowSums(counts(sce)>0)>0.01*ncol(sce)
  lc              = t(log1p(counts(sce)[ind1,]))
  lc              = lc/Matrix::rowMeans(lc)
  vrs             = apply(lc,2,stats::var)
  hvg             = order(vrs,decreasing=TRUE)[1:ntop]
  lc              = lc[,hvg]
  lc              = lc/Matrix::rowMeans(lc)

  #- "simulated" doublets
  if(verb) cat("-> simulating doublets\n")
  p1  = sample(1:ncol(sce),srat*ncol(sce),replace=TRUE)
  p2  = sample(1:ncol(sce),srat*ncol(sce),replace=TRUE)
  lc2 = t(log1p(counts(sce)[ind1,p1][hvg,] + counts(sce)[ind1,p2][hvg,]))
  lc2 = lc2/Matrix::rowMeans(lc2)
  X   = rbind(lc,lc2)

  #- learn classifier
  if(verb) cat("-> training classifier\n")
  colnames(X) = paste("GEN",1:ncol(X),sep="_") #- ranger picky
  y           = c(rep(-1,nrow(lc)),rep(1,nrow(lc2)))
  mm  = xgb.DMatrix(X,label=(y+1)/2)
  #- fixed rounds:
  if(nmax != "tune"){
    res    = NA
    varImp = FALSE
    retRes = FALSE
    pre    = xgboost(mm,nrounds=nmax,tree_method="hist",
                     nthread = 2, early_stopping_rounds = 2, subsample=0.5,objective = "binary:logistic")
    sce$bcds_score = stats::predict(pre, newdat= mm[1:ncol(sce),])
    #- learning rounds with CV:
  } else {
    res = xgb.cv(data =mm, nthread = 2, nrounds = 500, objective = "binary:logistic",
                 nfold=5,metrics=list("error"),prediction=TRUE,early_stopping_rounds=2, tree_method="hist",subsample=0.5)
    ni  = res$best_iteration
    ac  = res$evaluation_log$test_error_mean[ni] + 1*res$evaluation_log$test_error_std[ni]
    ni  = min(which( res$evaluation_log$test_error_mean <= ac  ))
    nmax = ni
    pre = xgboost(mm,nrounds=nmax,tree_method="hist",
                  nthread = 2, early_stopping_rounds = 2, subsample=0.5,objective = "binary:logistic")
    sce$bcds_score = res$pred[1:ncol(sce)]
  }
  #- variable importance
  if(varImp){
    if(verb) cat("-> calculaing variable importance\n")
    vimp = xgb.importance(model=pre)
    vimp$col_index = match(vimp$Feature,colnames(X))
  }
  #- result
  if(retRes){
    hvg_bool                    = (1:nrow(sce)) %in% which(ind1)
    hvg_bool[hvg_bool]          = (1:sum(hvg_bool)) %in% hvg
    hvg_ord                     = rep(NA,nrow(sce))
    hvg_ord[hvg_bool]           = which(ind1)[hvg]
    rowData(sce)$bcds_hvg_bool  = hvg_bool
    rowData(sce)$bcds_hvg_ordr  = hvg_ord
    metadata(sce)$bcds_res_cv   = res
    metadata(sce)$bcds_res_all  = pre
    metadata(sce)$bcds_nmax     = nmax
    vimp$gene_index              = hvg_ord[hvg_bool][vimp$col_index]
    metadata(sce)$bcds_vimp     = vimp[1:100,-c(1,5)]
  }

  if(verb) cat("-> done.\n\n")
  return(sce)

}
