#' Find doublets/multiplets in UMI scRNA-seq data;
#'
#' Annotates doublets/multiplets using a binary classification approach to
#' discriminate artificial doublets from original data.
#'
#' @param sce single cell experiment (\code{SingleCellExperiment}) object to
#' analyze; needs \code{counts} in assays slot.
#' @param ntop integer, indicating number of top variance genes to consider.
#' Default: 500
#' @param srat numeric, indicating ratio between orginal number of "cells" and
#' simulated doublets; Default: 1
#' @param verb progress messages. Default: FALSE
#' @param retRes logical, should the trained classifier be returned?
#' Default: FALSE
#' @param nmax maximum number of training rounds; integer or "tune".
#'
#' Default: "tune"
#' @param varImp logical, should variable (i.e., gene) importance be returned?
#' Default: FALSE
#' @param estNdbl logical, should the numer of doublets be estimated from the data. Enables doublet calls. Default:FALSE. Use with caution.
#' @return sce input sce object \code{SingleCellExperiment} with doublet scores
#' added to colData as "bcds_score" column, and possibly more (details)
#' @importFrom Matrix Matrix rowSums rowMeans t
#' @importFrom stats var predict
#' @importFrom xgboost xgboost xgb.cv xgb.importance xgb.DMatrix
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assay assay<- assays assays<- assayNames
#' rowData rowData<- colData colData<-
#' @importFrom S4Vectors SimpleList
#' @importFrom S4Vectors metadata 'metadata<-'
#' @importFrom methods is
#' @export
#' @examples
#' data("sce_chcl")
#' ## create small data set using only 100 cells
#' sce_chcl_small = sce_chcl[, 1:100]
#' sce_chcl_small = bcds(sce_chcl_small)
bcds <- function(sce, ntop=500, srat=1, verb=FALSE, retRes=FALSE,
                 nmax="tune", varImp=FALSE, estNdbl = FALSE){

  #- check first argument (sce)
  if(!is(sce,"SingleCellExperiment"))
      stop('First argument (sce) needs to be of class SingleCellExperiment')
  if( !("counts" %in% names(assays(sce))) )
      stop('First argument (sce) needs to have a "counts" assay')
  if(!is(counts(sce),"sparseMatrix"))
      counts(sce) = Matrix::Matrix(counts(sce),sparse=TRUE)

  #- select variable genes
  if(verb) message("-> selecting genes\n")
  ind1            = Matrix::rowSums(counts(sce)>0)>0.01*ncol(sce)
  lc              = Matrix::t(log1p(counts(sce)[ind1,]))
  lc              = lc/Matrix::rowMeans(lc)
  vrs             = apply(lc,2,stats::var)
  hvg             = order(vrs,decreasing=TRUE)[seq_len(ntop)]
  lc              = lc[,hvg]
  lc              = lc/Matrix::rowMeans(lc)

  #- "simulated" doublets
  if(verb) message("-> simulating doublets\n")
  p1  = sample(seq_len(ncol(sce)),srat*ncol(sce),replace=TRUE)
  p2  = sample(seq_len(ncol(sce)),srat*ncol(sce),replace=TRUE)
  lc2 = Matrix::t(log1p(counts(sce)[ind1,p1][hvg,] + counts(sce)[ind1,p2][hvg,]))
  lc2 = lc2/Matrix::rowMeans(lc2)
  X   = rbind(lc,lc2)

  #- learn classifier
  if(verb) message("-> training classifier\n")
  colnames(X) = paste("GEN",seq_len(ncol(X)),sep="_") #- ranger picky
  y           = c(rep(-1,nrow(lc)),rep(1,nrow(lc2)))
  mm  = xgb.DMatrix(X,label=(y+1)/2)
  #- fixed rounds:
  if(nmax != "tune"){
    res    = NA
    varImp = FALSE
    retRes = FALSE
    pre    = xgboost(mm,nrounds=nmax,tree_method="hist",
                     nthread = 2, early_stopping_rounds = 2, subsample=0.5,
                     objective = "binary:logistic",verbose=0)
    sce$bcds_score = stats::predict(pre, newdat= mm[seq_len(ncol(sce)),])
    #- get doublet calls
    if(estNdbl){
      dbl.pre = stats::predict(pre, newdat= mm[seq(ncol(sce)+1,nrow(X)),])
      est_dbl = get_dblCalls_ALL(sce$bcds_score,dbl.pre,rel_loss=srat)
      if(is.null(metadata(sce)$cxds)) metadata(sce)$cxds = list()
      metadata(sce)$bcds$ndbl = est_dbl
      metadata(sce)$bcds$sim_scores = dbl.pre
      sce$bcds_call = sce$bcds_score >= est_dbl["balanced","threshold"]
    }
  #- learning rounds with CV:
  } else {
    res = xgb.cv(data =mm, nthread = 2, nrounds = 500, objective = "binary:logistic",
                 nfold=5,metrics=list("error"),prediction=TRUE,
                 early_stopping_rounds=2, tree_method="hist",subsample=0.5,verbose=0)
    ni  = res$best_iteration
    ac  = res$evaluation_log$test_error_mean[ni] + 1*res$evaluation_log$test_error_std[ni]
    ni  = min(which( res$evaluation_log$test_error_mean <= ac  ))
    nmax = ni
    pre = xgboost(mm,nrounds=nmax,tree_method="hist",
                  nthread = 2, early_stopping_rounds = 2, subsample=0.5,
                  objective = "binary:logistic",verbose=0)
    sce$bcds_score = res$pred[seq_len(ncol(sce))]
    #- get doublet calls
    if(estNdbl){
      dbl.pre = stats::predict(pre, newdat= mm[seq(ncol(sce)+1,nrow(X)),])
      est_dbl = get_dblCalls_ALL(sce$bcds_score,dbl.pre,rel_loss=srat)
      if(is.null(metadata(sce)$cxds)) metadata(sce)$cxds = list()
      metadata(sce)$bcds$ndbl = est_dbl
      metadata(sce)$bcds$sim_scores = dbl.pre
      sce$bcds_call = sce$bcds_score >= est_dbl["balanced","threshold"]

    }
  }
  #- variable importance
  if(varImp){
    if(verb) message("-> calculating variable importance\n")
    vimp = xgb.importance(model=pre)
    vimp$col_index = match(vimp$Feature,colnames(X))
  }
  #- result
  if(retRes){
    hvg_bool                    = (seq_len(nrow(sce))) %in% which(ind1)
    hvg_bool[hvg_bool]          = (seq_len(sum(hvg_bool))) %in% hvg
    hvg_ord                     = rep(NA,nrow(sce))
    hvg_ord[hvg_bool]           = which(ind1)[hvg]
    rowData(sce)$bcds_hvg_bool  = hvg_bool
    rowData(sce)$bcds_hvg_ordr  = hvg_ord
    metadata(sce)$bcds_res_cv   = res
    metadata(sce)$bcds_res_all  = pre
    metadata(sce)$bcds_nmax     = nmax
    if(varImp){
      vimp$gene_index             = hvg_ord[hvg_bool][vimp$col_index]
      metadata(sce)$bcds_vimp     = vimp[seq_len(100),-c(1,5)]
    }
  }

  if(verb) message("-> done.\n\n")
  return(sce)

}
