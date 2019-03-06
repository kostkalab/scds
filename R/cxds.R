#' Find doublets/multiplets in UMI scRNA-seq data;
#'
#' Annotates doublets/multiplets using co-expression based approach
#' 
#' @param sce single cell experiment (\code{SingleCellExperiment}) object to analyze; needs \code{counts} in assays slot.
#' @param ntop integer, indicating number of top variance genes to consider. Default: 500
#' @param binThresh integer, minimum counts to consider a gene "present" in a cell. Default: 0
#' @param verb progress messages. Default: FALSE
#' @param retRes logical, whether to return gene pair scores & top-scoring gene pairs? Default: FALSE.
#' @return sce input sce object \code{SingleCellExperiment} with doublet scores added to colData as "cxds_score" column.
#' @importFrom Matrix Matrix rowSums rowMeans t
#' @import SingleCellExperiment
#' @importFrom SummarizedExperiment assay assay<- assays assays<- assayNames rowData rowData<- colData colData<-
#' @importFrom S4Vectors SimpleList
#' @importFrom S4Vectors metadata 'metadata<-'
#' @importFrom stats pbinom
#' @export
#' @examples 
#' data("sce_chcl")
#' sce_chcl = cxds(sce_chcl)
cxds <- function( sce, ntop=500, binThresh=0, verb=FALSE, retRes=FALSE){
#===========================================================================================

  #- 1. Binarize counts
  if(verb) cat("\n-> binarizing counts\n")
  B = counts(sce) > binThresh

  #- 2. Select high-variance genes
  if(verb) cat("-> selecting genes\n")
  ps  = Matrix::rowMeans(B)
  vrs = ps*(1-ps)
  hvg = order(vrs,decreasing=TRUE)[1:ntop]
  Bp  = B[hvg,]

  #- score all pairs
  #- Simple model for (01) and (10) observations per gene pair
  if(verb) cat("-> scoring gene pairs\n")
  ps    = Matrix::rowMeans(Bp)
  prb   = outer(ps,1-ps)
  prb   = prb + t(prb)
  obs   = Bp %*% (1-Matrix::t(Bp))
  obs   = obs + Matrix::t(obs)
  #- log p-vals for the observed numbers of 01 and 10
  pvmat = stats::pbinom(as.matrix(obs)-1,prob=prb,size=ncol(Bp),lower.tail=FALSE,log=TRUE)
  #- scores
  if(verb) cat("-> calcluating cell scores\n")
  scores = Matrix::diag(Matrix::t(Bp)%*%pvmat%*%Bp)

  if(retRes){
    if(verb) cat("-> prioritizing gene pairs\n")
    res = list(scores=-scores, S=-pvmat, hvg=hvg, binThresh=binThresh)
    #- rank gene pairs by wighted average doublet contribution
    tmp = Matrix::t(Matrix::t(Bp) *(-scores))
    tmp = tmp %*% Matrix::t(Bp)
    tmp = tmp * pvmat
    pvmat = as.matrix(tmp)

    colnames(pvmat) = 1:ncol(pvmat)
    rownames(pvmat) = colnames(pvmat)
    topPrs = matrix(NA,nrow=100,ncol=2)
    for(i in 1:100){
      pr = which(pvmat == min(pvmat) ,arr.ind=TRUE)[1:2]
      topPrs[i,] = as.integer(colnames(pvmat)[pr])
      pvmat = pvmat[-pr,]
      pvmat = pvmat[,-pr]
    }
    res$topPairs = topPrs

    #- put result in sce
    if(verb) cat("-> done.\n\n")
    hvg_bool                     = (1:nrow(sce)) %in% res$hvg
    hvg_ord                      = rep(NA,nrow(sce))
    hvg_ord[hvg_bool]            = res$hvg
    rowData(sce)$cxds_hvg_bool   = hvg_bool
    rowData(sce)$cxds_hvg_ordr   = hvg_ord
    metadata(sce)$cxds_S         = res$S
    metadata(sce)$cxds_topPairs  = res$topPairs
    metadata(sce)$cxds_binThresh = res$binThresh
    sce$cxds_score               = res$scores
  } else{
    if(verb) cat("-> done.\n\n")
    sce$cxds_score = -scores;
  }
  return(sce)
}

#' Extract top-scoring gene pairs from an SingleCellExperiment where cxds has been run
#'
#' @param sce single cell experiment to analyze; needs "counts" in assays slot.
#' @param n integer. The number of gene pairs to extract. Default: 100
#' @importFrom Matrix t
#' @return matrix Matrix with two colulmns, each containing gene indexes for gene pairs (rows).
cxds_getTopPairs <- function(sce,n=100){
#=======================================

  ind = rowData(sce)$cxds_hvg_ordr[!is.na(rowData(sce)$cxds_hvg_ordr)]
  Bp  = counts(sce)[ind,] > metadata(sce)$cxds_binThresh
  imp = Matrix::t(Matrix::t(Bp) * sce$cxds_score)
  imp = imp  %*% Matrix::t(Bp)
  imp = imp * metadata(sce)$cxds_S
  imp = as.matrix(imp)
  res = list(imp=imp)

  colnames(imp) = 1:ncol(imp)
  rownames(imp) = colnames(imp)
  topPrs = matrix(NA,nrow=n,ncol=2)
  for(i in 1:n){
    pr = which(imp == max(imp) ,arr.ind=TRUE)[1:2]
    topPrs[i,] = as.integer(colnames(imp)[pr])
    imp = imp[-pr,]
    imp = imp[,-pr]
  }
  res$topPairs = topPrs
}
