#' Find doublets/multiplets in UMI scRNA-seq data;
#'
#' Annotates doublets/multiplets using co-expression based approach
#'
#' @param sce single cell experiment (\code{SingleCellExperiment}) object to analyze; needs \code{counts} in assays slot.
#' @param ntop integer, indimessageing number of top variance genes to consider. Default: 500
#' @param binThresh integer, minimum counts to consider a gene "present" in a cell. Default: 0
#' @param verb progress messages. Default: FALSE
#' @param retRes logical, whether to return gene pair scores & top-scoring gene pairs? Default: FALSE.
#' @param estNdbl logical, should the numer of doublets be estimated from the data. Enables doublet calls. Default:FALSE. Use with caution.
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
#' ## create small data set using only 100 cells
#' sce_chcl_small = sce_chcl[, 1:100]
#' sce_chcl_small = cxds(sce_chcl_small)
cxds <- function( sce, ntop=500, binThresh=0, verb=FALSE, retRes=FALSE, estNdbl = FALSE){
#===========================================================================================

  #- check sce argument
  if(!is(sce,"SingleCellExperiment"))
      stop('First argument (sce) needs to be of class SingleCellExperiment')
  if( !("counts" %in% names(assays(sce))) )
      stop('First argument (sce) needs to have a "counts" assay')

  #- 1. Binarize counts
  if(verb) message("\n-> binarizing counts\n")
  B = counts(sce) > binThresh

  #- 2. Select high-variance genes
  if(verb) message("-> selecting genes\n")
  ps  = Matrix::rowMeans(B)
  vrs = ps*(1-ps)
  hvg = order(vrs,decreasing=TRUE)[seq_len(ntop)]
  Bp  = B[hvg,]

  #- score all pairs
  #- Simple model for (01) and (10) observations per gene pair
  if(verb) message("-> scoring gene pairs\n")
  ps    = Matrix::rowMeans(Bp)
  prb   = outer(ps,1-ps)
  prb   = prb + t(prb)
  obs   = Bp %*% (1-Matrix::t(Bp))
  obs   = obs + Matrix::t(obs)
  #- log p-vals for the observed numbers of 01 and 10
  S = stats::pbinom(as.matrix(obs)-1,prob=prb,size=ncol(Bp),
                        lower.tail=FALSE,log=TRUE)
  #- scores
  if(verb) message("-> calcluating cell scores\n")
  scores = Matrix::colSums(Bp * (S%*%Bp))
  sce$cxds_score = -scores
  #- estimate number of doublets
  if(estNdbl){
    if(verb) message("-> estimating number of doublets\n")
    nsamp       = ncol(sce)
    p1          = sample(seq_len(ncol(sce)),nsamp,replace=TRUE)
    p2          = sample(seq_len(ncol(sce)),nsamp,replace=TRUE)
    Bp_sim      = Bp[,p1] + Bp[,p2] #- still logical
    scores_sim  = as.numeric(Matrix::colSums(Bp_sim * (S%*%Bp_sim)))
    est_dbl     = get_dblCalls_ALL(-scores,-scores_sim, rel_loss=1)
    if(is.null(metadata(sce)$cxds)) metadata(sce)$cxds = list()
    metadata(sce)$cxds$ndbl       = est_dbl
    metadata(sce)$cxds$sim_scores = -scores_sim
    sce$cxds_call = sce$cxds_score >= est_dbl["balanced","threshold"]
  }

  if(retRes){
    if(verb) message("-> prioritizing gene pairs\n")
    res = list(scores=-scores, S=-S, hvg=hvg, binThresh=binThresh)
    #- rank gene pairs by wighted average doublet contribution
    tmp = Matrix::t(Matrix::t(Bp) *(-scores))
    tmp = tmp %*% Matrix::t(Bp)
    tmp = tmp * S
    S = as.matrix(tmp)

    colnames(S) = seq_len(ncol(S))
    rownames(S) = colnames(S)
    topPrs = matrix(NA,nrow=100,ncol=2)
    for(i in seq_len(100)){
      pr = which(S == min(S) ,arr.ind=TRUE)[c(1,2)]
      topPrs[i,] = as.integer(colnames(S)[pr])
      S = S[-pr,]
      S = S[,-pr]
    }
    res$topPairs = topPrs

    #- put result in sce
    if(verb) message("-> done.\n\n")
    hvg_bool                     = (seq_len(nrow(sce))) %in% res$hvg
    hvg_ord                      = rep(NA,nrow(sce))
    hvg_ord[hvg_bool]            = res$hvg
    rowData(sce)$cxds_hvg_bool   = hvg_bool
    rowData(sce)$cxds_hvg_ordr   = hvg_ord
    metadata(sce)$cxds$S         = res$S
    metadata(sce)$cxds$topPairs  = res$topPairs
    metadata(sce)$cxds$binThresh = res$binThresh
    #sce$cxds_score               = res$scores
  } else{
    if(verb) message("-> done.\n\n")
    #sce$cxds_score = -scores;
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
  Bp  = counts(sce)[ind,] > metadata(sce)$cxds$binThresh
  imp = Matrix::t(Matrix::t(Bp) * sce$cxds_score)
  imp = imp  %*% Matrix::t(Bp)
  imp = imp * metadata(sce)$cxds$S
  imp = as.matrix(imp)
  res = list(imp=imp)

  colnames(imp) = seq_len(ncol(imp))
  rownames(imp) = colnames(imp)
  topPrs = matrix(NA,nrow=n,ncol=2)
  for(i in seq_len(n)){
    pr = which(imp == max(imp) ,arr.ind=TRUE)[c(1,2)]
    topPrs[i,] = as.integer(colnames(imp)[pr])
    imp = imp[-pr,]
    imp = imp[,-pr]
  }
  res$topPairs = topPrs
}
