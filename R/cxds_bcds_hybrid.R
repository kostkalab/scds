#' Find doublets/multiples in UMI scRNA-seq data;
#'
#' Annotates doublets/multiplets using the hybrid approach
#'
#' @param sce single cell experiment (\code{SingleCellExperiment}) object to analyze; needs \code{counts} in assays slot.
#' @param cxdsArgs list, arguments for cxds function in list form. Default: NULL
#' @param bcdsArgs list, arguments for bcds function in list form. Default: NULL
#' @param verb logical, switch on/off progress messages
#' @param force logical, force a (re)run of \code{cxds} and \code{bcds}. Default: FALSE
#' @param estNdbl logical, should the numer of doublets be estimated from the data. Enables doublet calls. Default:FALSE. Use with caution.
#' @return sce input sce object \code{SingleCellExperiment} with doublet scores added to colData as "hybrid_score" column.
#' @importFrom dplyr %>%
#' @import SingleCellExperiment
#' @export
#' @examples
#' data("sce_chcl")
#' ## create small data set using only 100 cells
#' sce_chcl_small = sce_chcl[, 1:100]
#' sce_chcl_small = cxds_bcds_hybrid(sce_chcl_small)
cxds_bcds_hybrid <- function(sce, cxdsArgs=NULL, bcdsArgs=NULL, verb=FALSE, estNdbl=FALSE,force=FALSE){
#==============================================================================

  #- enforce consistency wrt doublet calls
  cxdsArgs$estNdbl = estNdbl
  bcdsArgs$estNdbl = estNdbl

  #- re-run only if scores are missing or estNdbl has not been performed.
  #  FIXME: cxdsArgs and bcdsArgs might be inconsistent between what was passed to this function and what was used in the \code{sce} that is passed.

  run_cxds  = is.null(colData(sce)$cxds_score)
  run_bcds  = is.null(colData(sce)$bcds_score)
  run_cxds  = run_cxds | (estNdbl & (is.null(metadata(sce)$cxds$ndbl)))
  run_bcds  = run_bcds | (estNdbl & (is.null(metadata(sce)$bcds$ndbl)))

  if(force) {
    run_bcds = TRUE
    run_cxds = TRUE
  }

  if(run_cxds) {
    if(verb) message("-> running cxds\n")
    sce = do.call(cxds,c(list(sce=sce),cxdsArgs))
  }
  if(run_bcds) {
    if(verb) message("-> running bcds\n")
    sce = do.call(bcds,c(list(sce=sce),bcdsArgs))
  }
  #- sqish and average scores
  squish  = function(x) { x = x-min(x) ; x = x/max(x) ; return(x) }

  s.cxds   = sce$cxds_score %>% squish
  s.bcds   = sce$bcds_score %>% squish
  s.hybrid = s.cxds + s.bcds
  #- put in sce
  sce$hybrid_score = s.hybrid

  if(estNdbl){
    #- hybrid score of artificial doublets
    s.cxds.sim = metadata(sce)$cxds$sim_scores
    s.bcds.sim = metadata(sce)$bcds$sim_scores
    if(length(s.cxds.sim) > length(s.bcds.sim)){
        s.cxds.sim = s.cxds.sim[sample(seq_len(length(s.bcds.sim)))]
    }
    if(length(s.cxds.sim) < length(s.bcds.sim)){
        s.bcds.sim = s.bcds.sim[sample(seq_len(length(s.cxds.sim)))]
    }
    s.hybrid.sim = (s.cxds.sim %>% squish) + (s.bcds.sim %>% squish)

    #- annotate
    est_dbl = get_dblCalls_ALL(sce$hybrid_score,s.hybrid.sim)
    if(is.null(metadata(sce)$hybrid)) metadata(sce)$hybrid = list()
    metadata(sce)$hybrid$ndbl = est_dbl
    metadata(sce)$hybrid$sim_scores = s.hybrid.sim
    sce$hybrid_call = sce$hybrid_score >= est_dbl["balanced","threshold"]
  }

  return(sce)
}
