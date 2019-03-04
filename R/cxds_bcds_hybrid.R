#' Find doublets/multiples in UMI scRNA-seq data;
#'
#' Annotates doublets/multiplets using the hybrid approach
#' 
#' @param sce single cell experiment (\code{SingleCellExperiment}) object to analyze; needs \code{counts} in assays slot.
#' @param cxds_args list, arguments for cxds function in list form. Default: NULL
#' @param bcds_args list, arguments for bcds function in list form. Default: NULL
#' @return sce input sce object \code{SingleCellExperiment} with doublet scores added to colData as "hybrid_score" column.
#' @import SingleCellExperiment
#' @export
#' @examples 
#' data("sce")
#' sce = cxds_bcds_hybrid(sce)
cxds_bcds_hybrid <- function(sce, cxds_args=NULL, bcds_args=NULL){
#=================================================================

  #- check if we need to run cxds/bcds
  if(is.null(colData(sce)$cxds_score)) sce = do.call(cxds,c(list(sce=sce),cxds_args))
  if(is.null(colData(sce)$bcds_score)) sce = do.call(bcds,c(list(sce=sce),bcds_args))

  #- sqish and average scores
  s.cxds = sce$cxds_score - min(sce$cxds_score)
  s.bcds = sce$bcds_score - min(sce$bcds_score)
  s.cxds = s.cxds/max(s.cxds)
  s.bcds = s.bcds/max(s.bcds)

  #- put in sce
  sce$hybrid_score = s.cxds + s.bcds

  return(sce)
}
