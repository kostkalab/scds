

#' Derive doublet calls from classification probabilities
#'
#' Given class probabilities (or scores) discriminating real data from artificial doublets, derive doublet calls. Based on selecting a ROC cutoff, see \emph{The Inconsistency of ‘‘Optimal’’ Cutpoints Obtained using Two Criteria basedon the Receiver Operating Characteristic Curve}, \href{https://dx.doi.org/10.1093/aje/kwj063}{(doi)}.
#'
#' @param scrs_real numeric vector, the scores for the real/original data
#' @param scrs_sim numeric vector, the scores for the artificial doublets
#' @param rel_loss numeric scalar, relative weight of a false positive classification compared with a false negative. Default:1 (same loss for fp and fn).
#' @return numeric, vector containing the (estimated) number of doublets, the score threshold and the fraction of artificial doublets missed (false negative rate, of sorts)
#' @importFrom pROC roc
get_dblCalls_ROC <- function(scrs_real, scrs_sim, rel_loss=1){
#=============================================================

    p    = length(scrs_sim)/(length(scrs_sim)+length(scrs_real))
    rc   = roc(response=c(rep(0,length(scrs_real)),rep(1,length(scrs_sim))),predictor=c(scrs_real,scrs_sim))
    sens = rc$sensitivities
    spec = rc$specificities
    #- use youden to find ROC cutoff
    r       = p/rel_loss/(1-p)
    cut_ind = which.max(sens+r*spec-1)
    thresh  = rc$thresholds[cut_ind]

    ndbl    = sum(scrs_real >= thresh)
    fram    = mean(scrs_sim < thresh) #- fraction of sim doublets missed

    res = c(ndbl,thresh,fram)
    names(res) = c("number","threshold","fnr") #- false negative rate
    return(res)
}

#' Derive doublet calls from doublset scores
#'
#' Given score vectors for real data and artificial doubles, derive doublet calls based on determining doublet score cutoffs.
#'
#' @param scrs_real numeric vector, the scores for the real/original data
#' @param scrs_sim numeric vector, the scores for the artificial doublets
#' @param type character or numeric, describes how the score threshold for calling doublets is determined. Either \code{"balanced"} or a number between zero and one that indicates the fraction of artificial doublets missed when making calls. Default: \code{"balanced"}.
#' @return numeric, vector containing the (estimated) number of doublets, the score threshold and the fraction of artificial doublets missed (false negative rate, of sorts)
#' @importFrom stats optimize ecdf uniroot quantile
get_dblCalls_dist <- function(scrs_real,scrs_sim, type="balanced"){
#==================================================================

  #- do "balanced errpr"

  if(type == "balanced"){
  #======================

    es  = ecdf(scrs_sim)
    er  = ecdf(scrs_real)
    rtf = function(thresh) 1-er(thresh) #- right tail; decreases mono with arg
    ltf = function(thresh) es(thresh)   #- left tail; increases mono with arg

    res        = uniroot(function(x)ltf(x) - rtf(x), c(min(scrs_real),max(scrs_real)))
    res_val    = res$root
  } else {
    if(!is.numeric(type)) stop("invalid type argument\n")
    res_val = quantile(scrs_sim,prob=type)
  }
  res_ndl    = sum(scrs_real >= res_val)
  res_fnr    = mean(scrs_sim < res_val)
  res        = c(res_ndl,res_val,res_fnr)
  names(res) = c("number","threshold","fnr") #- false negative rate
  return(res)
}

#' Wrapper for getting doublet calls
#'
#' @param scrs_real numeric vector, the scores for the real/original data
#' @param scrs_sim numeric vector, the scores for the artificial doublets
#' @param rel_loss numeric scalar, relative weight of a false positive classification compared with a false negative. Default:1 (same loss for fp and fn).
#' @return numeric, matrix containing the (estimated) number of doublets, the score threshold and the fraction of artificial doublets missed (false negative rate, of sorts) as columns and four types of estimating: "youden", "balanced" and a false negative rate of artificial doublets of 0.1 and 0.01, respecitvely.
get_dblCalls_ALL <- function(scrs_real,scrs_sim,rel_loss=1){
#=========================================================
  est_dbl = rbind( get_dblCalls_ROC( scrs_real,scrs_sim,rel_loss),
                   get_dblCalls_dist(scrs_real,scrs_sim,"balanced"),
                   get_dblCalls_dist(scrs_real,scrs_sim,0.1),
                   get_dblCalls_dist(scrs_real,scrs_sim,0.01))
  rownames(est_dbl)  = c("youden","balanced","0.1","0.01")
  return(est_dbl)
}
