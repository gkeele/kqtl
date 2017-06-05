#' Returns a significance threshold based on fitting max LODs or max -log10p to a generalized extreme value (GEV) distribution
#'
#' This function takes an scan.h2lmm() object, and returns a specified number of outcome samples, either permutations or
#' from the null model of no locus effect.
#'
#' @param extreme.statistic A vector of max LODs or min p-values from null scans of the data.
#' @param use.lod DEFAULT: FALSE. "FALSE" specifies LOD scores. "TRUE" specifies p-values.
#' @param percentile DEFAULT: 0.95. The desired alpha level (false positive probability) from the GEV distribution.
#' @export
#' @examples get.gev.thresholds()
get.gev.thresholds <- function(extreme.statistic, use.lod=FALSE, percentile=0.95){
  if(!use.lod){
    evd.pars <- as.numeric(evir::gev(-log10(extreme.statistic))$par.est)
  }
  if(use.lod){
    evd.pars <- as.numeric(evir::gev(extreme.statistic)$par.est)
  }
  thresh <- evir::qgev(p=percentile, xi=evd.pars[1], sigma=evd.pars[2], mu=evd.pars[3])
  return(thresh)
}