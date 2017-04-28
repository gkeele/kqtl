#' @export
get.gev.thresholds <- function(min.pval, percentile=0.95){
  evd.pars <- as.numeric(evir::gev(-log10(min.pval))$par.est)
  thresh <- evir::qgev(p=percentile, xi=evd.pars[1], sigma=evd.pars[2], mu=evd.pars[3])
  return(thresh)
}