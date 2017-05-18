#' Returns a matrix of parametric bootstrap samples from the null model of no locus effect
#'
#' This function takes an scan.h2lmm() object, and returns a specified number of paratmetric bootstrap samples
#' from the null model of no locus effect.
#'
#' @param scan.object A scan.h2lmm() object.
#' @param use.REML DEFAULT: TRUE. Determines whether the variance components for the parametric sampling are 
#' based on maximizing the likelihood (ML) or the residual likelihood (REML).
#' @param num.samples The number of parametric bootstrap samples to return.
#' @param seed DEFAULT: 1. The sampling process is random, thus a seed must be set for samples to be the same
#' across machines.
#' @export
#' @examples generate.null.bootstrap.matrix()
generate.null.bootstrap.matrix <- function(scan.object, use.REML=TRUE, num.samples, seed=1){
  if(class(scan.object$fit0) != "lmerMod"){
    Xb <- scan.object$fit0$x %*% scan.object$fit0$coefficients
    n <- nrow(scan.object$fit0$x)
    K <- scan.object$fit0$K
    weights <- scan.object$fit0$weights
    if(use.REML){
      if(is.null(K)){
        tau2 <- 0
        sigma2 <- scan.object$fit0$sigma2.mle*(n/(n - 1))
      }
      else{
        tau2 <- scan.object$fit0.REML$tau2.mle
        sigma2 <- scan.object$fit0.REML$sigma2.mle
      }
    }
    else{
      tau2 <- scan.object$fit0$tau2.mle
      sigma2 <- scan.object$fit0$sigma2.mle  
    }
    set.seed(seed)
    
    sim.y.matrix <- matrix(NA, nrow=n, ncol=num.samples)
    for(i in 1:num.samples){
      if(is.null(K)){
        u <- rep(0, n)
      }
      else{
        u <- c(mnormt::rmnorm(1, mean=rep(0, n), varcov=K*tau2))
      }
      if(is.null(weights)){
        e <- rnorm(n=n, mean=0, sd=sqrt(sigma2))
      }
      else{
        e <- c(mnormt::rmnorm(1, mean=rep(0, n), varcov=diag(1/weights)*sigma2))
      }
      sim.y.matrix[,i] <- Xb + u + e
    }
  }
  else{
    stop("Need to add lmer-based functionality!!")
  }
  
  sim.threshold.object <- list(y.matrix=sim.y.matrix,
                               formula=scan.object$formula,
                               weights=weights,
                               K=K)
  return(sim.threshold.object)
}

#' @export
generate.parametric.perm.matrix <- function(scan.object, use.REML, num.samples, seed=1){
  if(use.REML){
    null.fit <- scan.object$fit0.REML
  }
  else{
    null.fit <- scan.object$fit0
  }
  set.seed(seed)
  
  Xb <- null.fit$x %*% null.fit$coefficients
  perm.y.matrix <- matrix(NA, nrow=length(Xb), ncol=num.samples)
  for(i in 1:num.samples){
    u <- c(mnormt::rmnorm(1, mean=rep(0, nrow(null.fit$x)), varcov=null.fit$K*null.fit$tau2.mle))
    if(is.null(null.fit$weights)){
      e <- rnorm(n=nrow(null.fit$x), mean=0, sd=sqrt(null.fit$sigma2.mle))
    }
    else{
      e <- c(mnormt::rmnorm(1, mean=rep(0, nrow(null.fit$x)), varcov=diag(1/null.fit$weights)*null.fit$sigma2.mle))
    }
    perm.y.ranks <- order(Xb + u + e)
    perm.y.matrix[,i] <- null.fit$y[perm.y.ranks]
  }
  sim.threshold.object <- list(y.matrix=perm.y.matrix,
                               formula=scan.object$formula,
                               weights=null.fit$weights[perm.y.ranks], # Keeping the weights matched with outcomes
                               K=null.fit$K)
  return(sim.threshold.object)
}

#' @export
run.threshold.scans <- function(sim.threshold.object, keep.full.scans=TRUE,
                                genomecache, data,
                                model=c("additive", "full", "diplolasso"),
                                use.par="h2", use.multi.impute=TRUE, num.imp=10, brute=TRUE, use.fix.par=FALSE, 
                                scan.seed=1, do.augment=FALSE, chr="all",
                                use.augment.weights=FALSE, use.full.null=FALSE, added.data.points=1, ...){
  y.matrix <- sim.threshold.object$y.matrix
  formula <- sim.threshold.object$formula
  weights <- sim.threshold.object$weights
  K <- sim.threshold.object$K
  
  num.scans <- ncol(y.matrix)
  
  h <- DiploprobReader$new(genomecache)
  loci <- h$getLoci()
  loci.chr <- h$getChromOfLocus(loci)
  if(chr != "all"){
    loci.chr <- h$getChromOfLocus(loci)
    loci <- loci[loci.chr %in% chr]
  }
  
  full.results <- these.chr <- these.pos <- NULL
  if(keep.full.scans){
    full.results <- matrix(NA, nrow=num.scans, ncol=length(loci))
    colnames(full.results) <- loci
    these.chr <- h$getChromOfLocus(loci)
    these.pos <- list(Mb=h$getLocusStart(loci, scale="Mb"),
                      cM=h$getLocusStart(loci, scale="cM"))
  }
  max.results <- rep(NA, num.scans)
  
  iteration.formula <- formula(paste0("new.y ~ ", unlist(strsplit(formula, split="~"))[-1]))
  for(i in 1:num.scans){
    new.y <- data.frame(new.y=y.matrix[,i], SUBJECT.NAME=colnames(K))
    this.data <- merge(x=new.y, y=data, by="SUBJECT.NAME", all.x=TRUE)
    
    this.scan <- scan.h2lmm(genomecache=genomecache, data=this.data, formula=iteration.formula, K=K, model=model,
                            use.par=use.par, use.multi.impute=use.multi.impute, num.imp=num.imp, chr=chr, brute=brute, use.fix.par=use.fix.par, seed=scan.seed, do.augment=do.augment, 
                            weights=weights, use.augment.weights=use.augment.weights, use.full.null=use.full.null, added.data.points=added.data.points,
                            ...)
    if(keep.full.scans){
      full.results[i,] <- this.scan$p.value
    }
    max.results[i] <- min(this.scan$p.value)
    cat("threshold scan:", i, "\n")
  }
  return(list(full.results=list(p.values=full.results, 
                                chr=these.chr, 
                                pos=these.pos), 
              max.statistics=max.results))
}



