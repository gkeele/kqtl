#' Returns a matrix of parametric bootstrap samples from the null model of no locus effect
#'
#' This function takes an scan.h2lmm() object, and returns a specified number of paratmetric bootstrap samples
#' from the null model of no locus effect.
#'
#' @param scan.object A scan.h2lmm() object.
#' @param use.REML DEFAULT: TRUE. Determines whether the variance components for the parametric sampling are 
#' based on maximizing the likelihood (ML) or the residual likelihood (REML).
#' @param num.samples The number of parametric bootstrap samples to return.
#' @param seed DEFAULT: 1. The sampling process is random, thus a seed must be set for samples to be consistent
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
    sim.y.matrix <- matrix(NA, nrow=n, ncol=num.samples)
    
    set.seed(seed)
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
      rownames(sim.y.matrix) <- names(scan.object$fit0$y)
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

#' Returns a matrix of (parametric) permutaions of the outcome
#'
#' This function takes an scan.h2lmm() object, and returns a specified number of permutation samples. When the
#' the null model contains covariates or there is a genetic relationship matrix, permutations are returned based
#' the maintains the structure of the data according to the null model.
#'
#' @param scan.object A scan.h2lmm() object.
#' @param use.REML DEFAULT: TRUE. Determines whether the variance components for the parametric sampling are 
#' based on maximizing the likelihood (ML) or the residual likelihood (REML).
#' @param num.samples The number of permutation samples to return.
#' @param seed DEFAULT: 1. The sampling process is random, thus a seed must be set for samples to be consistent
#' across machines.
#' @export
#' @examples generate.perm.matrix()
generate.perm.matrix <- function(scan.object, use.REML=TRUE, num.samples, seed=1){
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
    perm.y.matrix <- matrix(NA, nrow=n, ncol=num.samples)
    
    set.seed(seed)  
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
      perm.y.ranks <- order(Xb + u + e)
      perm.y.matrix[,i] <- scan.object$fit0$y[perm.y.ranks]
      rownames(perm.y.matrix) <- names(scan.object$fit0$y)
    }
  }
  else{
    stop("Need to add lmer-based functionality!!")
  }
  sim.threshold.object <- list(y.matrix=perm.y.matrix,
                               formula=scan.object$formula,
                               #weights=weights[perm.y.ranks],
                               weights=weights,
                               K=K)
  return(sim.threshold.object)
}

#' Runs threshold scans from a matrix of outcomes, either parametric bootstraps from the null model or permutations
#'
#' This function takes an object produced from either generate.null.bootstrap.matrix() or generate.perm.matrix(), and 
#' runs genome scans on the outcomes contained in them.
#'
#' @param sim.threshold.object An object created by either generate.null.bootstrap.matrix() or generate.perm.matrix().
#' @param keep.full.scans DEFAULT: TRUE. Returns full genome scans for every outcome sample in the sim.threshold.object. Can be used
#' for visualization of the procedure, but greatly increases the size of the output object.
#' @param genomecache The path to the genome cache directory. The genome cache is a particularly structured
#' directory that stores the haplotype probabilities/dosages at each locus. It has an additive model
#' subdirectory and a full model subdirectory. Each contains subdirectories for each chromosome, which then
#' store .RData files for the probabilities/dosages of each locus.
#' @param data A data frame with outcome and potential covariates. Should also have IDs
#' that link to IDs in the genome cache, often the individual-level ID named "SUBJECT.NAME".
#' @param model DEFAULT: additive. Specifies how to model the founder haplotype probabilities. The additive options specifies
#' use of haplotype dosages, and is most commonly used. The full option regresses the phenotype on the actual
#' diplotype probabilities. The diplolasso option specifies the DiploLASSO model.
#' @param use.multi.impute DEFAULT: TRUE. This option specifies whether to use ROP or multiple imputations.
#' @param num.imp DEFAULT: 11. IF multiple imputations are used, this specifies the number of imputations to perform.
#' @param scan.seed DEFAULT: 1. The sampling process is random, thus a seed must be set for samples to be consistent
#' across machines.
#' @export
#' @examples run.threshold.scans()
run.threshold.scans <- function(sim.threshold.object, keep.full.scans=TRUE,
                                genomecache, data,
                                model=c("additive", "full", "diplolasso"),
                                use.multi.impute=TRUE, num.imp=11, 
                                scan.seed=1, ...){
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
    new.y <- data.frame(new.y=y.matrix[,i], SUBJECT.NAME=rownames(y.matrix))
    this.data <- merge(x=new.y, y=data, by="SUBJECT.NAME", all.x=TRUE)

    this.scan <- scan.h2lmm(genomecache=genomecache, data=this.data, formula=iteration.formula, K=K, model=model,
                            use.multi.impute=use.multi.impute, num.imp=num.imp, seed=scan.seed,
                            weights=weights,
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



