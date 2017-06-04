#' Returns a matrix of outcome samples, either permutations or from the null model of no locus effect
#'
#' This function takes an scan.h2lmm() object, and returns a specified number of outcome samples, either permutations or
#' from the null model of no locus effect.
#'
#' @param scan.object A scan.h2lmm() object.
#' @param method DEFAULT: "bootstrap". "bootstrap" specifies parametric bootstraps from the null model. "permutation" specifies
#' parametric permutations that respect the structure of the data. Permutations are more appropriate if the data have highly
#' influential data points.
#' @param use.REML DEFAULT: TRUE. Determines whether the variance components for the parametric sampling are 
#' based on maximizing the likelihood (ML) or the residual likelihood (REML).
#' @param num.samples The number of parametric bootstrap samples to return.
#' @param seed DEFAULT: 1. The sampling process is random, thus a seed must be set for samples to be consistent
#' across machines.
#' @export
#' @examples generate.sample.outcomes.matrix()
generate.sample.outcomes.matrix <- function(scan.object, model.type=c("null", "alt"), 
                                            method=c("bootstrap", "permutation"), use.REML=TRUE, 
                                            use.BLUP=FALSE, num.samples, seed=1){
  model.type <- model.type[1]
  method <- method[1]
  
  if(model.type == "null"){ fit <- scan.object$fit0; locus <- NULL }
  if(model.type == "alt"){ fit <- scan.object$fit1; locus <- scan.object$loci }
  fit0.REML <- scan.object$fit0.REML
  if(class(fit) != "lmerMod"){
    Xb <- fit$x %*% fit$coefficients
    n <- nrow(fit$x)
    K <- fit$K
    weights <- fit$weights
    return.weights <- weights
    if(is.null(weights)){ 
      weights <- rep(1, nrow(K)) 
    }
    if(use.REML){
      if(is.null(K)){
        tau2 <- 0
        sigma2 <- fit$sigma2.mle*(n/(n - 1))
      }
      else{
        tau2 <- fit0.REML$tau2.mle
        sigma2 <- fit0.REML$sigma2.mle
      }
    }
    else{
      tau2 <- fit$tau2.mle
      sigma2 <- fit$sigma2.mle  
    }
    sim.y.matrix <- matrix(NA, nrow=n, ncol=num.samples)
    
    if(is.null(K)){
      u <- rep(0, n)
    }
    else{
      original.K <- K
      impute.map <- scan.object$impute.map
      K <- reduce.large.K(large.K=K, impute.map=impute.map)
      if(use.BLUP){
        X <- fit$x
        Sigma <- original.K*tau2 + diag(1/weights)*sigma2
        inv.Sigma <- solve(Sigma)
        u.BLUP <- (original.K*tau2) %*% inv.Sigma %*% (diag(nrow(original.K)) - X %*% solve(t(X) %*% inv.Sigma %*% X) %*% t(X) %*% inv.Sigma) %*% fit$y  
      }
    }
    
    set.seed(seed)
    for(i in 1:num.samples){
      if(!is.null(K)){
        ## Handling potential replicates
        if(use.BLUP){
          u <- u.BLUP
        }
        else{
          u <- c(mnormt::rmnorm(1, mean=rep(0, nrow(K)), varcov=K*tau2))
        }
        names(u) <- unique(impute.map[,2])
        u <- u[impute.map[,2]]
      }
      if(is.null(weights)){
        e <- rnorm(n=n, mean=0, sd=sqrt(sigma2))
      }
      else{
        e <- c(mnormt::rmnorm(1, mean=rep(0, n), varcov=diag(1/weights)*sigma2))
      }
      y.sample <- Xb + u + e
      if(method == "bootstrap"){
        sim.y.matrix[,i] <- y.sample
      }
      if(method == "permutation"){
        perm.y.ranks <- order(y.sample)
        sim.y.matrix[,i] <- fit$y[perm.y.ranks]
      }
    }
    rownames(sim.y.matrix) <- names(fit$y)
  }
  else{
    stop("Need to add lmer-based functionality!!")
  }
  sim.threshold.object <- list(y.matrix=sim.y.matrix,
                               formula=scan.object$formula,
                               weights=return.weights,
                               K=K,
                               method=method,
                               impute.map=scan.object$impute.map,
                               locus=locus)
  return(sim.threshold.object)
}

reduce.large.K <- function(large.K, impute.map){
  map.order <- match(impute.map[,1], table=colnames(large.K))
  impute.map <- impute.map[map.order,]
  colnames(large.K) <- rownames(large.K) <- impute.map[,2]
  K <- large.K[unique(as.character(impute.map[,2])), unique(as.character(impute.map[,2]))]
  return(K)
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
#' @param chr DEFAULT: "all". The chromosomes to conduct scans over.
#' @param scan.seed DEFAULT: 1. The sampling process is random, thus a seed must be set for samples to be consistent
#' across machines.
#' @export
#' @examples run.threshold.scans()
run.threshold.scans <- function(sim.threshold.object, keep.full.scans=TRUE,
                                genomecache, data,
                                model=c("additive", "full", "diplolasso"),
                                use.multi.impute=TRUE, num.imp=11, chr="all", 
                                scan.seed=1, ...){
  y.matrix <- sim.threshold.object$y.matrix
  formula <- sim.threshold.object$formula
  weights <- sim.threshold.object$weights
  K <- sim.threshold.object$K
  pheno.id <- names(sim.threshold.object$impute.map)[1]
  geno.id <- names(sim.threshold.object$impute.map)[2]
  
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
  
  iteration.formula <- formula(paste0("new_y ~ ", unlist(strsplit(formula, split="~"))[-1]))
  for(i in 1:num.scans){
    new.y <- data.frame(y.matrix[,i], rownames(y.matrix))
    names(new.y) <- c("new_y", pheno.id)
    this.data <- merge(x=new.y, y=data, by=pheno.id, all.x=TRUE)

    this.scan <- scan.h2lmm(genomecache=genomecache, data=this.data, 
                            formula=iteration.formula, K=K, model=model,
                            use.multi.impute=use.multi.impute, num.imp=num.imp, 
                            pheno.id=pheno.id, geno.id=geno.id, seed=scan.seed,
                            weights=weights, chr=chr,
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



