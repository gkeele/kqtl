#' @export
scan.EM <- function(genomecache, data, formula, model="additive",
                    chr="all", 
                    seed=1, oracle=FALSE,
                    print.locus.finished=TRUE, convergence.limit=0.0001, step.limit=1000
                    do.augment=FALSE, just.these.loci=NULL, print.locus.fit=FALSE, ...){
 
  h <- DiploprobReader$new(genomecache)
  founders <- h$getFounders()
  num.founders <- length(founders)
  loci <- h$getLoci()
  cache.subjects <- rownames(h$getLocusMatrix(loci[1], model="additive"))
  
  data.and.K <- make.processed.data(formula=formula, data=data, cache.subjects=cache.subjects, K=NULL, impute.on="SUBJECT.NAME")
  data <- data.and.K$data

  loci.chr <- h$getChromOfLocus(loci)
  if(chr != "all"){
    loci.chr <- h$getChromOfLocus(loci)
    loci <- loci[loci.chr %in% chr]
  }
  if(!is.null(just.these.loci)){
    loci <- loci[loci %in% just.these.loci]
    loci.chr <- loci.chr[loci %in% just.these.loci]
  }
  
  augment.indicator <- NULL
  formula.string <- Reduce(paste, deparse(formula))
  null.formula <- make.null.formula(formula, do.augment=do.augment)
  original.n <- nrow(data)
  old.data <- data
  
  ###### Augmentation
  if(do.augment){
    augment.n <- ifelse(model=="additive", num.founders, num.founders + choose(num.founders, 2))
    data <- make.simple.augment.data(data=data, augment.n=augment.n)
  }
  
  ###### Null model
  ## No kinship effect - weights or no weights
  fit0 <- lmmbygls(null.formula, data=data, eigen.K=NULL, K=NULL, use.par="h2", fix.par=0, weights=NULL, brute=FALSE)
  
  LOD.vec <- p.vec <- df <- rep(NA, length(loci))
  null.data <- data
  
  for(i in 1:length(loci)){
    X.add <- h$getLocusMatrix(loci[i], model="additive", subjects=data$SUBJECT.NAME[1:original.n])
    colnames(X.add) <- gsub(pattern="/", replacement=".", x=colnames(X.add), fixed=TRUE)
    locus.formula <- make.EM.alt.formula(formula=formula, X=X.add)
    
    X.full <- h$getLocusMatrix(loci[i], model="full", subjects=data$SUBJECT.NAME[1:original.n])
    colnames(X.full) <- gsub(pattern="/", replacement=".", x=colnames(X.full), fixed=TRUE)
    
    if(do.augment){
      X.names <- rownames(X.add)
      if(model=="additive"){
        X.add <- rbind(X.add, 2*diag(augment.n))
      }
      if(model=="full"){
        X.add <- rbind(X.add, diag(augment.n))
      }
      rownames(X.add) <- c(X.names, paste0("augment.obs", 1:augment.n))
    }
      
    ROP.data <- cbind(null.data, X.add)
    ## Using ROP to get starting values
    ROP.results <- lmmbygls(formula=locus.formula, data=ROP.data, 
                            eigen.K=NULL, K=NULL,
                            use.par="h2", fix.par=0)
    coefficients <- ROP.results$coefficients
    coefficients[is.na(coefficients)] <- 0
    sigma2 <- ROP.results$sigma2.mle
    
    expanded.diplotype.data <- get.expanded.X(X.full, coefficients=NULL)
    expanded.dosage.data <- rotate.full.to.add.data(diplotype.data=expanded.diplotype.data)
    X <- as.matrix(expanded.dosage.data[,-1])
    colnames(X) <- names(expanded.dosage.data)[-1]
    rownames(X) <- expanded.dosage.data[,1]
    X <- apply.contrast(X)
    expanded.y <- get.expanded.y(y=data[,all.vars(formula)[1]], num.groups=36)
    start.weights <- get.weights.from.diplotype.prob.data(X.full)
    
    steps <- 0
    has.converged <- FALSE
    max.dif.vec <- converge.vec <- NULL
    set.seed(seed)
    while(!has.converged){
      previous.coefficients <- coefficients
      previous.sigma2 <- sigma2
      weights <- E.step(start.weights=start.weights, y=expanded.y, X=X, coefficients=coefficients, sigma2=sigma2)
      WLS.fit <- M.step(weights=weights, y=expanded.y, X=X)
      coefficients <- WLS.fit$coefficients
      sigma2 <- WLS.fit$sigma2
      
      coefficient.dif <- get.coef.dif(coefficients=coefficients, previous.coefficients=previous.coefficients)
      max.dif <- max(abs(c(coefficient.dif, sigma2 - previous.sigma2)))
      max.dif.vec <- c(max.dif.vec, max.dif)
      if(max.dif < convergence.limit){ 
        has.converged <- TRUE 
        converge.vec <- c(converge.vec, TRUE)
      }
      steps <- steps + 1
      
      cat(paste0("step: ", steps, ", max diff: ", max.dif, "\n"))
      if(steps == step.limit){
        has.converged <- TRUE
        converge.vec <- c(converge.vec, FALSE)
      }
    }
    alt.logLik <- get.mixture.likelihood(start.weights=start.weights, expanded.X=X, 
                                         expanded.y=expanded.y, coefficients=coefficients, sigma2=sigma2)
    im.rank <- sum(!is.na(coefficients))
    
    df[i] <- im.rank
    LOD.vec[i] <- log10(exp(alt.logLik - fit0$logLik))
    p.vec[i] <- pvalue.per.locus.im(im.logLik=alt.logLik, im.rank=im.rank, fit0=fit0)
    if(print.locus.fit){ cat(paste("locus", i, "out of", length(loci)), "\n") }
  }
  names(LOD.vec) <- names(p.vec) <- names(df) <- loci
  output <- list(LOD=LOD.vec,
                 p.value=p.vec,
                 MI.LOD=NULL,
                 MI.p.value=NULL,
                 df=df,
                 pos=list(Mb=h$getMarkerLocation(loci, scale="Mb"), cM=h$getMarkerLocation(loci, scale="cM")),
                 loci=loci, 
                 chr=h$getChromOfLocus(loci),
                 fit0=fit0,
                 fit0.REML=NULL,
                 y=data$y,
                 formula=formula.string,
                 model.type=model)
  return(output)
}

rotate.full.to.add.data <- function(diplotype.data){
  rotate.diplotypes.to.dosages <- straineff.mapping.matrix()
  dosage.data <- as.matrix(diplotype.data[,-1]) %*% rotate.diplotypes.to.dosages
  dosage.data <- data.frame(SUBJECT.NAME=diplotype.data$SUBJECT.NAME, dosage.data)
  names(dosage.data)[-1] <- LETTERS[1:8]
  return(dosage.data)
}
get.weights.from.diplotype.prob.data <- function(diplotype.data){
  probs.diplotype <- as.matrix(diplotype.data[,-1])
  weights <- as.vector(t(probs.diplotype))
  return(weights)
}
get.expanded.X <- function(diplotype.data, coefficients=NULL){
  if(!is.null(coefficients)){
    diplotype.data <- diplotype.data[,c("SUBJECT.NAME", names(coefficients))]
  }
  
  n <- nrow(diplotype.data)
  p <- ncol(diplotype.data) - 1
  X <- matrix(rep(diag(p), n), ncol=p, byrow=TRUE)
  colnames(X) <- names(diplotype.data)[-1]
  
  expand.subjects <- paste(rep(diplotype.data[,1], each=p), rep(1:p, p), sep=".")
  expand.data <- data.frame(SUBJECT.NAME=expand.subjects, X)
  return(expand.data)
}
get.expanded.y <- function(y, num.groups){
  y <- rep(y, each=num.groups)
  return(y)
}
apply.contrast <- function(X){
  C <- create.treatment.contrast.matrix(num.par=ncol(X), remove.column.index=ncol(X))
  new.X <- cbind(rep(1, nrow(X)), X %*% C)
  colnames(new.X) <- c("(Intercept)", colnames(X)[-ncol(X)])
  return(new.X)
}
create.treatment.contrast.matrix <- function(num.par, remove.column.index){
  C <- contr.treatment(n=num.par, base=remove.column.index)
  return(C)
}
E.step <- function(start.weights, y, X, coefficients, sigma2, p=36){
  mu <- X %*% coefficients
  pheno.prob.matrix <- matrix(dnorm(y, mean=mu, sd=sqrt(sigma2)), ncol=p, byrow=TRUE)
  dip.prob.matrix <- matrix(start.weights, ncol=p, byrow=TRUE)
  product.matrix <- dip.prob.matrix * pheno.prob.matrix
  scale.matrix <- matrix(rep(rowSums(product.matrix), each=p), ncol=p, byrow=TRUE)
  weight.matrix <- product.matrix/scale.matrix
  weights <- as.vector(t(weight.matrix))
  return(weights)
}
M.step <- function(weights, y, X, p=36){
  fit <- lm.wfit(x=X, y=y, w=weights)
  # fit$residuals appear to be not be weighted
  #fit$sigma2 <- sum((fit$residuals^2))/sum(weights)
  fit$sigma2 <- sum(fit$residuals * weights * fit$residuals)/sum(weights)
  return(fit)
}
remove.weights <- function(weights, previous.coefficients, coefficients){
  weight.matrix <- matrix(weights, ncol=length(previous.coefficients), byrow=TRUE)
  remove.col <- is.na(coefficients[names(previous.coefficients)])
  if(any(remove.col)){
    weight.matrix <- weight.matrix[,-which(remove.col)]
  }
  weights <- as.vector(t(weight.matrix))
  return(weights)
}
reshape.X <- function(expanded.X, previous.coefficients, coefficients, n){
  remove.col.logical <- is.na(coefficients[names(previous.coefficients)])
  remove.col <- which(remove.col.logical)
  if(any(remove.col.logical)){
    remove.row <- NULL
    for(i in 1:length(remove.col)){
      remove.row <- c(remove.row, (0:(n-1))*length(previous.coefficients) + remove.col[i])
    }
    remove.row <- sort(remove.row)
    expanded.X <- expanded.X[-remove.row, -remove.col]
  }
  return(expanded.X)
}
get.coef.dif <- function(coefficients, previous.coefficients){
  coef.dif <- coefficients[names(previous.coefficients)] - previous.coefficients
  coef.dif <- coef.dif[!(is.na(coef.dif) | is.nan(coef.dif))]
  return(coef.dif)
}
get.mixture.likelihood <- function(start.weights, expanded.X, expanded.y, coefficients, sigma2, p=36){
  #coefficients[is.na(coefficients)] <- 0
  mu <- expanded.X %*% coefficients
  #p <- length(coefficients)
  pheno.prob.matrix <- matrix(dnorm(expanded.y, mean=mu, sd=sqrt(sigma2)), ncol=p, byrow=TRUE)
  dip.prob.matrix <- matrix(start.weights, ncol=p, byrow=TRUE)
  product.matrix <- dip.prob.matrix * pheno.prob.matrix
  logLik <- sum(log(rowSums(product.matrix)))
  return(logLik)
}
pvalue.per.locus.im <- function(im.logLik, im.rank, fit0){
  stat <- -2*(fit0$logLik - im.logLik)
  stat <- ifelse(stat < 0, 0, stat)
  pvalue <- pchisq(q=stat, df=im.rank-fit0$rank, lower.tail=FALSE)
  return(pvalue)
}


