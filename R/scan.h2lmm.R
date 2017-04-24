#' @export
scan.h2lmm <- function(genomecache, data, formula, K=NULL,
                       model=c("additive", "full", "diplolasso"), diplolasso.refit=FALSE,
                       use.par="h2", use.multi.impute=TRUE, num.imp=10, chr="all", brute=TRUE, use.fix.par=TRUE, 
                       seed=1, impute.on="SUBJECT.NAME",
                       weights=NULL, do.augment=FALSE, use.augment.weights=FALSE, use.full.null=FALSE, added.data.points=1, 
                       just.these.loci=NULL,
                       print.locus.fit=FALSE,
                       ...){ # diplolasso and diplolasso refit not work correctly
  model <- model[1]
  # Defaults if DiploLASSO is not used
  fit0.glmnet <- diplolasso.penalty.factor <- NULL
  
  h <- bagpipe.backend::DiploprobReader$new(genomecache)
  founders <- h$getFounders()
  num.founders <- length(founders)
  loci <- h$getLoci()
  cache.subjects <- rownames(h$getLocusMatrix(loci[1], model="additive"))
  
  data.and.K <- make.processed.data(formula=formula, data=data, cache.subjects=cache.subjects, K=K, impute.on=impute.on)
  data <- data.and.K$data
  K <- data.and.K$K
  if(!is.null(weights)){ weights <- weights[data$SUBJECT.NAME] }

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
    augment.indicator <- c(rep(0, original.n), rep(1, augment.n))
    if(!use.full.null){
      data <- make.simple.augment.data(data=data, augment.n=augment.n)
      K <- make.simple.augment.K(K=K, augment.n=augment.n)
    }
    if(use.full.null){
      no.augment.K <- K
      K <- make.full.null.augment.K(K=no.augment.K, original.n=original.n, augment.n=augment.n)
      data <- make.full.null.augment.data(formula=formula, data=data, no.augment.K=no.augment.K, use.par=use.par, brute=brute,
                                          original.n=original.n, augment.n=augment.n, weights=weights)
    }
    if(use.augment.weights){
      weights <- make.augment.weights(data=data, augment.n=augment.n, added.data.points=added.data.points)
    }
  }
  
  ###### Null model
  ## check for LMER notation
  use.lmer <- check.for.lmer.formula(null.formula)
  if(use.lmer & !is.null(K)){
    stop("Cannot use LMER sparse random effects AND a non-sparse random effect", call.=FALSE)
  }
  if(use.lmer & model == "diplolasso"){
    stop("Cannot use DiploLASSO with LMER currently", call.=FALSE)
  }
  
  if(use.lmer){
    fit0 <- lmmbylmer(formula=null.formula, data=data, REML=FALSE, weights=weights)
    fit0.REML <- lmmbylmer(formula=null.formula, data=data, REML=TRUE, weights=weights)
    fix.par <- NULL
  }
  else{
    ## No kinship effect - weights or no weights
    if(is.null(K)){
      fit0 <- lmmbygls(null.formula, data=data, eigen.K=NULL, K=NULL, use.par="h2", fix.par=0, weights=weights, brute=brute)
      fit0.REML <- lmmbygls(null.formula, data=data, eigen.K=NULL, K=NULL, use.par="h2.REML", fix.par=0, weights=weights, brute=brute)
    }
    ## Kinship effect - weights or no weights
    else{
      ###### This function is made for handling constant weights at all loci
      if(!is.null(weights)){
        J <- weights^(-1/2) * t(weights^(-1/2) * K)
        eigen.J <- eigen(J)
        fit0 <- lmmbygls(null.formula, data=data, eigen.K=eigen.J, K=J, use.par=use.par, weights=weights, brute=brute)
        fit0.REML <- lmmbygls(null.formula, data=data, eigen.K=eigen.J, K=J, use.par="h2.REML", weights=weights, brute=brute)
      }
      else{
        eigen.K <- eigen(K)
        fit0 <- lmmbygls(null.formula, data=data, eigen.K=eigen.K, K=K, use.par=use.par, weights=weights, brute=brute)
        fit0.REML <- lmmbygls(null.formula, data=data, eigen.K=eigen.K, K=K, use.par="h2.REML", weights=weights, brute=brute)
      }
    }
    ## Null model for diplolasso
    if(model == "diplolasso"){
      if(fit0$h2 == 0 & is.null(weights)){
        this.y <- fit0$y
        this.x <- fit0$x
      }
      else{
        this.y <- fit0$M %*% fit0$y
        this.x <- fit0$M %*% fit0$x
      }
      fit0.glmnet <- glmnet::glmnet(y=this.y,
                                    x=cbind(this.x, h$getLocusMatrix(locus=h$getLoci()[1], model="additive")[1:length(this.y),]),
                                    nlambda=3, lambda.min=.9, family="gaussian",
                                    penalty.factor=c(rep(0, ncol(fit0$x) - 1), rep(1, choose(num.founders, 2))), intercept=TRUE,
                                    standardize=FALSE)
      diplolasso.penalty.factor <- c(rep(0, ncol(fit0$x) - 1 + num.founders - 1), rep(1, choose(num.founders, 2)))
    }
    ####### EMMA or EMMAX  
    if(use.fix.par){
      fix.par <- fit0$h2
    }
    if(!use.fix.par){
      fix.par <- NULL
    }
  }

  MI.LOD <- MI.p.value <- NULL
  LOD.vec <- p.vec <- df <- rep(NA, length(loci))
  null.data <- data
  
  ## Prepping imputation in multiple imputations
  if(use.multi.impute){
    impute.map <- data.frame(SUBJECT.NAME=data$SUBJECT.NAME, impute.on=data[,impute.on])
  }
  
  for(i in 1:length(loci)){
    if(use.multi.impute){
      if(i == 1){ # only at the beginning
        MI.LOD <- MI.p.value <- matrix(NA, nrow=num.imp, ncol=length(loci))
      }
      diplotype.prob.matrix <- h$getLocusMatrix(loci[i], model="full", subjects=old.data$SUBJECT.NAME)
      if(do.augment){
        if(model=="additive"){
          augment.matrix <- matrix(0, nrow=augment.n, ncol=choose(augment.n, 2) + augment.n)
          for(k in 1:augment.n){
            augment.matrix[k, k] <- 1
          }
        }
        if(model=="full"){
          augment.matrix <- diag(augment.n)
        }
        sample.names <- rownames(diplotype.prob.matrix)
        diplotype.prob.matrix <- rbind(diplotype.prob.matrix, augment.matrix)
        rownames(diplotype.prob.matrix) <- c(sample.names, paste0("augment.obs", 1:augment.n))
      }
      fit1 <- multi.imput.lmmbygls(num.imp=num.imp, data=data, formula=formula, founders=founders,
                                   diplotype.probs=diplotype.prob.matrix, 
                                   model=model, use.lmer=use.lmer, impute.map=impute.map,
                                   diplolasso.refit=diplolasso.refit, diplolasso.penalty.factor=diplolasso.penalty.factor,
                                   use.par=use.par, fix.par=fix.par, fit0=fit0, fit0.glmnet=fit0.glmnet, do.augment=do.augment, 
                                   brute=brute, seed=seed, weights=weights) 
      MI.LOD[,i] <- fit1$LOD
      MI.p.value[,i] <- fit1$p.value
      LOD.vec[i] <- median(fit1$LOD)
      p.vec[i] <- median(fit1$p.value)
    }
    if(!use.multi.impute){
      if(model %in% c("additive", "full")){
        X.check <- h$getLocusMatrix(loci[i], model=model)
        X <- h$getLocusMatrix(loci[i], model=model, subjects=data$SUBJECT.NAME[1:original.n])
        max.column <- which.max(colSums(X, na.rm=TRUE))[1]
        X <- X[,-max.column]
      }
      if(model == "diplolasso"){
        X.dosages <- h$getLocusMatrix(loci[i], model="additive", subjects=data$SUBJECT.NAME[1:original.n])
        max.column <- which.max(colSums(X.dosages))[1]
        X <- cbind(X.dosages[,-max.column],
                   h$getLocusMatrix(loci[i], model="full", subjects=data$SUBJECT.NAME[1:original.n])[,-(1:num.founders)])
      }
      colnames(X) <- gsub(pattern="/", replacement=".", x=colnames(X), fixed=TRUE)
      locus.formula <- make.alt.formula(formula=formula, X=X, do.augment=do.augment)
      if(do.augment){
        X.names <- rownames(X)
        if(model=="additive"){
          X <- rbind(X, 2*diag(augment.n)[,-max.column])
        }
        if(model=="full"){
          X <- rbind(X, diag(augment.n)[,-max.column])
        }
        rownames(X) <- c(X.names, paste0("augment.obs", 1:augment.n))
      }
      
      data <- cbind(null.data, X)
      if(use.lmer){
        fit1 <- lmmbylmer(formula=locus.formula, data=data, REML=FALSE, weights=weights)
        LOD.vec[i] <- log10(exp(as.numeric(logLik(fit1)) - as.numeric(logLik(fit0))))
        p.vec[i] <- pchisq(q=-2*(as.numeric(logLik(fit0)) - as.numeric(logLik(fit1))), df=length(fixef(fit1))-length(fixef(fit0)), lower.tail=FALSE)
      }
      else{
        if(model != "diplolasso"){
          fit1 <- lmmbygls(formula=locus.formula, data=data, 
                           eigen.K=fit0$eigen.K, K=fit0$K, 
                           use.par="h2", fix.par=fix.par, M=fit0$M, logDetV=fit0$logDetV,
                           brute=brute, 
                           weights=weights)
          LOD.vec[i] <- log10(exp(fit1$logLik - fit0$logLik))
          p.vec[i] <- pchisq(q=-2*(fit0$logLik - fit1$logLik), df=fit1$rank-fit0$rank, lower.tail=FALSE)
        }
        else{
          fit1 <- lmmbygls.diplolasso(locus.formula, data=data, 
                                      eigen.K=fit0$eigen.K, K=fit0$K,
                                      use.par="h2", fix.par=fix.par, M=fit0$M, logDetV=fit0$logDetV,
                                      brute=brute, 
                                      diplolasso.refit=diplolasso.refit, diplolasso.penalty.factor=diplolasso.penalty.factor,
                                      weights=weights, founders=founders)
          if(diplolasso.refit){
            LOD.vec[i] <- log10(exp(fit1$logLik - fit0$logLik))
            p.vec[i] <- pchisq(q=-2*(fit0$logLik - fit1$logLik), df=fit1$rank-fit0$rank, lower.tail=FALSE)
            df[i] <- fit1$rank
          }
          else{
            SSR.0 <- deviance(fit0.glmnet)[1] - 0.5*fit0$logDetV
            df.0 <- fit0.glmnet$df[1]
            LOD.vec[i] <- SSR.0 - fit1$SSR
            p.vec[i] <- pchisq(q=SSR.0 - fit1$SSR, df=fit1$rank - df.0, lower.tail=FALSE)
            df[i] <- fit1$rank
          }
        }
      }
      names(LOD.vec) <- names(p.vec) <- names(df) <- loci
    }
    if(print.locus.fit){
      cat("locus", i, "\n")
    }
  }
  output <- list(LOD=LOD.vec,
                 p.value=p.vec,
                 MI.LOD=MI.LOD,
                 MI.p.value=MI.p.value,
                 df=df,
                 pos=list(Mb=h$getMarkerLocation(loci, scale="Mb"), cM=h$getMarkerLocation(loci, scale="cM")),
                 loci=loci, 
                 chr=h$getChromOfLocus(loci),
                 fit0=fit0,
                 fit0.REML=fit0.REML,
                 y=data$y,
                 formula=formula.string,
                 model.type=model)
  return(output)
}

#' @export
single.locus.fit.h2lmm <- function(genomecache, data, formula, K, locus,
                                   model=c("additive", "full", "diplolasso"), diplolasso.refit=FALSE,
                                   use.par="h2", use.multi.impute=TRUE, num.imp=10, brute=TRUE, use.fix.par=FALSE, 
                                   seed=1, 
                                   weights=NULL, do.augment=FALSE, use.augment.weights=FALSE, use.full.null=FALSE, added.data.points=1, 
                                   ...) # diplolasso and diplolasso refit not work correctly
{
  model <- model[1]
  
  h <- bagpipe.backend::DiploprobReader$new(genomecache)
  founders <- h$getFounders()
  num.founders <- length(founders)
  cache.subjects <- rownames(h$getLocusMatrix(locus, model="additive"))
  
  data.and.K <- make.processed.data(formula=formula, data=data, cache.subjects=cache.subjects, K=K)
  data <- data.and.K$data
  K <- data.and.K$K
  if(!is.null(weights)){ weights <- weights[data$SUBJECT.NAME] }
  
  locus.chr <- h$getChromOfLocus(locus)
  
  augment.indicator <- NULL
  formula.string <- Reduce(paste, deparse(formula))
  null.formula <- make.null.formula(formula, do.augment=do.augment)
  original.n <- nrow(data)
  old.data <- data
  
  ###### Augmentation
  if(do.augment){
    augment.n <- ifelse(model=="additive", num.founders, num.founders + choose(num.founders, 2))
    augment.indicator <- c(rep(0, original.n), rep(1, augment.n))
    if(!use.full.null){
      data <- make.simple.augment.data(data=data, augment.n=augment.n)
      K <- make.simple.augment.K(K=K, augment.n=augment.n)
    }
    if(use.full.null){
      no.augment.K <- K
      K <- make.full.null.augment.K(K=no.augment.K, original.n=original.n, augment.n=augment.n)
      data <- make.full.null.augment.data(formula=formula, data=data, no.augment.K=no.augment.K, use.par=use.par, brute=brute,
                                          original.n=original.n, augment.n=augment.n, weights=weights)
    }
    if(use.augment.weights){
      weights <- make.augment.weights(data=data, augment.n=augment.n, added.data.points=added.data.points)
    }
  }
  
  ###### Null model
  ## No kinship effect - weights or no weights
  if(is.null(K)){
    fit0 <- lmmbygls(null.formula, data=data, eigen.K=NULL, K=NULL, use.par="h2", fix.par=0, weights=weights, brute=brute)
    fit0.REML <- lmmbygls(null.formula, data=data, eigen.K=NULL, K=NULL, use.par="h2.REML", fix.par=0, weights=weights, brute=brute)
  }
  ## Kinship effect - weights or no weights
  else{
    ###### This function is made for handling constant weights at all loci
    if(!is.null(weights)){
      J <- weights^(-1/2) * t(weights^(-1/2) * K)
      eigen.J <- eigen(J)
      fit0 <- lmmbygls(null.formula, data=data, eigen.K=eigen.J, K=J, use.par=use.par, weights=weights, brute=brute)
      fit0.REML <- lmmbygls(null.formula, data=data, eigen.K=eigen.J, K=J, use.par="h2.REML", weights=weights, brute=brute)
    }
    else{
      eigen.K <- eigen(K)
      fit0 <- lmmbygls(null.formula, data=data, eigen.K=eigen.K, K=K, use.par=use.par, weights=weights, brute=brute)
      fit0.REML <- lmmbygls(null.formula, data=data, eigen.K=eigen.K, K=K, use.par="h2.REML", weights=weights, brute=brute)
    }
  }
  ## Null model for diplolasso
  fit0.glmnet <- diplolasso.penalty.factor <- NULL
  if(model == "diplolasso"){
    if(fit0$h2 == 0 & is.null(weights)){
      this.y <- fit0$y
      this.x <- fit0$x
    }
    else{
      this.y <- fit0$M %*% fit0$y
      this.x <- fit0$M %*% fit0$x
    }
    fit0.glmnet <- glmnet::glmnet(y=this.y,
                                  x=cbind(this.x, h$getLocusMatrix(locus=h$getLoci()[1], model="additive")[1:length(this.y),]),
                                  nlambda=3, lambda.min=.9, family="gaussian",
                                  penalty.factor=c(rep(0, ncol(fit0$x) - 1), rep(1, choose(num.founders, 2))), intercept=TRUE,
                                  standardize=FALSE)
    diplolasso.penalty.factor <- c(rep(0, ncol(fit0$x) - 1 + num.founders - 1), rep(1, choose(num.founders, 2)))
  }
  
  ####### EMMA or EMMAX  
  if(use.fix.par){
    fix.par <- fit0$h2
  }
  if(!use.fix.par){
    fix.par <- NULL
  }
  
  null.data <- data
  
  if(use.multi.impute){
    diplotype.matrix <- h$getLocusMatrix(locus, model="full", subjects=old.data$SUBJECT.NAME)
    if(do.augment){
      if(model=="additive"){
        augment.matrix <- matrix(0, nrow=augment.n, ncol=choose(augment.n, 2) + augment.n)
        for(k in 1:augment.n){
          augment.matrix[k, k] <- 1
        }
      }
      if(model=="full"){
        augment.matrix <- diag(augment.n)
      }
      sample.names <- rownames(diplotype.matrix)
      diplotype.matrix <- rbind(diplotype.matrix, augment.matrix)
      rownames(diplotype.matrix) <- c(sample.names, paste0("augment.obs", 1:augment.n))
    }
    fit1 <- multi.imput.lmmbygls(num.imp=num.imp, data=data, formula=formula, founders=founders,
                                 diplotypes=diplotype.matrix, 
                                 model=model, 
                                 diplolasso.refit=diplolasso.refit, diplolasso.penalty.factor=diplolasso.penalty.factor,
                                 use.par=use.par, fix.par=fix.par, fit0=fit0, fit0.glmnet=fit0.glmnet, do.augment=do.augment, 
                                 brute=brute, seed=seed, weights=weights) 
      LOD <- median(fit1$LOD)
      p.val <- median(fit1$p.value)
    }
    else{
      if(model %in% c("additive", "full")){
        X.check <- h$getLocusMatrix(locus, model=model)
        X <- h$getLocusMatrix(locus, model=model, subjects=data$SUBJECT.NAME[1:original.n])
        max.column <- which.max(colSums(X, na.rm=TRUE))[1]
        X <- X[,-max.column]
      }
      if(model == "diplolasso"){
        X.dosages <- h$getLocusMatrix(locus, model="additive", subjects=data$SUBJECT.NAME[1:original.n])
        max.column <- which.max(colSums(X.dosages))[1]
        X <- cbind(X.dosages[,-max.column],
                   h$getLocusMatrix(locus, model="full", subjects=data$SUBJECT.NAME[1:original.n])[,-(1:num.founders)])
      }
      colnames(X) <- gsub(pattern="/", replacement=".", x=colnames(X), fixed=TRUE)
      locus.formula <- make.alt.formula(formula=formula, X=X, do.augment=do.augment)
      if(do.augment){
        X.names <- rownames(X)
        if(model=="additive"){
          X <- rbind(X, 2*diag(augment.n)[,-max.column])
        }
        if(model=="full"){
          X <- rbind(X, diag(augment.n)[,-max.column])
        }
        rownames(X) <- c(X.names, paste0("augment.obs", 1:augment.n))
      }
      
      data <- cbind(null.data, X)
      if(model != "diplolasso"){
        fit1 <- lmmbygls(locus.formula, data=data, 
                         eigen.K=fit0$eigen.K, K=fit0$K, 
                         use.par="h2", fix.par=fix.par, M=fit0$M, logDetV=fit0$logDetV,
                         brute=brute, 
                         weights=weights)
        LOD <- log10(exp(fit1$logLik - fit0$logLik))
        p.val <- pchisq(q=-2*(fit0$logLik - fit1$logLik), df=fit1$rank-fit0$rank, lower.tail=FALSE)
      }
      else{
        fit1 <- lmmbygls.diplolasso(locus.formula, data=data, 
                                    eigen.K=fit0$eigen.K, K=fit0$K,
                                    use.par="h2", fix.par=fix.par, M=fit0$M, logDetV=fit0$logDetV,
                                    brute=brute, 
                                    diplolasso.refit=diplolasso.refit, diplolasso.penalty.factor=diplolasso.penalty.factor,
                                    weights=weights, founders=founders)
        if(diplolasso.refit){
          LOD <- log10(exp(fit1$logLik - fit0$logLik))
          p.val <- pchisq(q=-2*(fit0$logLik - fit1$logLik), df=fit1$rank-fit0$rank, lower.tail=FALSE)
          df <- fit1$rank
        }
        else{
          SSR.0 <- deviance(fit0.glmnet)[1] - 0.5*fit0$logDetV
          df.0 <- fit0.glmnet$df[1]
          LOD <- SSR.0 - fit1$SSR
          p.val <- pchisq(q=SSR.0 - fit1$SSR, df=fit1$rank - df.0, lower.tail=FALSE)
          df <- fit1$rank
        }
    }
  }
  output <- list(LOD=LOD,
                 p.value=p.val,
                 df=df,
                 pos=list(Mb=h$getMarkerLocation(locus, scale="Mb"), cM=h$getMarkerLocation(locus, scale="cM")),
                 locus=locus, 
                 chr=h$getChromOfLocus(locus),
                 fit0=fit0,
                 fit0.REML=fit0.REML,
                 fit1=fit1,
                 formula=formula.string,
                 model.type=model)
  return(output)
}
