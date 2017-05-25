multi.imput.lmmbygls <- function(num.imp, data, formula,
                                 founders=founders, diplotype.probs, K=NULL, fit0=NULL, fit0.glmnet=NULL,
                                 use.par, fix.par=NULL, model=c("additive", "full"),
                                 use.lmer, impute.map,
                                 brute=TRUE, seed=1, do.augment,
                                 weights=NULL){
  model <- model[1]  
  eigen.K <- logDetV <- M <- NULL
  if(is.null(fit0)){
    null.formula <- make.null.formula(formula=formula, do.augment=do.augment)
    if(use.lmer){
      fit0 <- lmmbylmer(null.formula, data=data, REML=FALSE, weights=weights)
    }
    else{
      fit0 <- lmmbygls(null.formula, data=data, K=K, use.par=use.par, brute=brute, weights=weights)
      K <- fit0$K
    }
  }
  if(is.null(weights) & !use.lmer){
    eigen.K <- fit0$eigen.K
  }
  if(!is.null(fix.par) & !use.lmer){
    M <- fit0$M
    logDetV <- fit0$logDetV
  }
  full.to.dosages <- straineff.mapping.matrix()
  
  imp.logLik <- imp.h2 <- imp.df <- rep(0, num.imp)

  null.data <- data
  set.seed(seed)
  for(i in 1:num.imp){
    if(model == "additive"){
      if(any(diplotype.probs < 0)){
        diplotype.probs[diplotype.probs < 0] <- 0
        diplotype.probs <- t(apply(diplotype.probs, 1, function(x) x/sum(x)))
      }
      X <- run.imputation(diplotype.probs=diplotype.probs, impute.map=impute.map) %*% full.to.dosages
      max.column <- which.max(colSums(X))[1]
      X <- X[,-max.column]
      colnames(X) <- gsub(pattern="/", replacement=".", x=founders, fixed=TRUE)[-max.column]
    }
    if(model == "full"){
      X <- run.imputation(diplotype.probs=diplotype.probs, impute.map=impute.map)
      max.column <- which.max(colSums(X))[1]
      X <- X[,-max.column]
      colnames(X) <- gsub(pattern="/", replacement=".", x=colnames(diplotype.probs), fixed=TRUE)[-max.column]
    }

    locus.formula <- make.alt.formula(formula=formula, X=X, do.augment=do.augment)
    data <- cbind(null.data, X)
    if(use.lmer){
      fit1 <- lmmbylmer(formula=locus.formula, data=data, REML=FALSE, weights=weights)
      imp.logLik[i] <- as.numeric(logLik(fit1))
      imp.h2[i] <- NA
      imp.df[i] <- length(fixef(fit1))
    }
    else{
      fit1 <- lmmbygls(locus.formula, data=data, eigen.K=eigen.K, K=K,
                       logDetV=logDetV, M=M, 
                       use.par="h2", fix.par=fix.par,
                       brute=brute, weights=weights)
      imp.logLik[i] <- fit1$logLik
      imp.h2[i] <- fit1$h2
      imp.df[i] <- fit1$rank
    }
  }
  ## Summarizing over imputations
  if(use.lmer){
    imp.LOD <- log10(exp(imp.logLik - as.numeric(logLik(fit0))))
    imp.p.value <- pchisq(q=-2*(as.numeric(logLik(fit0)) - imp.logLik), df=imp.df - length(fixef(fit0)), lower.tail=FALSE)
  }
  else{
    imp.LOD <- log10(exp(imp.logLik - fit0$logLik))
    imp.p.value <- pchisq(q=-2*(fit0$logLik - imp.logLik), df=imp.df - fit0$rank, lower.tail=FALSE)
  }
  return(list(h2=imp.h2, 
              LOD=imp.LOD,
              p.value=imp.p.value))
}

multi.imput.lmmbygls.random <- function(num.imp, data, this.formula,
                                 diplotypes, K=NULL, fit0,
                                 use.par, model=c("additive", "full"),
                                 seed=1, do.augment){
  
  model <- model[1]
  K <- fit0$K
  eigen.K <- fit0$eigen.K

  full.to.dosages <- straineff.mapping.matrix()
  
  imp.logLik <- imp.h2 <- imp.p.value <- rep(0, num.imp)

  null.data <- data
  for(i in 1:num.imp) {
    set.seed(seed)
    if(model=="additive"){
      X <- t(apply(diplotypes, 1, function(x) rmultinom(1, 1, x))) %*% full.to.dosages
      colnames(X) <- h$getFounders()
    }
    if(model=="full"){
      X <- t(apply(diplotypes, 1, function(x) rmultinom(1, 1, x)))
      colnames(X) <- colnames(diplotypes)
    }
    fit1 <- random.lmmbygls(formula=this.formula, data=data,
                            use.par=use.par, fit0=fit0, Z=X,
                            null.h2=fit0$h2)
    
    chi.sq <- -2*(fit0$REML.logLik - fit1$REML.logLik)
    imp.p.value[i] <- ifelse(chi.sq == 0, 1, 0.5*pchisq(q=chi.sq, df=1, lower.tail=FALSE))
    imp.logLik[i] <- fit1$REML.logLik
    imp.h2[i] <- fit1$h2
  }

  return(list(logLik=imp.logLik, 
              h2=imp.h2, 
              p.value=imp.p.value))
}
