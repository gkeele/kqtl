gls.fit.diplolasso <- function(X, y,  M, logDetV, rotate=TRUE, diplolasso.refit=FALSE, penalty.factor, nlambda=500, lambda.min.ratio=0.01, ...){
  n  <- nrow(X)
  # Tranform to uncorrelated Sigma form
  if(!is.null(M)){
    MX <- M %*% X
    My <- M %*% y
  }
  if(is.null(M)){ # more efficient than MX <- diag(nrow(X)) %*% X
    MX <- X
    My <- y
  }
  
  ## Actual LASSO fit
  fit <- glmnet::glmnet(y=My, x=MX, intercept=TRUE, nlambda=nlambda, lambda.min.ratio=lambda.min.ratio, family="gaussian", penalty.factor=penalty.factor, standardize=F)
  deviance.fit <- deviance(fit)
  n <- length(y)
  bic <- n*log(deviance.fit/n) + log(n)*fit$df
  #bic <- deviance.fit + log(n)*fit$df # technically BIC + constant (constant = 2*log(saturated likelihood))
  my.model <- which.min(bic)
  include.betas <- which(fit$beta[, my.model] != 0)
  
  if(diplolasso.refit){
    MX <- MX[, c(1, include.betas), drop=FALSE]
    fit <- lm.fit(x=MX, y=My, ...)
    fit$rss <- sum(fit$residuals^2)
    fit$sigma2.mle <- fit$rss/n
    fit$uncorr.logLik <- -0.5*n*log(2*pi) - 0.5*n*log(fit$sigma2.mle) - 0.5*n
    fit$logLik <- fit$uncorr.logLik - 0.5*logDetV
    fit$logDetV <- logDetV
    fit$M <- M
  }
  else{
    fit$SSR <- deviance.fit[which.min(bic)] - 0.5*logDetV # Effectively likelihood, with the saturated model as a constant
    fit$rank <- fit$df[which.min(bic)]
  }
  # Correct logLik to value under correlated Sigma
  return(fit)
}

lmmbygls.diplolasso <- function(formula, data, K=NULL, eigen.K=NULL, fit0, fix.par=NULL,
                                M=NULL, logDetV=NULL, weights,
                                use.par=c("h2", "h2.REML"), 
                                brute=TRUE, 
                                diplolasso.refit=FALSE, diplolasso.penalty.factor,
                                founders,
                                subset, na.action,
                                method = "qr",
                                model = TRUE, 
                                contrasts = NULL,
                                verbose = FALSE,
                                ...) 
{
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$W <- m$K <- m$eigen.K <- m$fix.par <- m$use.par <- NULL
  m$M <- m$logDetV <- m$weights <- NULL
  m$method <- m$model <- m$x <- m$y <- m$contrasts <- m$verbose <- NULL
  m$brute <- m$diplolasso.refit <- m$diplolasso.penalty.factor <- m$founders <- m$... <- NULL
  
  m[[1L]] <- quote(stats::model.frame)
  m <- eval.parent(m)
  if (method == "model.frame"){
    return(m)    
  }
  Terms <- attr(m, "terms")
  y <- model.response(m)
  X <- model.matrix(Terms, m, contrasts)
  n <- nrow(X)
  q <- ncol(X)
  remove.for.diplolasso <- which(colnames(X) %in% c("(Intercept)", founders[1]))

  if(is.null(fix.par)){
    if(is.null(eigen.K)){
      if(is.null(K)){ ## No kinship effect setting: K - NULL, eigen.K - NULL
        fix.par <- 0
      }
      else{
        eigen.K <- eigen(K)
      }
    }
    Ut <- t(eigen.K$vectors) # a small optimization
  }
  ## Define local objective function/closure for Brent's optimization
  ### Optimize functions
  h2.fit <- function(h2, logLik.only=TRUE, verbose=FALSE, ...){
    if(is.null(fix.par)){ # EMMA or first time optimization
      if(is.null(weights)){
        d <- h2*eigen.K$values + (1-h2)
        M <- d^-0.5 * Ut 
        logDetV <- sum(log(d))
      }
      else{
        d <- h2*eigen.K$values + (1-h2)
        M <- d^-0.5 * t((1/sqrt(weights)) * t(Ut))
        logDetV <- sum(log(1/weights)) + sum(log(d)) # maybe right
      }
    }
    else{
      if(fix.par == 0 & is.null(weights)){
        M <- NULL # more efficient than I
        logDetV <- 0
      }
      if(fix.par == 0 & !is.null(weights)){
        M <- diag(1/sqrt(weights))
        logDetV <- sum(log(1/weights))
      }
    }
    fit <- gls.fit.diplolasso(X=X, y=y, M=M, logDetV=logDetV, diplolasso.refit=diplolasso.refit, penalty.factor=diplolasso.penalty.factor, ...)
    if(logLik.only){
      if(verbose){
        cat(sep="", "h2 = ", h2, " : logLik = ", fit$logLik, "\n")
      }
      return(fit$logLik)
    }
    fit$h2 <- h2
    return(fit)
  }
  h2.fit.REML <- function(h2, logLik.only=TRUE, verbose=FALSE, ...){
    if(is.null(fix.par)){
      if(is.null(weights)){
        d <- h2*eigen.K$values + (1-h2)
        M <- d^-0.5 * Ut 
        logDetV <- sum(log(d))
      }
      else{
        d <- h2*eigen.K$values + (1-h2)
        M <- d^-0.5 * t(sqrt(weights) * Ut)
        logDetV <- sum(log(weights)) + sum(log(d)) # maybe right
      }
    }
    else{
      if(fix.par == 0 & is.null(weights)){
        M <- NULL # more efficient than I
        logDetV <- 0
      }
      if(fix.par == 0 & !is.null(weights)){
        M <- diag(sqrt(weights))
        logDetV <- sum(log(weights))
      }
    }
    fit <- gls.fit.diplolasso(X=X, y=y, M=M, logDetV=logDetV, diplolasso.refit=diplolasso.refit, penalty.factor=diplolasso.penalty.factor, ...)
    adjusted.logLik <- -0.5*n*log(2*pi) - 0.5*n*log(fit$sigma2) - 0.5*(n-ncol(X)) - 0.5*logDetV 
    REML.logLik <- adjusted.logLik + 0.5*(ncol(X)*log(2*pi*fit$sigma2) + log(det(t(X)%*%X)) - log(det(t(X)%*%t(Ut)%*%diag(1/d)%*%Ut%*%X)))
    fit$REML.logLik <- REML.logLik
    if (logLik.only){
      if (verbose){
        cat(sep="", "h2 = ", h2, " : logLik = ", fit$REML.logLik, "\n")
      }
      return (fit$REML.logLik)
    }
    fit$h2 <- h2
    return(fit)
  }
  ## Optimize using objective function, or otherwise given value of h2
  fit <- NULL
  if(is.null(fix.par)){
    if(verbose){
      cat("Optimizing logLik for parameter\n")
    }
    if(use.par[1] == "h2"){
      peak <- optimize(f=h2.fit, logLik.only=TRUE, verbose=verbose, ..., interval=c(0,1), maximum=TRUE)
      fit  <- h2.fit(h2=peak$maximum, logLik.only=FALSE, verbose=FALSE)
      if(brute){
        fit.h2.0 <- h2.fit(h2=0, logLik.only=FALSE, verbose=FALSE)
        if(fit$logLik < fit.h2.0$logLik){
          fit <- fit.h2.0
        }
        fit.h2.1 <- h2.fit(h2=1, logLik.only=FALSE, verbose=FALSE)
        if(fit$logLik < fit.h2.1$logLik){
          fit <- fit.h2.1
        }
      }
    }
    if(use.par[1] == "h2.REML"){
      peak <- optimize(f=h2.fit.REML, ..., interval=c(0,1), maximum=TRUE)
      fit  <- h2.fit.REML(h2=peak$maximum, logLik.only=FALSE, verbose=FALSE)
    }
    fit$h2.optimized <- TRUE
  } 
  if(!is.null(fix.par)) {
    if(use.par[1] == "h2"){
      fit <- h2.fit(h2=fix.par, logLik.only=FALSE, verbose=FALSE)
    }
    fit$h2.optimized <- FALSE
  }
  fit$gls.sigma2.mle <- fit$sigma2.mle
  fit$gls.sigma2     <- fit$sigma2
  fit$sigma2.mle <- (1 - fit$h2)*fit$gls.sigma2.mle
  fit$tau2.mle <- fit$h2*fit$gls.sigma2.mle
  
  fit$terms <- Terms
  fit$call <- call
  if (model){
    fit$model <- m
  }
  fit$na.action <- attr(m, "na.action")
  fit$x <- X
  fit$y <- y
  fit$eigen.K <- eigen.K
  fit$K <- K
  fit$xlevels <- .getXlevels(Terms, m)
  fit$contrasts <- attr(X, "contrasts")
  class(fit) <- "lmmbygls"
  return(fit)
}


### Old version, glmnet has issues when fitting a no intercept model 
###   - which unfortunately makes us use a bit of a hack
# gls.fit.diplolasso <- function(X, y,  M, logDetV, rotate=TRUE, diplolasso.refit=FALSE, penalty.factor, nlambda=500, lambda.min.ratio=0.01, ...){
#   n  <- nrow(X)
#   # Tranform to uncorrelated Sigma form
#   if(!is.null(M)){
#     MX <- M %*% X
#     My <- M %*% y
#   }
#   if(is.null(M)){ # more efficient than MX <- diag(nrow(X)) %*% X
#     MX <- X
#     My <- y
#   }
#   
#   ## Annoying calculation because glmnet seems to only produce deviances - find the log likelihood of the saturated model (LL_sat)
#   if(is.null(M)){
#     fit.null.glmnet <- glmnet(y=My, x=MX[,-1], intercept=TRUE, nlambda=3, lambda.min.ratio=0.9, family="gaussian", penalty.factor=rep(1, ncol(MX) - 1), standardize=F)
#   }
#   else{
#     fit.null.glmnet <- glmnet(y=My, x=MX, intercept=FALSE, nlambda=3, lambda.min.ratio=0.9, family="gaussian", penalty.factor=c(0, rep(1, ncol(MX) - 1)), standardize=F)
#   }
#   deviance.null <- deviance(fit.null.glmnet)[1]
#   this.fit0 <- lm.fit(y=My, x=MX[,1, drop=FALSE])
#   this.fit0$rss <- sum(this.fit0$residuals^2)
#   this.fit0$sigma2.mle <- this.fit0$rss/n
#   this.fit0$uncorr.logLik <- -0.5*n*log(2*pi) - 0.5*n*log(this.fit0$sigma2.mle) - 0.5*n
#   LL_sat <- deviance.null/2 + this.fit0$uncorr.logLik
#   
#   ## Actual LASSO fit
#   if(is.null(M)){
#     fit <- glmnet(y=My, x=MX[,-1], intercept=TRUE, nlambda=nlambda, lambda.min.ratio=lambda.min.ratio, family="gaussian", penalty.factor=penalty.factor[-1], standardize=F)
#   }
#   else{
#     fit <- glmnet(y=My, x=MX, intercept=FALSE, nlambda=nlambda, lambda.min.ratio=lambda.min.ratio, family="gaussian", penalty.factor=penalty.factor, standardize=F)
#   }
#   deviance.fit <- deviance(fit)
#   
#   n <- length(My)
#   LL <- LL_sat - deviance.fit/2
#   bic <- LL + log(n)*fit$df
#   my.model <- which.min(bic)
#   include.betas <- fit$beta[, my.model] != 0
#   #SSR <- deviance(fit)
#   #n <- length(y)
#   #bic <- n*log(SSR/n) + log(n)*fit$df
#   #my.model <- which.min(bic)
#   #include.betas <- fit$beta[, my.model] != 0
#   #browser()
#   if(diplolasso.refit){
#     if(is.null(M)){
#       
#     }
#     else{
#       
#     }
#     MX <- M %*% X[, include.betas, drop=FALSE]
#     fit <- lm.fit(x=MX, y=My, ...)
#     fit$rss <- sum(fit$residuals^2)
#     fit$sigma2 <- fit$rss/fit$df.residual
#   }
#   else{
#     #SSr <- SSR[which.min(bic)]
#     #df.bic <- fit$df[which.min(bic)]
#     if(is.null(M)){
#       fit$rss <- sum((My - cbind(MX[,1, drop=FALSE], MX[,-1][,include.betas, drop=FALSE]) %*% c(fit$a0[my.model], fit$beta[,my.model][include.betas]))^2)
#       fit$rank <- fit$df[my.model] + 1
#     }
#     else{
#       fit$rss <- sum((My - MX[, include.betas, drop=FALSE] %*% fit$beta[,my.model][include.betas])^2)
#       fit$rank <- fit$df[my.model]
#     }
#     fit$sigma2 <- fit$rss/fit$rank
#   }
#   fit$sigma2.mle <- fit$rss/n
#   fit$uncorr.logLik <- -0.5*n*log(2*pi) - 0.5*n*log(fit$sigma2.mle) - 0.5*n
#   # Correct logLik to value under correlated Sigma
#   fit$logLik <- fit$uncorr.logLik - 0.5*logDetV
#   fit$logDetV <- logDetV
#   fit$M <- M
#   return(fit)
# }


