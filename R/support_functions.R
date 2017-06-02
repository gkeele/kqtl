#' Returns the rank-based inverse normal transformation
#'
#' This function takes a phenotype vector and returns the rank-based inverse normal transformation.
#'
#' @param phenotype A vector of phenotype values for which the rank-based inverse normal transformation is output.
#' @param prop DEFAULT: 0.5. This allows Inf to not be returned for the maximum of phenotype.
#' @export
#' @examples rint()
rint <- function(phenotype, prop=0.5){
  rint_phenotype <- qnorm((rank(phenotype)-prop)/length(phenotype))
  return(rint_phenotype)
}

predict.lmmbygls <- function(fit0.no.augment, original.n, augment.n, covariates, weights){
  e <- rnorm(augment.n, 0, sd=sqrt(fit0.no.augment$sigma2.mle))
  if(!is.null(weights)){
    e <- sqrt(weights[-(1:original.n)]) * e
  }
  u <- rmvnorm(n=1, mean=rep(0, original.n + augment.n), sigma=fit0.no.augment$tau2.mle*K, method="chol")[-(1:original.n)]
  covariate.matrix <- rep(1, augment.n)
  if(!is.null(covariates)){
    for(i in 1:length(covariates)){
      if(is.factor(data[,covariates[i]])){
        covariate.matrix <- cbind(covariate.matrix, matrix(0, nrow=augment.n, ncol=nlevels(data[,covariates[i]])-1))
      }
      if(is.numeric(data[,covariates[i]])){
        covariate.matrix <- cbind(covariate.matrix, rep(mean(data[,covariates[i]]), augment.n))
      }
    }
  }
  null.mean <- (rbind(fit0.no.augment$x, covariate.matrix) %*% fit0.no.augment$coefficients)[-(1:original.n),]
  null.y.hat <- null.mean + u + e
  return(null.y.hat)
}

#### Not in use currently #####
calc.LRT.mean.coef <- function(){
  mean.coef <- colMeans(imp.coef, na.rm=TRUE)
  mean.varcomps <- colMeans(imp.varcomps)
  lambda <- sum(mean.varcomps)
  sample.size <- nrow(imp.X[[1]])
  H.inv <- t(fit1$M)%*%fit1$M
  imp.constr.logLik <- rep(0, length(imp.X))
  y <- fit1$y
  if(is.null(weights)){
    for(i in 1:length(imp.X)){
      #this.X <- imp.X[[i]][,-ncol(imp.X[[1]])]
      this.X <- imp.X[[i]]
      imp.constr.logLik[i] <- -0.5*sample.size*(log(2*pi) + log(lambda)) - 0.5*(1/lambda)*t(y - this.X[, !is.nan(mean.coef)]%*%mean.coef[!is.nan(mean.coef)]) %*% H.inv %*% (y - this.X[, !is.nan(mean.coef)]%*%mean.coef[!is.nan(mean.coef)]) - 0.5*fit1$logDetV
    }
  }
  if(!is.null(weights)){
    for(i in 1:length(X.list)){
      this.X <- imp.X[[i]][,-ncol(imp.X[[i]])]
      imp.constr.logLik[i] <- -0.5*sample.size*log(2*pi) - 0.5*t(y - this.X[, !is.nan(mean.coef)]%*%mean.coef[!is.nan(mean.coef)]) %*% H.inv %*% (y - this.X[, !is.nan(mean.coef)]%*%mean.coef[!is.nan(mean.coef)]) - 0.5*fit1$logDetV
    }
  }
  imp.constr.LRT <- 2*(imp.constr.logLik - fit0$logLik)
  return(imp.constr.LRT)
}

calc.mi.LRT <- function(){
  num.imp <- length(imp.constr.LRT)
  k <- df1 - df0
  
  rL <- (num.imp + 1)*(mean(imp.LRT) - mean(imp.constr.LRT))/(k*(num.imp - 1))
  
  DL <- mean(imp.constr.LRT)/k*(1 + rL)
  
  v <- k*(num.imp - 1)
  w.rL <- ifelse(v > 4, 4 + (v-4)*(1 + (1 - 2/v)*(1/rL))^2, (1/2)*v*(1 + 1/k)*(1 + (1/rL))^2)
  p.val <- pf(DL, df1=k, df2=w.rL, lower.tail=FALSE)
  return(p.val)
}
###########################

#' @export
straineff.mapping.matrix <- function(M=8){
  T <- M*(M+1)/2
  mapping <- matrix(rep(0, T*M), M, T)
  idx <- 1;
  for (i in 1:M){
    mapping[i, idx] <- mapping[i, idx] + 2
    idx <- idx + 1;
  }
  for (i in 2:M){
    for (j in 1:(i-1)){
      mapping[i, idx] <- mapping[i, idx] + 1;
      mapping[j, idx] <- mapping[j, idx] + 1;
      idx <- idx + 1;
    }
  }
  return(t(mapping))
}

run.imputation <- function(diplotype.probs, impute.map){
  pheno.id <- names(impute.map)[1]
  geno.id <- names(impute.map)[2]
  diplotype.probs <- data.frame(1:nrow(diplotype.probs), rownames(diplotype.probs), diplotype.probs, stringsAsFactors=FALSE)
  names(diplotype.probs)[1:2] <- c("original.order", geno.id)
  diplotype.probs <- merge(x=diplotype.probs, y=impute.map, by=geno.id)
  diplotype.probs <- diplotype.probs[order(diplotype.probs$original.order),]
  diplotype.probs <- diplotype.probs[, names(diplotype.probs) != "original.order"]
  rownames(diplotype.probs) <- diplotype.probs[,geno.id]
  
  imputable.diplotype.probs <- diplotype.probs[(unique(as.character(rownames(diplotype.probs)))),]
  imputable.diplotype.probs <- imputable.diplotype.probs[,!(names(imputable.diplotype.probs) %in% names(impute.map))]
  imputation <- t(apply(imputable.diplotype.probs, 1, function(x) rmultinom(1, 1, x)))
  full.imputation <- imputation[as.character(impute.map[, geno.id]),]
  rownames(full.imputation) <- impute.map[, pheno.id]
  return(full.imputation)
}
  
process_eigen_decomposition <- function(eigen.decomp, tol=1e-6){
  # from Robert, who took it from MASS::mvrnorm()
  if(!all(eigen.decomp$values >= -tol * abs(eigen.decomp$values[1L]))){
    stop("K is not positive definite")
  }
  if(any(eigen.decomp$values < 0)){
    if(any(eigen.decomp$values < -tol)){
      message("Zeroing negative eigenvalues: smallest eigenvalue was ", min(eigen.decomp$values), "\n")
    }
    eigen.decomp$values <- pmax(eigen.decomp$values, 0)
  }
  return(eigen.decomp)
}

replicates.eigen <- function(Z, K) {
  eigen <- eigen(K %*% crossprod(Z,Z ), symmetric=FALSE)
  return(list(values=eigen$values,
              vectors=qr.Q(qr(Z %*% eigen$vectors))))
}

get.p.value <- function(fit0, fit1, method=c("LRT", "ANOVA")){
  method <- method[1]
  if(method == "LRT"){
    p.value <- pchisq(q=-2*(fit0$logLik - fit1$logLik), df=fit1$rank-fit0$rank, lower.tail=FALSE)
  }
  if(method == "ANOVA"){
    class(fit0) <- class(fit1) <- "lm"
    p.value <- anova(fit0, fit1)$`Pr(>F)`[2]
  }
  return(p.value)
}




