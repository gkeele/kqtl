
make.null.formula <- function(formula, do.augment){
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula.string <- paste0("y ~ ", unlist(strsplit(this.formula.string, split="~"))[-1])
  this.formula <- as.formula(ifelse(do.augment, paste0(this.formula.string, " + augment.indicator"), this.formula.string))
  return(this.formula)
}
make.alt.formula <- function(formula, X, do.augment){
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula.string <- paste0("y ~ ", unlist(strsplit(this.formula.string, split="~"))[-1])
  this.formula.string <- ifelse(do.augment, paste0(this.formula.string, " + augment.indicator"), this.formula.string)
  this.formula <- as.formula(paste(this.formula.string, paste(gsub(pattern="/", replacement=".", x=colnames(X), fixed=TRUE), collapse=" + "), sep=" + "))
  return(this.formula)
}
make.snp.alt.formula <- function(formula){
  this.formula.string <- Reduce(paste, deparse(formula))
  this.formula.string <- paste0("y ~ ", unlist(strsplit(this.formula.string, split="~"))[-1])
  this.formula <- as.formula(paste(this.formula.string, "SNP", sep=" + "))
  return(this.formula)
}
remove.whitespace.formula <- function(formula){
  formula.string <- paste0(Reduce(paste, deparse(formula)))
  formula.string <- gsub("[[:space:]]", "", formula.string)
  return(as.formula(formula.string))
}
check.for.lmer.formula <- function(formula){
  formula.string <- paste0(Reduce(paste, deparse(formula)))
  formula.string <- gsub("[[:space:]]", "", formula.string)
  use.lmer <- grepl(pattern="\\([a-zA-Z0-9\\.]+\\|[a-zA-Z0-9\\.]+\\)", x=formula.string, perl=TRUE)
  return(use.lmer)
}
make.processed.data <- function(formula, data, cache.subjects, K, impute.on){
  all.variables <- all.vars(formula)
  covariates <- all.variables[-1]
  lh.formula.string <- unlist(strsplit(Reduce(paste, deparse(formula)), split="~"))[1]
  lh.formula.string <- gsub("[[:space:]]", "", lh.formula.string)
  covariates <- c(covariates, unique(c("SUBJECT.NAME", impute.on)))
  formula.string <- paste(lh.formula.string,
                          paste(covariates, collapse="+"),
                          sep="~")
  data <- model.frame(formula(formula.string), data=data)
  names(data) <- c("y", covariates)
  # selecting those in both data and cache
  #data <- data[as.character(data$SUBJECT.NAME) %in% cache.subjects,]
  cache.subjes <- cache.subjects[cache.subjects %in% as.character(data$SUBJECT.NAME)]
  data <- data[match(x=as.character(data$SUBJECT.NAME), table=cache.subjects, nomatch=0),]
  if(!is.null(K)){
    #data <- data[as.character(data$SUBJECT.NAME) %in% colnames(K),]
    K <- K[as.character(data$SUBJECT.NAME), as.character(data$SUBJECT.NAME)]
  }
  if(length(covariates) > 0){
    covariate.matrix <- matrix(NA, nrow=nrow(data), ncol=length(covariates))
    for(i in 1:length(covariates)){
      if(is.factor(data[,covariates[i]])){
        factor.counts <- table(data[,covariates[i]])
        data[,covariates[i]] <- gdata::reorder.factor(x=data[,covariates[i]], new.order=names(sort(factor.counts[factor.counts != 0], decreasing=TRUE)))
      }
    }
  }
  return(list(data=data, K=K))
}

make.simple.augment.K <- function(K, augment.n){
  if(!is.null(K)){
    original.K.names <- colnames(K)
    K <- as.matrix(bdiag(K, diag(augment.n)))
    rownames(K) <- colnames(K) <- c(original.K.names, paste0("augment.obs", 1:augment.n))
  }
  return(K)
}
make.simple.augment.data <- function(data, formula, augment.n){
  real.y.names <- data$SUBJECT.NAME
  all.variables <- all.vars(formula)
  covariates <- all.variables[-1]
  augment.y <- rep(mean(data$y), augment.n)
  augment.y.names <- paste0("augment.obs", 1:augment.n)
  covariate.matrix <- NULL
  if(length(all.variables) > 1){
    covariate.matrix <- matrix(NA, nrow=augment.n, ncol=length(covariates))
    for(i in 1:length(covariates)){
      if(is.factor(data[,covariates[i]])){
        covariate.matrix[,i] <- rep(levels(data[,covariates[i]])[1], augment.n)
      }
      if(is.numeric(data[,covariates[i]])){
        covariate.matrix[,i] <- rep(mean(data[,covariates[i]]), augment.n)
      }
    }
  }
  if(is.null(covariate.matrix)){
    partial.augment.data <- data.frame(augment.y, augment.y.names)
    names(partial.augment.data) <- c("y", "SUBJECT.NAME")
  }
  if(!is.null(covariate.matrix)){
    partial.augment.data <- data.frame(augment.y, covariate.matrix, augment.y.names)
    names(partial.augment.data) <- c("y", covariates, "SUBJECT.NAME")
  }
  data <- rbind(data, partial.augment.data)
  return(data)
}
make.augment.weights <- function(data, augment.n, added.data.points){
  weights <- c(rep(1, nrow(data) - augment.n), 
               rep(augment.n/added.data.points, augment.n))
  return(weights)
}

make.full.null.augment.K <- function(K, augment.n, original.n){
  if(!is.null(K)){
    original.K.names <- colnames(K)
    K <- as.matrix(bdiag(K, diag(augment.n)))
    K[-(1:original.n), 1:original.n] <- K[1:original.n, -(1:original.n)] <- 0.5
    K[-(1:original.n), -(1:original.n)][K[-(1:original.n), -(1:original.n)] == 0] <- 0.5
    rownames(K) <- colnames(K) <- c(original.K.names, paste0("augment.obs", 1:augment.n))
  }
  return(K)
}

make.full.null.augment.data <- function(formula, data, no.augment.K, use.par, brute,
                                        original.n, augment.n, weights){
  all.variables <- all.vars(formula)
  covariates <- all.variables[-1]
  
  null.formula.no.augment <- make.null.formula(formula=formula, is.augmented=FALSE)
  fit0.no.augment <- lmmbygls(formula=null.formula.no.augment, data=data, covariates=covariates, K=no.augment.K,
                              use.par=use.par, brute=brute, null.test=TRUE)
  set.seed(seed)
  y.null.hat <- predict.lmmbygls(fit0.no.augment=fit0.no.augment, original.n=original.n, augment.n=augment.n, 
                                 covariates=covariates, weights=weights)
  real.y.names <- data$SUBJECT.NAME
  augment.y.names <- paste0("augment.obs", 1:augment.n)
  
  covariate.matrix <- NULL
  if(!is.null(covariates)){
    covariate.matrix <- matrix(NA, nrow=augment.n, ncol=length(covariates))
    for(i in 1:length(covariates)){
      if(is.factor(data[,covariates[i]])){
        covariate.matrix[,i] <- rep(levels(data[,covariates[i]])[1], augment.n)
      }
      if(is.numeric(data[,covariates[i]])){
        covariate.matrix[,i] <- rep(mean(data[,covariates[i]]), augment.n)
      }
    }
  }
  partial.augment.data <- data.frame(y.null.hat, covariate.matrix, augment.y.names)
  names(partial.augment.data) <- c(outcome, covariates, "SUBJECT.NAME")
  data <- rbind(data, partial.augment.data)
  data <- cbind(data, data.frame(augment.indicator=c(rep(0, original.n), rep(1, augment.n))))
  return(data)
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
  diplotype.probs <- data.frame(original.order=1:nrow(diplotype.probs), SUBJECT.NAME=rownames(diplotype.probs), diplotype.probs, stringsAsFactors=FALSE)
  diplotype.probs <- merge(x=diplotype.probs, y=impute.map, by="SUBJECT.NAME")
  diplotype.probs <- diplotype.probs[order(diplotype.probs$original.order),]
  diplotype.probs <- diplotype.probs[, names(diplotype.probs) != "original.order"]
  
  imputable.diplotype.probs <- diplotype.probs[as.integer(rownames(unique(data.frame(diplotype.probs[,"impute.on"])))),]
  rownames(imputable.diplotype.probs) <- imputable.diplotype.probs[, "impute.on"]
  imputable.diplotype.probs <- imputable.diplotype.probs[,!(names(imputable.diplotype.probs) %in% names(impute.map))]
  imputation <- t(apply(imputable.diplotype.probs, 1, function(x) rmultinom(1, 1, x)))
  full.imputation <- imputation[as.character(impute.map[, "impute.on"]),]
  rownames(full.imputation) <- impute.map[, "SUBJECT.NAME"]
  return(full.imputation)
}
  

