#' Returns a matrix of outcome samples, either permutations or from the null model of no locus effect
#'
#' This function takes an scan.h2lmm() object, and returns a specified number of outcome samples, either permutations or
#' from the null model of no locus effect.
#'
#' @param scan.object A scan.h2lmm() object.
#' @param method DEFAULT: "bootstrap". "bootstrap" specifies parametric bootstraps from the null model. "permutation" specifies
#' permutations that can respect the structure of the data. Permutations are more appropriate if the data have highly
#' influential data points.
#' @param use.REML DEFAULT: TRUE. Determines whether the variance components for the parametric sampling are 
#' based on maximizing the likelihood (ML) or the residual likelihood (REML).
#' @param num.samples The number of parametric bootstrap samples to return.
#' @param seed DEFAULT: 1. The sampling process is random, thus a seed must be set for samples to be consistent
#' across machines.
#' @export
#' @examples generate.simple.sample.outcomes.matrix()
generate.simple.sample.outcomes.matrix <- function(formula, data, pheno.id="SUBJECT.NAME", 
                                                   method=c("bootstrap", "permutation"), use.REML=TRUE, 
                                                   num.samples, seed=1){
  method <- method[1]
  
  use.data <- model.frame(formula(paste0(paste0(Reduce(paste, deparse(formula))), "+", pheno.id)), data=data)
  names(use.data)[1] <- "y"
  set.seed(seed)
  if(method == "bootstrap"){
    null.formula.string <- paste0("y ~", unlist(strsplit(paste0(Reduce(paste, deparse(formula))), split="~"))[-1])
    fit0 <- lm(formula(null.formula.string), data=use.data)
    num.samples <- length(fit0$residuals)
    sigma2 <- ifelse(use.REML, sum(fit0$residuals^2)/(num.samples - 1), sum(fit0$residuals^2)/num.samples)
    new.y <- data.frame(sapply(1:num.samples, function(x) fit0$fitted.values + rnorm(n=num.samples, mean=0, sd=sqrt(sigma2))))
  }
  if(method == "permutation"){
    new.y <- data.frame(sapply(1:num.samples, function(x) sample(use.data$y)))
  }
  names(new.y) <- paste0("y", 1:num.samples)
  new.y$SUBJECT.NAME <- use.data$SUBJECT.NAME
  return(new.y)
}

#' @export
instability.lm.scan <- function(formula, data, num.bs.scans=100,
                                genomecache,
                                model=c("additive", "full"),
                                seed=1,
                                pheno.id="SUBJECT.NAME",
                                use.ROP=TRUE, num.imp=10, chr="all", just.these.loci=NULL, print.locus.fit=TRUE,
                                ...){

  h <- DiploprobReader$new(genomecache)
  loci <- h$getLoci()
  if(!is.null(just.these.loci)){
    loci <- just.these.loci
  }
  
  full.results <- matrix(NA, nrow=num.bs.scans, ncol=length(loci))
  colnames(full.results) <- loci
  these.chr <- h$getChromOfLocus(loci)
  these.pos <- list(cM=h$getLocusStart(loci, scale="cM"),
                    Mb=h$getLocusStart(loci, scale="Mb"))
  founders <- gsub(pattern="/", replacement=".", x=h$getFounders(), fixed=TRUE)
  ##### closures (functions defined within function)
  get.qr.null <- function(){
    qr.null <- qr.alt
    num.gen.columns <- ifelse(model=="additive", length(founders)-1, length(founders)+choose(length(founders)-1, 2))
    total.columns <- ncol(design.matrix)
    qr.null$qr <- qr.null$qr[,1:(total.columns-num.gen.columns), drop=FALSE]
    qr.null$rank <- length(1:(total.columns-num.gen.columns))
    qr.null$qraux <- qr.null$qraux[1:(total.columns-num.gen.columns)]
    qr.null$pivot <- qr.null$pivot[1:(total.columns-num.gen.columns)]
    return(qr.null)
  }

  run.locus.bs.fits <- function(){
    get.f.stat.p.val <- function(qr.alt, qr.null, y){
      rss0 <- sum(qr.resid(qr.null, y)^2)
      rss1 <- sum(qr.resid(qr.alt, y)^2)
      df1 <- qr.alt$rank - qr.null$rank
      df2 <- num.samples - qr.alt$rank
      
      mst <- (rss0 - rss1)/df1
      mse <- rss1/df2
      f.stat <- mst/mse
      p.val <- pf(q=f.stat, df1=df1, df2=df2, lower.tail=FALSE)
      return(p.val)
    }
    this.locus <- rep(NA, num.bs.scans)
    for(bs in 1:num.bs.scans){
      this.bs <- data.frame(y=new.y[,paste0("X", bs)], SUBJECT.NAME=new.y$SUBJECT.NAME)
      bs.data <- merge(this.bs, full.real.data[,-2], by="SUBJECT.NAME", all=FALSE)

      this.locus[bs] <- get.f.stat.p.val(qr.alt=qr.alt, qr.null=qr.null, y=bs.data$y)
      if(print.locus.fit){
        if(use.ROP){
          cat(paste0("bs: ", bs, ", locus: ", i, " of ", length(loci), "\n"))
        }
        if(!use.ROP){
          cat(paste0("imp: ", imp, ", bs: ", bs, ", locus: ", i, " of ", length(loci), "\n"))
        }
      }
    }
    return(this.locus)
  }
  impute.results <- array(NA, dim=c(num.imp, num.bs.scans, length(loci)))
  for(i in 1:length(loci)){
    this.locus.scan <- rep(NA, num.bs.scans)
    if(use.ROP){
      # ROP
      geno.data <- h$getLocusMatrix(loci[i], model=model)
      geno.names <- gsub(pattern="/", replacement=".", x=colnames(geno.data), fixed=TRUE)
      geno.data <- data.frame(SUBJECT.NAME=rownames(geno.data), geno.data)
      full.real.data <- merge(x=use.data, y=geno.data, by="SUBJECT.NAME", all=FALSE)
      set.to.intercept <- which.max(colSums(full.real.data[,geno.names]))
      full.real.data <- full.real.data[,-which(names(full.real.data) == geno.names[set.to.intercept])]
      alt.formula.string <- paste0(null.formula.string, " + ", paste(geno.names[-set.to.intercept], collapse=" + "))
      design.matrix <- model.matrix(formula(alt.formula.string), data=full.real.data)
      
      qr.alt <- qr(design.matrix)
      # get qr.null from qr.alt
      qr.null <- get.qr.null()
      this.locus <- run.locus.bs.fits()
      full.results[,i] <- this.locus
      impute.results <- NULL
    }
    if(!use.ROP){
      set.seed(seed)
      prob.geno.data <- h$getLocusMatrix(loci[i], model="full")
      prob.geno.data[prob.geno.data < 0] <- 0
      for(imp in 1:num.imp){
        geno.names <- colnames(prob.geno.data)
        geno.data <- t(apply(prob.geno.data, 1, function(x) rmultinom(1, 1, x)))
        colnames(geno.data) <- geno.names
        if(model == "additive"){
          full.to.dosages <- straineff.mapping.matrix()
          geno.data <- geno.data %*% full.to.dosages
          colnames(geno.data) <- founders
          geno.names <- founders
        }
        geno.data <- data.frame(SUBJECT.NAME=rownames(geno.data), geno.data)
        full.real.data <- merge(x=use.data, y=geno.data, by="SUBJECT.NAME", all=FALSE)
        set.to.intercept <- which.max(colSums(full.real.data[,geno.names]))
        full.real.data <- full.real.data[,-which(names(full.real.data) == geno.names[set.to.intercept])]
        alt.formula.string <- paste0(null.formula.string, " + ", paste(geno.names[-set.to.intercept], collapse=" + "))
        design.matrix <- model.matrix(formula(alt.formula.string), data=full.real.data)
        qr.alt <- qr(design.matrix)
        qr.null <- get.qr.null()
        this.locus <- run.locus.bs.fits()
        impute.results[imp,,i] <- this.locus
      }
      full.results[,i] <- apply(impute.results[,,i, drop=FALSE], 2, function(x) median(x))
    }
  }
  return(list(impute.results=impute.results, full.results=list(p.values=full.results, chr=these.chr, pos=these.pos),
              formula=paste0(Reduce(paste, deparse(formula))),
              model=model))
}



sim.instability.lm.scan <- function(formula, data, num.bs.scans=100,
                                    genomecache,
                                    model=c("additive", "full"),
                                    num.founders=2,
                                    seed=1,
                                    scale="cM",
                                    use.ROP=TRUE, num.imp=10,
                                    dir.par=1,
                                    num.sim=1000,
                                    ...){
  
  use.data <- model.frame(formula(paste0(paste0(Reduce(paste, deparse(formula))), "+ SUBJECT.NAME")), data=data)
  names(use.data)[1] <- "y"
  null.formula.string <- paste0("y ~", unlist(strsplit(paste0(Reduce(paste, deparse(formula))), split="~"))[-1])
  fit0 <- lm(formula(null.formula.string), data=use.data)
  sigma2 <- sum(fit0$residuals^2)/(length(fit0$residuals) - 1)
  set.seed(seed)
  new.y <- data.frame(sapply(1:num.bs.scans, function(x) fit0$coefficients + rnorm(n=num.sim, mean=0, sd=sqrt(sigma2))))
  SUBJECT.NAME <- paste0("sim.", 1:num.sim)
  
  h <- DiploprobReader$new(genomecache)
  loci <- h$getLoci()
  
  full.results <- matrix(NA, nrow=num.bs.scans, ncol=length(loci))
  colnames(full.results) <- loci
  these.chr <- h$getChromOfLocus(loci)
  these.pos <- h$getLocusStart(loci, scale=scale)
  
  founders <- LETTERS[1:num.founders]
  
  ##### closures (functions defined within function)
  get.qr.null <- function(){
    qr.null <- qr.alt
    num.gen.columns <- ifelse(model=="additive", length(founders), length(founders)+choose(length(founders), 2))
    total.columns <- ncol(design.matrix)
    qr.null$qr <- qr.null$qr[,1:(total.columns-num.gen.columns), drop=FALSE]
    qr.null$rank <- length(1:(total.columns-num.gen.columns))
    qr.null$qraux <- qr.null$qraux[1:(total.columns-num.gen.columns)]
    qr.null$pivot <- qr.null$pivot[1:(total.columns-num.gen.columns)]
    return(qr.null)
  }
  
  run.locus.bs.fits <- function(){
    get.f.stat.p.val <- function(qr.alt, qr.null, y){
      rss0 <- sum(qr.resid(qr.null, y)^2)
      rss1 <- sum(qr.resid(qr.alt, y)^2)
      df1 <- qr.alt$rank - qr.null$rank
      df2 <- num.sim - qr.alt$rank
      
      mst <- (rss0 - rss1)/df1
      mse <- rss1/df2
      f.stat <- mst/mse
      p.val <- pf(q=f.stat, df1=df1, df2=df2, lower.tail=FALSE)
      return(p.val)
    }
    this.locus <- rep(NA, num.bs.scans)
    for(bs in 1:num.bs.scans){
      this.locus[bs] <- get.f.stat.p.val(qr.alt=qr.alt, qr.null=qr.null, y=new.y[,bs])
      if(use.ROP){
        cat(paste0("bs: ", bs, ", locus: ", i, " of ", length(loci), "\n"))
      }
      if(!use.ROP){
        cat(paste0("imp: ", imp, ", bs: ", bs, ", locus: ", i, " of ", length(loci), "\n"))
      }
    }
    return(this.locus)
  }
  
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
  
  full.to.dosages <- straineff.mapping.matrix(M=num.founders)
  dip.names <- as.vector(apply(full.to.dosages, 1, function(x) paste(rev(LETTERS[1:num.founders][which(x == 1)]), collapse=".")))
  dip.names[1:num.founders] <- paste(LETTERS[1:num.founders], LETTERS[1:num.founders], sep=".")

  alt.formula.string <- ifelse(model=="additive",
                               paste0(null.formula.string, " + ", paste(founders, collapse="+")),
                               paste0(null.formula.string, " + ", paste(dip.names, collapse="+")))
  
  set.seed(seed)
  for(i in 1:length(loci)){
    this.locus.scan <- rep(NA, num.bs.scans)
    
    require(gtools)
    pre.geno.data <- matrix(NA, nrow=num.sim, ncol=num.founders + choose(num.founders, 2))
    dip.row <- 1
    while(dip.row <= num.sim){
      this.row <- rdirichlet(n=1, rep(dir.par, num.founders + choose(num.founders, 2)))
      if(!(any(is.nan(this.row)) | any(is.na(this.row)) | any(is.infinite(this.row)))){
        pre.geno.data[dip.row,] <- this.row
        dip.row <- dip.row + 1
      }
    }
    #pre.geno.data <- rdirichlet(n=num.sim, rep(dir.par, num.founders + choose(num.founders, 2)))
    colnames(pre.geno.data) <- dip.names
    
    if(use.ROP){
      # ROP
      if(model == "additive"){
        geno.data <- pre.geno.data %*% full.to.dosages
        colnames(geno.data) <- founders
      }
      design.matrix <- cbind(rep(1, num.sim), geno.data)
      colnames(design.matrix)[1] <- "(Intercept)"
      qr.alt <- qr(design.matrix)
      # get qr.null from qr.alt
      qr.null <- get.qr.null()
      this.locus <- run.locus.bs.fits()
      full.results[,i] <- this.locus
    }
    if(!use.ROP){
      set.seed(seed)
      impute.results <- matrix(NA, nrow=num.imp, ncol=num.bs.scans)
      for(imp in 1:num.imp){
        repeat.this <- TRUE
        while(repeat.this){
          geno.data <- t(apply(pre.geno.data, 1, function(x) rmultinom(1, 1, x)))
          if(!(any(is.na(geno.data)) | any(is.nan(geno.data)) | any(is.infinite(geno.data)))){
            repeat.this <- FALSE
          }
        }
        colnames(geno.data) <- dip.names
        if(model == "additive"){
          #geno.data <- geno.data %*% full.to.dosages #has issues for some reason (makes NaN's occasionally)
          geno.data <- t(apply(geno.data, 1, function(x) x %*% full.to.dosages))
          colnames(geno.data) <- founders
        }
        design.matrix <- cbind(rep(1, num.sim), geno.data)
        colnames(design.matrix)[1] <- "(Intercept)"
        qr.alt <- qr(design.matrix)
        qr.null <- get.qr.null()
        this.locus <- run.locus.bs.fits()
        #browser()
        #this.ci <- eexp(-log(this.locus), ci=TRUE, conf.level=0.95)$interval$limits[c("LCL", "UCL")]
        impute.results[imp,] <- this.locus
      }
      full.results[,i] <- apply(impute.results, 2, function(x) median(x))
    }
  }    
  return(list(full.results=list(p.values=full.results, chr=these.chr, pos=these.pos),
                formula=paste0(Reduce(paste, deparse(formula))),
                model=model,
                num.sim=num.sim,
                num.founders=num.founders))
}
    
    
# temp.scan <- instability.lm.scan(log(BMI_Butt)~1+Dissector, data=Solberg.pheno.data, num.bs.scans=1000,
#                                 genomecache="~/Documents/SolbergHS/caches/DOQTL_genomecache_989rats/",
#                                 model="additive",
#                                 seed=5,
#                                 scale="cM",
#                                 use.ROP=TRUE, num.imp=10, just.these.loci=c("WKYc89c09_s1_485", "WKYc69b04_s1_762"))
# temp.scan <- instability.lm.scan(log(BMI_Butt)~1, data=Solberg.pheno.data, num.bs.scans=1000,
#                                  genomecache="~/Documents/SolbergHS/caches/DOQTL_genomecache_989rats/",
#                                  model="additive",
#                                  seed=1,
#                                  scale="cM",
#                                  use.ROP=TRUE, num.imp=10, just.these.loci=c("WKYc89c09_s1_485", "WKYc69b04_s1_762"))
# 
# BMI.woDissector.ROP.instability.scan <- instability.lm.scan(log(BMI_Butt)~1, data=Solberg.pheno.data, num.bs.scans=1000,
#                                  genomecache="~/Documents/SolbergHS/caches/DOQTL_genomecache_989rats/",
#                                  model="additive",
#                                  seed=1,
#                                  scale="cM",
#                                  use.ROP=TRUE)
# BMI.wDissector.ROP.instability.scan <- instability.lm.scan(log(BMI_Butt)~1+Dissector, data=Solberg.pheno.data, num.bs.scans=1000,
#                                                            genomecache="~/Documents/SolbergHS/caches/DOQTL_genomecache_989rats/",
#                                                            model="additive",
#                                                            seed=1,
#                                                            scale="cM",
#                                                            use.ROP=TRUE)
# 
# cairo_pdf("~/Documents/SolbergHS/instability/balance/BMI_comparison_instability.pdf", height=7, width=7)
# par(mfrow=c(2,1))
# null.pvalue.ci.plot(par.bs.scans=BMI.woDissector.ROP.instability.scan, cache.title="Solberg DO-QTL genomecache")
# null.pvalue.ci.plot(par.bs.scans=BMI.wDissector.ROP.instability.scan, cache.title="Solberg DO-QTL genomecache")
# dev.off()
# 
# BMI.ROP.instability.scan <- instability.lm.scan(sqrt(Ins0)~1+number+Handler, data=Solberg.pheno.data, num.bs.scans=1000,
#                                                 genomecache="~/Documents/SolbergHS/caches/DOQTL_genomecache_989rats/",
#                                                 model="additive",
#                                                 seed=1,
#                                                 scale="cM",
#                                                 use.ROP=TRUE, num.imp=10, just.these.loci=c("WKYc89c09_s1_485", "WKYc69b04_s1_762"))