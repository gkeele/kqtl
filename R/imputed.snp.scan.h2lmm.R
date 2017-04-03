imputed.snp.scan.h2lmm <- function(data, formula, K,
                                   allele.dir, genomecache,
                                   model=c("additive", "full"),
                                   use.par="h2", chr="all", brute=TRUE, use.fix.par=FALSE, 
                                   use.chol=FALSE,
                                   just.these.loci=NULL,
                                   print.locus.fit=FALSE,
                                   map.file=NULL,
                                   exclusion.freq=.Machine$double.eps,
                                   X.list=NULL, return.X.list=FALSE, # Makes multiple scans more efficient
                                   ...){
  model <- model[1]
  
  if(is.null(X.list)){ rm("X.list") }
  
  mapping.matrix <- straineff.mapping.matrix()

  these.chr <- chr
  if(chr == "all"){ these.chr <- c(1:20, "X") }
  # Genomecache for imputation
  h <- bagpipe.backend::DiploprobReader$new(genomecache)
  
  if(!exists("X.list")){
    loci <- h$getLoci(chr=these.chr)
    loci <- loci[-grep(pattern="^c[X0-9]+\\.loc", x=loci, perl=TRUE)] # Remove pseudomarkers created by qtl2geno
  }
  else{
    loci <- names(X.list)
  }
  founders <- h$getFounders()
  cache.subjects <- rownames(h$getLocusMatrix(locus=loci[1], model="additive"))
  loci.chr <- h$getChromOfLocus(loci)
  
  data.and.K <- make.processed.data(formula=formula, data=data, K=K, cache.subjects=cache.subjects)
  data <- data.and.K$data
  K <- data.and.K$K
  
  if(!is.null(just.these.loci)){
    loci <- loci[loci %in% just.these.loci]
    loci.chr <- loci.chr[loci %in% just.these.loci]
  }
  
  formula.string <- Reduce(paste, deparse(formula))
  null.formula <- make.null.formula(formula=formula, do.augment=FALSE)
  locus.formula <- make.snp.alt.formula(formula=formula)
  original.n <- nrow(data)
  
  ## Fitting null model
  eigen.K <- eigen(K)
  fit0 <- lmmbygls(null.formula, data=data, eigen.K=eigen.K, K=K, use.par=use.par, weights=NULL, brute=brute)
  fit0.REML <- lmmbygls(null.formula, data=data, eigen.K=eigen.K, K=K, use.par="h2.REML", weights=NULL, brute=brute)
  
  if(use.fix.par){
    fix.par <- fit0$h2
  }
  if(!use.fix.par){
    fix.par <- NULL
  }
  
  null.data <- data
  null.K <- K
  
  if(!exists("X.list")){
    X.list <- make.imputed.design.matrix.list.for.all.loci(loci=loci, loci.chr=loci.chr, n=nrow(data), model=model, h=h, 
                                                           allele.dir=allele.dir, mapping.matrix=mapping.matrix,
                                                           founders=founders, exclusion.freq=exclusion.freq)
    keep.loci <- loci %in% names(X.list)
    loci <- loci[keep.loci]
    loci.chr <- loci.chr[keep.loci]
  }
  LOD.vec <- p.vec <- h2.record <- rep(0, length(loci))
  for(i in 1:length(loci)){
    X <- X.list[[i]]
    data <- cbind(null.data, X[null.data$SUBJECT.NAME,, drop=FALSE])
    fit1 <- lmmbygls(locus.formula, data=data, 
                     eigen.K=fit0$eigen.K, K=fit0$K, 
                     use.par="h2", fix.par=fix.par, M=fit0$M, logDetV=fit0$logDetV,
                     brute=brute, 
                     weights=NULL, calc.h2.range.by=NULL, use.chol=use.chol)
    LOD.vec[i] <- log10(exp(fit1$logLik - fit0$logLik))
    p.vec[i] <- pchisq(q=-2*(fit0$logLik - fit1$logLik), df=fit1$rank-fit0$rank, lower.tail=FALSE)
    h2.record[i] <- fit1$h2
    if(print.locus.fit){
      cat("locus", i, "\n")
    }
  }
  if(is.null(map.file)){
    pos <- NULL
  }
  else{
    map.data <- read.table(map.file, header=TRUE)
    rownames(map.data) <- map.data$marker
    pos <- list(cM=map.data[loci,]$pos, Mb=map.data[loci,]$bp/1000000)
  }
  if(!return.X.list){ X.list <- NULL }
  output <- list(LOD=LOD.vec,
                 p.value=p.vec,
                 pos=pos,
                 loci=loci, 
                 chr=loci.chr,
                 h2=h2.record,
                 fit0=fit0,
                 fit0.REML=fit0.REML,
                 formula=formula.string,
                 model.type=model,
                 X.list=X.list,
                 null.data=null.data)
  return(output)
}

read.in.genotype.dir <- function(genotype.dir, genotype.file.prefix, chr){
  if(chr == "all"){ chr <- c(1:20, "X") }
  for(i in 1:length(chr)){
    this.genotype.data <- read.table(paste(genotype.dir, paste(genotype.file.prefix, chr[i], "geno", sep="."), sep="/"), header=TRUE)
    colnames(this.genotype.data)[1] <- "SUBJECT.NAME"
    if(i == 1){
      genotype.data <- this.genotype.data
      chr.vec <- rep(chr[i], ncol(this.genotype.data)-1)
    }
    if(i > 1){
      genotype.data <- merge(genotype.data, this.genotype.data, by="SUBJECT.NAME")
      chr.vec <- c(chr.vec, rep(chr[i], ncol(this.genotype.data)-1))
    }
  }
  return(list(genotype.data=genotype.data, chr=chr.vec))
}

extract.imputed.design.matrix.from.doqtl.genotype <- function(probs, allele.dir, 
                                                              snp, snp.chr, model, 
                                                              founders, mapping.matrix){
  if(model == "additive"){
    grep.command <- paste0("grep -A 3 '", snp, "' ", 
                           paste0(allele.dir, "/chr", snp.chr, ".alleles"))
    founder.alleles.table <- system(grep.command, intern=TRUE)
    founder.alleles.table <- matrix(unlist(strsplit(x=founder.alleles.table[-1], split="\t", fixed=TRUE)), nrow=3, byrow=TRUE)[,-1] # Remove column of "allele"
    founder.alleles.table <- founder.alleles.table[founder.alleles.table[,1] != "NA",] # Remove NA row
    founder.alleles <- founder.alleles.table[,1][apply(founder.alleles.table[,-1], 2, function(x) which.max(x))]
    ref.allele <- founder.alleles.table[1,1]
    ref.allele.founder.count <- as.numeric(founder.alleles == ref.allele)
    full.ref.allele.count <- as.vector(mapping.matrix %*% matrix(ref.allele.founder.count, ncol=1))
    snp.count <- probs %*% matrix(full.ref.allele.count, ncol=1)
    # Converting to count of minor allele
    if(sum(snp.count)/(2*length(snp.count)) > 0.5){
      snp.count <- 2 - snp.count
    }
    X <- matrix(snp.count, ncol=1)
    rownames(X) <- rownames(probs)
    colnames(X) <- "SNP"
  }
  return(X)
}

make.imputed.design.matrix.list.for.all.loci <- function(loci, loci.chr, n, model, h, 
                                                         allele.dir, mapping.matrix,
                                                         founders, exclusion.freq){
  p <- length(loci)
  if(model=="additive"){
    X.list <- rep(list(matrix(NA, nrow=n, ncol=1)), p)
    for(i in 1:p){
      probs <- h$getLocusMatrix(locus=loci[i], model="full")
      X.list[[i]] <- extract.imputed.design.matrix.from.doqtl.genotype(probs=probs, allele.dir=allele.dir, 
                                                                       snp=loci[i], snp.chr=loci.chr[i], model=model, 
                                                                       founders=founders, mapping.matrix=mapping.matrix)
      rownames(X.list[[i]]) <- rownames(probs)
      colnames(X.list[[i]]) <- "SNP"
    }
    names(X.list) <- loci
    
    # Removing low frequency alleles
    keep.loci.index <- unlist(lapply(X.list, function(x) sum(x)/(2*length(x)))) > exclusion.freq
    X.list[which(!keep.loci.index)] <- NULL
    return(X.list)
  }
}

