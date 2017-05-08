#' Run a SNP-based genome scan from probabilities stored in a genome cache directory and founder strain
#' alleles stored .allele files
#'
#' This function primarily takes a formula, data frame, genome cache, and directory of founder strain
#' .alleles files to run a genome scan.
#'
#' @param data A data frame with outcome and potential covariates. Should also have individual IDs
#' that link to IDs in the genome cache  with a column named "SUBJECT.NAME".
#' @param formula An lm style formula with functions of outcome and covariates contained in data frame.
#' @param K DEFAULT: NULL. A positive semi-definite relationship matrix, usually a realized genetic relationship matrix (GRM)
#' based on SNP genotypes or the founder haplotype probabilities. Colnames and rownames should match
#' the SUBJECT.NAME column in the data frame. If no K matrix is specified, either lmer is used (if sparse random effects
#' are included in the formula) or a fixed effect model (equivalent to lm).
#' @param allele.dir The path to the directory of .allele files that specify which SNP alleles correspond
#' to which founder haplotype. .allele files for a format used by HAPPY.
#' @param genomecache The path to the genome cache directory. The genome cache is a particularly structured
#' directory that stores the haplotype probabilities/dosages at each locus. It has an additive model
#' subdirectory and a full model subdirectory. Each contains subdirectories for each chromosome, which then
#' store .RData files for the probabilities/dosages of each locus.
#' @param model DEFAULT: additive. Specifies how to model the founder haplotype probabilities. The additive options specifies
#' use of SNP dosages, and is most commonly used. The full option regresses the phenotype on the actual
#' genotype probabilities.
#' @param chr DEFAULT: "all". Specifies which chromosomes to scan.
#' @param brute DEFAULT: TRUE. During the optimization to find maximum likelihood parameters, this specifies checking the
#' boundaries of h2=0 and h2=1. Slightly less efficient, but otherwise the optimization procedure will not directly check
#' these values.
#' @param use.fix.par DEFAULT: TRUE. This specifies an approximate fitting of mixed effect model (Kang et al. 2009). Much
#' more efficient, as the optimization of h2 only needs to be performed once for the null model rather than every locus. 
#' Technically less powerful, though in practice it has proven to be almost equal to the exact procedure.
#' @param just.these.loci DEFAULT: NULL. Specifies a reduced set of loci to fit.
#' @param print.locus.fit DEFAULT: FALSE. If TRUE, prints out how many loci have been fit currently.
#' @param exclusion.freq DEFAULT: .Machine$double.eps. Loci with observed minor allele frequencies beneath
#' the specified value are removed from the scan.
#' @param X.list DEFAULT: NULL. This specifies the SNP-based design matrices for all the loci. If a scan
#' of the same population with the same markers has been performed, this option can save a lot of time.
#' @param return.X.list DEFAULT: FALSE. The scan procedure can return the list of design matrices for all
#' loci. It does increase the size of the output, though can then be used as an input for subsequent scans
#' of the same data.
#' @export
#' @examples imputed.snp.scan.h2lmm()
imputed.snp.scan.h2lmm <- function(data, formula, K,
                                   allele.dir, genomecache,
                                   model=c("additive", "full"),
                                   use.par="h2", chr="all", brute=TRUE, use.fix.par=FALSE, 
                                   just.these.loci=NULL,
                                   print.locus.fit=FALSE,
                                   exclusion.freq=.Machine$double.eps,
                                   X.list=NULL, return.X.list=FALSE, # Makes multiple scans more efficient
                                   ...){
  model <- model[1]
  
  if(is.null(X.list)){ rm("X.list") }
  
  mapping.matrix <- straineff.mapping.matrix()

  these.chr <- chr
  if(chr == "all"){ these.chr <- c(1:20, "X") }
  # Genomecache for imputation
  h <- DiploprobReader$new(genomecache)
  
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
  
  data.and.K <- make.processed.data(formula=formula, data=data, K=K, cache.subjects=cache.subjects, impute.on="SUBJECT.NAME")
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
    data <- cbind(null.data, X[as.character(null.data$SUBJECT.NAME),, drop=FALSE])
    fit1 <- lmmbygls(locus.formula, data=data, 
                     eigen.K=fit0$eigen.K, K=fit0$K, 
                     use.par="h2", fix.par=fix.par, M=fit0$M, logDetV=fit0$logDetV,
                     brute=brute, 
                     weights=NULL)
    LOD.vec[i] <- log10(exp(fit1$logLik - fit0$logLik))
    p.vec[i] <- pchisq(q=-2*(fit0$logLik - fit1$logLik), df=fit1$rank-fit0$rank, lower.tail=FALSE)
    h2.record[i] <- fit1$h2
    if(print.locus.fit){ cat(paste("locus", i, "out of", length(loci)), "\n") }
  }
  if(!return.X.list){ X.list <- NULL }
  output <- list(LOD=LOD.vec,
                 p.value=p.vec,
                 pos=list(Mb=h$getMarkerLocation(loci, scale="Mb"), cM=h$getMarkerLocation(loci, scale="cM")),
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

