pairwise.cor.snp.scan <- function(data, formula, K,
                                  allele.dir, genomecache,
                                  model=c("additive", "full"), chr, point.locus,
                                  just.these.loci=NULL,
                                  print.progress=FALSE,
                                  map.file=NULL,
                                  exclusion.freq=.Machine$double.eps,
                                  X.list,
                                  ...){
  model <- model[1]
  
  do.augment <- FALSE # Necessary to work with old support_functions.R
  
  mapping.matrix <- straineff.mapping.matrix()
  
  these.chr <- chr
  # Genomecache for imputation
  h <- bagpipe.backend::DiploprobReader$new(genomecache)
  
  loci <- names(X.list)
  
  founders <- h$getFounders()
  cache.subjects <- rownames(h$getLocusMatrix(locus=loci[1], model="additive"))
  loci.chr <- h$getChromOfLocus(loci)
  
  loci <- loci[loci.chr == these.chr]
  loci.chr <- loci.chr[loci.chr == these.chr]
  
  data.and.K <- make.processed.data(formula=formula, data=data, K=K, cache.subjects=cache.subjects)
  data <- data.and.K$data
  
  if(!is.null(just.these.loci)){
    loci <- loci[loci %in% just.these.loci]
    loci.chr <- loci.chr[loci %in% just.these.loci]
  }
  
  null.data <- data
  null.K <- K
  
  # Grabbing the design matrix (vector) of the point SNP
  point.X <- X.list[[point.locus]][null.data$SUBJECT.NAME,, drop=FALSE]
  
  r2.vec <- rep(0, length(loci))
  for(i in 1:length(loci)){
    X <- X.list[[loci[i]]][null.data$SUBJECT.NAME,, drop=FALSE]
    
    r2.vec[i] <- cor(point.X, X)^2
    
    if(print.progress){
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
  output <- list(r2=r2.vec,
                 point.locus=point.locus,
                 pos=pos,
                 loci=loci, 
                 chr=chr,
                 model.type=model)
  return(output)
}

extract.r2.interval <- function(scan.object, r2.scan.object, r2.level=0.6){
  high.r2 <- r2.scan.object$r2[which(r2.scan.object$r2 > r2.level)]
  high.r2.loci <- r2.scan.object$loci[which(r2.scan.object$r2 > r2.level)]
  high.r2.pvalue <- scan.object$p.value[which(scan.object$loci %in% high.r2.loci)]
  high.r2.pos.Mb <- scan.object$pos$Mb[which(scan.object$loci %in% high.r2.loci)]
  high.r2.pos.cM <- scan.object$pos$cM[which(scan.object$loci %in% high.r2.loci)]
  max.pos.index <- which.max(high.r2.pos.cM)
  min.pos.index <- which.min(high.r2.pos.cM)
  ## Interval
  interval.width.cM <- high.r2.pos.cM[max.pos.index] - high.r2.pos.cM[min.pos.index]
  interval.width.Mb <- high.r2.pos.Mb[max.pos.index] - high.r2.pos.Mb[min.pos.index]
  ## Boundary markers
  ub.marker <- high.r2.loci[max.pos.index]
  lb.marker <- high.r2.loci[min.pos.index]
  ub.cM <- high.r2.pos.cM[max.pos.index]
  lb.cM <- high.r2.pos.cM[min.pos.index]
  ub.Mb <- high.r2.pos.Mb[max.pos.index]
  lb.Mb <- high.r2.pos.Mb[min.pos.index]
  
  results <- data.frame(peak.marker=r2.scan.object$point.locus, r2.level=r2.level,
                        ub.marker=ub.marker, ub.cM=round(ub.cM, 2), ub.Mb=round(ub.Mb, 2),
                        lb.marker=lb.marker, lb.cM=round(lb.cM, 2), lb.Mb=round(lb.Mb, 2),
                        interval.width.cM=round(interval.width.cM, 2), interval.width.Mb=round(interval.width.Mb, 2))
  return(results)
}
