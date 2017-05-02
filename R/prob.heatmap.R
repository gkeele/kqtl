#' Plot founder haplotype dosage/probabilities, ordered by phenotype
#'
#' This function produces a probability heatmap plot, ordered by the phenotype. This gives an idea of what the regression
#' procedure is actually being handed as inputs.
#'
#' @param marker A marker that is contained in the genome cache. In general this should be a marker of interest, such as one
#' beneath a putative QTL peak.
#' @param p.value DEFAULT: NULL. Includes the observed p-value in the plot title.
#' @param genomecache The path to the genome cache that contains founder haplotype information.
#' @param model DEFAULT: "additive". If "additive", dosages are plotted. If "full", probabilities are plotted.
#' @param phenotype The name of the phenotype column in the data set.
#' @param phenotype.data A data.frame object that contains the phenotype information. Should also have a column that matches
#' genomes in the genome cache.
#' @param merge.by DEFAULT: "SUBJECT.NAME". Specifies the columns to merge phenotype and haplotype data.
#' @param founder.labels DEFAULT: "NULL". If NULL, will default to the labels in the genome cache.
#' @export
#' @examples prob.heatmap()
prob.heatmap = function(marker, p.value=NULL, genomecache, model="additive",
                        phenotype, phenotype.data, merge.by="SUBJECT.NAME", founder.labels=NULL){
  if(!is.null(p.value)){
    p.value <- round(-log10(p.value), 4)
  }
  
  h <- DiploprobReader$new(genomecache)
  X <- h$getLocusMatrix(locus=marker, model=model)
  if(is.null(founder.labels)){
    founder.labels <- colnames(X)
  }
  num.col <- length(founder.labels)
  # if(model == "additive"){
  #   if(sum(X[1,]) == 2){
  #     X <- X/2
  #   }
  # }
  
  subjects <- h$getSubjects()
  X.data <- data.frame(rownames(X), X)
  names(X.data)[1] <- merge.by
  
  # allow function of phenotype
  phenotype.data <- model.frame(formula(paste(phenotype, "~ 1 +", merge.by)), data=phenotype.data)
  names(phenotype.data)[1] <- "y"
  final.data <- merge(x=phenotype.data, y=X.data, by=merge.by, all=FALSE)

  final.data <- final.data[order(final.data$y),] # sort by phenotypic value
  probs <- as.matrix(final.data[,-(1:2)]) # Just keep prob of 8 strains for a certain marker
  probs <- probs[,rev(1:ncol(probs))]

  s <- summary(final.data[,"y"])
  s1 <- as.character(round(s[1], 2))
  s2 <- as.character(round(s[2], 2))
  s3 <- as.character(round(s[3], 2))
  s5 <- as.character(round(s[5], 2))
  s6 <- as.character(round(s[6], 2))
  
  # plot heatmap
  ## save original par settings
  op <- par()
  cols <- rev(gray(10000:1/10000))
  par(plt=c(0.1,.75,0.1,.8))    ##set the margin  
  image(z=1-probs, axes=FALSE, col=cols)
  #heatmap(t(probs), Rowv=NA, Colv=NA, col = rev(cols), scale="column", margins=c(4,4), labCol="") 
  box()
  axis(2, at=seq(0, num.col, 1+1/num.col)/num.col, 
       labels=rev(founder.labels), 
       lty=0, srt=90, las=2) # add txt on the strain   
  axis(1, at=0.5, labels=phenotype, tick=FALSE)
  axis(3, at=c(0,0.25,0.5,0.75,1), labels=c(s1,s2,s3,s5,s6))
  this.title <- ifelse(is.null(p.value), marker, paste0(marker, ": -log10P=", p.value))
  title(this.title, line=2.5)
  
  ramp.label <- ifelse(model == "additive", "Dos", "Prob")
  par(fig=c(.9,.97,.3,.6), mai=c(.1,.1,.3,.1), new=TRUE)
  if(model == "additive"){ image(y=seq(from=0, to=2, length.out=length(cols)), z=matrix(seq(from=0, to=2, length.out=length(cols)), nrow=1), 
                                 zlim=c(0, 2), ylim=c(0, 2), axes=FALSE, col=rev(cols), main=ramp.label) } #for the legend 
  if(model == "full"){ image(y=seq(from=0, to=1, length.out=length(cols)), z=matrix(seq(from=0, to=1, length.out=length(cols)), nrow=1), 
                             zlim=c(0, 1), ylim=c(0, 1), axes=FALSE, col=rev(cols), main=ramp.label) }
  box()
  axis(2, las=1)          
  
  par(op)
}

#' @export
prob.image = function(marker.data, marker=NULL, p.value=NULL, 
                      phenotype, phenotype.data, column.labels=NULL){
  if(!is.null(p.value)){ p.value <- round(-log10(p.value), 4) }
  if(is.null(column.labels)){ column.labels <- colnames(marker.data)[-1] }
  if(is.null(marker)){ marker <- "locus" }
  num.col <- length(column.labels)
  if(sum(marker.data[1,-1], na.rm=TRUE) != 1){
    marker.data[,-1] <- marker.data[,-1]/2
  }
  # allow function of phenotype
  phenotype.data <- model.frame(formula(paste0(phenotype, " ~ 1 + SUBJECT.NAME")), data=phenotype.data)
  names(phenotype.data)[1] <- "y"
  final.data <- merge(x=phenotype.data, y=marker.data, by="SUBJECT.NAME", all=FALSE)
  
  final.data <- final.data[order(final.data$y),] # sort by phenotypic value
  probs <- as.matrix(final.data[,-(1:2)]) # Just keep prob of 8 strains for a certain marker
  
  s <- summary(final.data[,"y"])
  s1 <- as.character(round(s[1], 2))
  s2 <- as.character(round(s[2], 2))
  s3 <- as.character(round(s[3], 2))
  s5 <- as.character(round(s[5], 2))
  s6 <- as.character(round(s[6], 2))
  
  # plot image
  ## save original par settings
  oplt <- par()$plt
  cols <- rev(gray(10000:1/10000))
  par(plt=c(0.1,.75,0.1,.8))    ##set the margin  
  image(z=1-probs, axes=F, col=cols, zlim=c(0, 1))
  box()
  axis(2, at=seq(0, num.col, 1+1/num.col)/num.col, 
       labels=column.labels, 
       lty=0, srt=90, las=2) # add txt on the strain   
  axis(1, at=0.5, labels=phenotype, tick=FALSE)
  axis(3, at=c(0,0.25,0.5,0.75,1), labels=c(s1,s2,s3,s5,s6))
  this.title <- ifelse(is.null(p.value), marker, paste0(marker, ": -log10P=", p.value))
  title(this.title, line=2.5)
  par(plt <- oplt)
}
