#' @export
prob.heatmap = function(marker, p.value=NULL, genomecache, model, 
                        phenotype, phenotype.data, column.labels=NULL){
  if(!is.null(p.value)){
    p.value <- round(-log10(p.value), 4)
  }
  
  h <- bagpipe.backend::DiploprobReader$new(genomecache)
  X <- h$getLocusMatrix(locus=marker, model=model)
  if(is.null(column.labels)){
    column.labels <- colnames(X)
  }
  num.col <- length(column.labels)
  if(model == "additive"){
    if(sum(X[1,]) == 2){
      X <- X/2
    }
  }
  
  subjects <- h$getSubjects()
  X.data <- data.frame(SUBJECT.NAME=rownames(X), X)

  # allow function of phenotype
  phenotype.data <- model.frame(formula(paste0(phenotype, " ~ 1 + SUBJECT.NAME")), data=phenotype.data)
  names(phenotype.data)[1] <- "y"
  final.data <- merge(x=phenotype.data, y=X.data, by="SUBJECT.NAME", all=FALSE)

  #final.data <- final.data[!is.na(final.data[,phenotype]),]
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
  image(z=1-probs, axes=F, col=cols)
  #heatmap(t(probs), Rowv=NA, Colv=NA, col = rev(cols), scale="column", margins=c(4,4), labCol="") 
  box()
  axis(2, at=seq(0, num.col, 1+1/num.col)/num.col, 
       labels=rev(column.labels), 
       lty=0, srt=90, las=2) # add txt on the strain   
  axis(1, at=0.5, labels=phenotype, tick=FALSE)
  axis(3, at=c(0,0.25,0.5,0.75,1), labels=c(s1,s2,s3,s5,s6))
  this.title <- ifelse(is.null(p.value), marker, paste0(marker, ": -log10P=", p.value))
  title(this.title, line=2.5)
  
  par(fig=c(.9,.97,.3,.6), mai=c(.1,.1,.3,.1), new=T)
  image(z=1-rbind(1:length(cols)), axes=F, col=cols, main="Prob") #for the legend
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
