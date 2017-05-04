#' Plot whole genome and single chromosome windows of haplotype-based genome scan in a PDF output document
#'
#' This function takes the genome scan output from scan.h2lmm() and plots the whole genome and single chromosome zoom-ins for 
#' all the specified chromosomes. When multiple imputations are used, includes the 95% confidence band on the median in the zoomed-in
#' plots.
#'
#' @param scan.object A scan.h2lmm() object (ROP or multiple imputations). If multiple imputations, median and confidence interval 
#' on median are plotted.
#' @param chr DEFAULT: c(1:19, "X"). The chromosomes to be plotted. DEFAULT is all the mouse chromosomes.
#' @param use.lod DEFAULT: TRUE. Plots either the LOD score or the -log10 p-value.
#' @param scale DEFAULT: "Mb". Specifies the scale of genomic position to be plotted. Either Mb or cM are expected.
#' @param main.colors DEFAULT: "black". The color of the main association score to be plotted.
#' @param median.band.col DEFAULT: "gray88". The color of the 95% confident band plotted around the median.
#' @param main DEFAULT: "". Adds a title above the model.
#' @param y.max.manual DEFAULT: NULL. Manually adds a max y-value. Allows multiple genome scans to easily be on the same scale.
#' @param hard.thresholds DEFAULT: NULL. Specify one or more horizontal threshold lines.
#' @param thresholds.col DEFAULT: "red". Set the colors of the specified thresholds.
#' @param pdf.output.path That path of the PDF file to be generated.
#' @param pdf.height DEFAULT: 5. The height of an individual pages of the PDF.
#' @param pdf.width DEFAULT: 9. The width of an individual pages of the PDF.
#' @export
#' @examples genome.plotter.to.pdf()
genome.plotter.to.pdf <- function(scan.object, chr=c(1:19, "X"), use.lod=TRUE,
                                  scale=c("Mb", "cM"), main.col="black", median.band.col="gray88", main="", y.max.manual=NULL,
                                  hard.thresholds=NULL, thresholds.col="red",
                                  pdf.output.path, pdf.height=5, pdf.width=9){
  scale <- scale[1]
  pdf(pdf.output.path, height=pdf.height, width=pdf.width)
  genome.plotter.whole(non.mi.scan.list=list(scan.object), use.lod=use.lod,
                       scale=scale, main.colors=main.col, use.legend=FALSE,
                       main=main,
                       y.max.manual=y.max.manual,
                       hard.thresholds=hard.thresholds, thresholds.col=thresholds.col)
  for(i in 1:length(chr)){
    genome.plotter.chr(scan.object=scan.object, chr=chr[i], use.lod=use.lod,
                       scale=scale, main.col=main.col, median.band.col=median.band.col,
                       main=main, y.max.manual=y.max.manual, hard.thresholds=hard.thresholds,
                       thresholds.col=thresholds.col)
  }
  dev.off()
}

#' Plot single chromosome windows of haplotype-based genome scan
#'
#' This function takes the genome scan output from scan.h2lmm() and plots the portion that corresponds to a single chromosome.
#' When multiple imputations are used, includes the 95/% confidence band on the median.
#'
#' @param scan.object A scan.h2lmm() object (ROP or multiple imputations). If multiple imputations, median and confidence interval 
#' on median are plotted.
#' @param chr The chromosome to be plotted.
#' @param use.lod DEFAULT: TRUE. Plots either the LOD score or the -log10 p-value.
#' @param scale DEFAULT: "Mb". Specifies the scale of genomic position to be plotted. Either Mb or cM can be used.
#' @param main.colors DEFAULT: "black". The color of the main association score to be plotted.
#' @param median.band.col DEFAULT: "gray88". The color of the 95% confident band plotted around the median.
#' @param main DEFAULT: "". Adds a title above the model.
#' @param y.max.manual DEFAULT: NULL. Manually adds a max y-value. Allows multiple genome scans to easily be on the same scale.
#' @param hard.thresholds DEFAULT: NULL. Specify one or more horizontal threshold lines.
#' @param thresholds.col DEFAULT: "red". Set the colors of the specified thresholds.
#' @export
#' @examples genome.plotter.chr()
genome.plotter.chr <- function(scan.object, chr, use.lod=TRUE,
                               scale=c("Mb", "cM"), main.col="black", median.band.col="gray88",
                               main="",
                               y.max.manual=NULL,
                               hard.thresholds=NULL, thresholds.col="red"){
  scale <- scale[1]
  MI <- all.CI <- CI <- NULL
  if(length(thresholds.col) < length(hard.thresholds)){ rep(thresholds.col, length(hard.thresholds)) }
  
  if(use.lod){
    all.outcome <- scan.object$LOD
    outcome <- all.outcome[scan.object$chr == chr]
    plot.this <- "LOD"
    this.ylab <- "LOD"
    if(!is.null(scan.object$MI.LOD)){
      all.MI <- scan.object$MI.LOD
      # Finding the 95% CI on the median
      all.CI <- apply(all.MI, 2, function(x) ci.median(x))
      CI <- all.CI[,scan.object$chr == chr]
    }
  }
  else{
    all.outcome <- -log10(scan.object$p.value)
    outcome <- all.outcome[scan.object$chr == chr]
    plot.this <- "p.value"
    this.ylab <- expression("-log"[10]*"P")
    if(!is.null(scan.object$MI.LOD)){
      all.MI <- -log10(scan.object$MI.p.value)
      # Finding the 95% CI on the median
      all.CI <- apply(all.MI, 2, function(x) ci.median(x, conf=0.95))
      CI <- all.CI[,scan.object$chr == chr]
    }
  }
  
  if(scale == "Mb"){
    pos <- scan.object$pos$Mb[scan.object$chr == chr]
  }
  else if(scale == "c<"){
    pos <- scan.object$pos$cM[scan.object$chr == chr]
  }
  
  order.i <- order(pos)
  
  outcome <- outcome[order.i]
  pos <- pos[order.i]
  
  min.pos <- min(pos)
  max.pos <- max(pos)

  
  # Finding max y of plot window
  y.max <- ceiling(max(all.outcome, hard.thresholds, all.CI[2,])) 
  if(!is.null(y.max.manual)){
    y.max <- y.max.manual
  }
  
  this.title <- c(main, paste0(scan.object$formula, " + locus (", scan.object$model.type, ")"))
  
  plot(pos, outcome, 
       xlim=c(0, max.pos), 
       ylim=c(0, y.max), 
       yaxt="n", xlab=paste("Chr", chr, paste0("(", scale, ")")), ylab=this.ylab, main=this.title,
       frame.plot=F, type="l", lwd=1.5, col=main.col)
  axis(side=2, at=0:y.max, las=2)
  if(!is.null(CI)){
    polygon(x=c(pos, rev(pos)), y=c(CI[1,], rev(CI[2,])), density=NA, col=median.band.col)
  }
  points(pos, outcome, type="l", pch=20, cex=0.5, lwd=1.5, col=main.col)

  if(!is.null(hard.thresholds)){
    for(i in 1:length(hard.thresholds)){
      abline(h=hard.thresholds[i], col=thresholds.col[i], lty=2)
    }
  }
}

#' Plot one or more haplotype-based genome scans flexibly (whole genomes or subset of chromosomes)
#'
#' This function takes the genome scan output from scan.h2lmm() and flexibly plots out the genome scan.
#'
#' @param non.mi.scan.list A list of scan.h2lmm() objects that are plotted within the same genome scan plot.
#' @param mi.scan DEFAULT: NULL. Multiple imputation scan.h2lmm() object.
#' @param use.lod DEFAULT: TRUE. Plots either the LOD score or the -log10 p-value.
#' @param scale DEFAULT: "Mb". Specifies the scale of genomic position to be plotted. Either Mb or cM are expected.
#' @param main.colors DEFAULT: "black". The color of the main association score to be plotted.
#' @param use.legend DEFAULT: TRUE. Include a legend for the different associations. If TRUE, the labels are the names of the non.mi.scan.list object.
#' @param main DEFAULT: NULL. Adds a title above the model.
#' @param my.legend.cex DEFAULT: 0.6. Specifies the size of the text in the legend.
#' @param y.max.manual DEFAULT: NULL. Manually adds a max y-value. Allows multiple genome scans to easily be on the same scale.
#' @param hard.thresholds DEFAULT: NULL. Specify one or more horizontal threshold lines.
#' @param thresholds.col DEFAULT: "red". Set the colors of the specified thresholds.
#' @export
#' @examples genome.plotter.whole()
genome.plotter.whole <- function(non.mi.scan.list, mi.scan=NULL, use.lod=TRUE, just.these.chr=NULL,
                                 scale="Mb", main.colors=c("black", "gray48", "blue"), mi.colors=c("purple", "cyan", "dark green"),
                                 use.legend=TRUE, main="",
                                 my.legend.cex=0.6,
                                 bs.max=NULL,
                                 y.max.manual=NULL,
                                 hard.thresholds=NULL, thresholds.col="red")
{
  if(length(thresholds.col) < length(hard.thresholds)){ rep(thresholds.col, length(hard.thresholds)) }
  main.object <- non.mi.scan.list[[1]]
  if(use.lod){
    outcome <- main.object$LOD
    plot.this <- "LOD"
    this.ylab <- "LOD"
  }
  if(!use.lod){
    outcome <- -log10(main.object$p.value)
    plot.this <- "p.value"
    this.ylab <- expression("-log"[10]*"P")
  }
  chr <- main.object$chr
  
  pos <- ifelse(rep(scale=="Mb", length(outcome)), main.object$pos$Mb, main.object$pos$cM)
  
  if(!is.null(just.these.chr)){
    keep.chr <- chr %in% just.these.chr
    chr <- chr[keep.chr]
    outcome <- outcome[keep.chr]
    pos <- pos[keep.chr]
  }
  
  has.X <- FALSE
  if(any(chr=="X")){
    has.X <- TRUE
    chr[chr=="X"] <- length(unique(chr))
  }
  
  pre.chr <- as.factor(as.numeric(chr))
  order.i <- order(pre.chr, pos)
  
  outcome <- outcome[order.i]
  pre.chr <- pre.chr[order.i]
  pos <- pos[order.i]
  
  min.pos <- tapply(pos, pre.chr, function(x) min(x))
  max.pos <- tapply(pos, pre.chr, function(x) max(x))
  chr.types <- levels(pre.chr)
  
  # Finding thresholds
  bs.threshold.95 <- bs.threshold.90 <- NULL
  if(!is.null(bs.max[[1]])){
    bs.threshold.95 <- get.gev.thresholds(bs.max[[1]])
    bs.threshold.90 <- get.gev.thresholds(bs.max[[1]], percentile=0.9)
  }
  
  mi.bs.threshold.95 <- mi.bs.threshold.90 <- NULL
  if(!is.null(bs.max[[2]])){
    mi.bs.threshold.95 <- get.gev.thresholds(bs.max[[2]])
    mi.bs.threshold.90 <- get.gev.thresholds(bs.max[[2]], percentile=0.9)
  }
  
  # Setting up multiple imputations
  if(is.null(mi.scan)){
    mi.mat <- NULL
  }
  if(!is.null(mi.scan)){
    if(use.lod){
      mi.mat <- mi.scan$mi.LOD[,order.i]
    }
    if(!use.lod){
      mi.mat <- -log10(mi.scan$mi.p.mat)[,order.i]
    }
  }
  
  # Finding max y of plot window
  y.max <- ceiling(max(outcome, bs.threshold.95, mi.bs.threshold.95, mi.mat, hard.thresholds)) 
  if(length(non.mi.scan.list) > 1){
    for(i in 2:length(non.mi.scan.list)){
      if(use.lod){
        y.max <- ceiling(max(y.max, unlist(non.mi.scan.list[[i]][plot.this])))
      }
      if(!use.lod){
        y.max <- ceiling(max(y.max, -log10(unlist(non.mi.scan.list[[i]][plot.this]))))
      }
    }
  }
  if(!is.null(y.max.manual)){
    y.max <- y.max.manual
  }
  
  shift.left <- min(pos[chr==chr.types[1]])

  this.title <- c(main, paste0(non.mi.scan.list[[1]]$formula, " + locus (", non.mi.scan.list[[1]]$model.type, ")"))
  
  plot(pos[pre.chr==chr.types[1]], outcome[pre.chr==chr.types[1]], 
       xlim=c(shift.left, sum(max.pos)+(length(chr.types)-1)), 
       ylim=c(-0.1, y.max), 
       xaxt="n", yaxt="n", xlab="", ylab=this.ylab, main=this.title,
       frame.plot=F, type="l", pch=20, cex=0.5, lwd=1.5, col=main.colors[1])
  axis(side=2, at=0:y.max, las=2)
  if(!is.null(mi.scan)){
    for(i in 1:nrow(mi.mat)){
      points(pos[pre.chr==chr.types[1]], mi.mat[i,][pre.chr==chr.types[1]], type="l", lwd=0.5, col="orange")
    }
    points(pos[pre.chr==chr.types[1]], outcome[pre.chr==chr.types[1]], type="l", lwd=1.5, col=main.colors[1])
    if(use.lod){
      mi.outcome.LOD.scale <- unlist(mi.scan["LOD"])[order.i]
      points(pos[pre.chr==chr.types[1]], mi.outcome.LOD.scale[pre.chr==chr.types[1]], type="l", lwd=0.75, col=scales::alpha(mi.colors[1], 0.5))
      points(pos[pre.chr==chr.types[1]], mi.outcome.LOD.exp.scale[pre.chr==chr.types[1]], type="l", lwd=0.75, col=scales::alpha(mi.colors[2], 0.5))
    }
    if(!use.lod){
      mi.outcome.pvalue.scale <- unlist(mi.scan["p.value"])[order.i]
      points(pos[pre.chr==chr.types[1]], mi.outcome.pvalue.scale[pre.chr==chr.types[1]], type="l", lwd=0.75, col=scales::alpha(mi.colors[1], 0.5))
    }
  }
  
  label.spots <- min.pos[1] + (max.pos[1] - min.pos[1])/2
  shift <- max.pos[1]
  if(length(chr.types) > 1){
    for(i in 2:length(chr.types)){
      this.pos <- pos[pre.chr==chr.types[i]] + shift
      if(i %% 2 == 0){
        polygon(x=c(min(this.pos), min(this.pos):max(this.pos), max(this.pos)), 
                y=c(0, rep(y.max, length(min(this.pos):max(this.pos))), 0), border=NA, col="gray88")
      }
      label.spots <- c(label.spots, min.pos[i] + shift + (max.pos[i] - min.pos[i])/2)
      points(this.pos, outcome[pre.chr==chr.types[i]], type="l", lwd=1.5, col=main.colors[1])
      if(!is.null(mi.scan)){
        for(j in 1:nrow(mi.mat)){
          points(pos[pre.chr==chr.types[i]] + shift, mi.mat[j,][pre.chr==chr.types[i]], type="l", lwd=0.5, col="orange")
        }
        points(pos[pre.chr==chr.types[i]] + shift, outcome[pre.chr==chr.types[i]], type="l", lwd=1.5, col=main.colors[1])
        if(use.lod){
          points(pos[pre.chr==chr.types[i]] + shift, mi.outcome.LOD.scale[pre.chr==chr.types[i]], type="l", lwd=0.75, col=scales::alpha(mi.colors[1], 0.5))
          points(pos[pre.chr==chr.types[i]] + shift, mi.outcome.LOD.exp.scale[pre.chr==chr.types[i]], type="l", lwd=0.75, col=scales::alpha(mi.colors[2], 0.5))
        }
        if(!use.lod){
          points(pos[pre.chr==chr.types[i]] + shift, mi.outcome.pvalue.LOD.scale[pre.chr==chr.types[i]], type="l", lwd=0.75, col=scales::alpha(mi.colors[1], 0.5))
          points(pos[pre.chr==chr.types[i]] + shift, mi.outcome.pvalue.exp.scale[pre.chr==chr.types[i]], type="l", lwd=0.75, col=scales::alpha(mi.colors[2], 0.5))
          points(pos[pre.chr==chr.types[i]] + shift, mi.outcome.mi.pvalue[pre.chr==chr.types[i]], type="l", lwd=0.75, col=scales::alpha(mi.colors[3], 0.5))
        }
      }
      shift <- shift + max.pos[i]
    }
  }
  
  # Plot other method's statistics
  if(length(non.mi.scan.list) > 1){
    for(i in 2:length(non.mi.scan.list)){
      if(use.lod){
        compare.shift <- 0
        compare.outcome <- unlist(non.mi.scan.list[[i]]["LOD"])[keep.chr][order.i]
        for(j in 1:length(chr.types)){
          points(pos[pre.chr==chr.types[j]] + compare.shift, compare.outcome[pre.chr==chr.types[j]], type="l", col=main.colors[i], lwd=1.5)
          compare.shift <- compare.shift + max.pos[j]
        }
      }
      if(!use.lod){
        compare.shift <- 0
        compare.outcome <- -log10(unlist(non.mi.scan.list[[i]]["p.value"]))[keep.chr][order.i]
        for(j in 1:length(chr.types)){
          points(pos[pre.chr==chr.types[j]] + compare.shift, compare.outcome[pre.chr==chr.types[j]], type="l", col=main.colors[i], lwd=1.5)
          compare.shift <- compare.shift + max.pos[j]
        }
      }
    }
  }
  # Permutation threshold
  if(!is.null(bs.max[[1]])){
    abline(h=bs.threshold.95, lty=2, col="red")
    abline(h=bs.threshold.90, lty=2, col="blue")
    if(!is.null(bs.max[[2]])){
      abline(h=mi.bs.threshold.95, lty=6, col="red")
      abline(h=mi.bs.threshold.90, lty=6, col="blue")
    }
  }
  if(has.X){
    axis.label <- c(chr.types[-length(chr.types)], "X")
  }
  if(!has.X){
    axis.label <- chr.types
  }
  axis(side=1, tick=F, line=NA, at=label.spots, labels=axis.label, cex.axis=0.7, padj=-1.5)
  if(use.legend){
    if(is.null(mi.scan)){
      legend("topright", legend=names(non.mi.scan.list), 
             lty=rep(1, length(non.mi.scan.list)), lwd=c(rep(1.5, length(non.mi.scan.list))), 
             col=main.colors[1:length(non.mi.scan.list)], bty="n", cex=my.legend.cex)
    }
    if(!is.null(mi.scan)){
      if(use.lod){
        legend("topright", legend=c(names(non.mi.scan.list), paste0("MI - ", nrow(mi.mat)), "MI: Mean LOD", "MI: Mean LOD (Lik scale)"), 
               lty=rep(1, length(non.mi.scan.list)+3), lwd=c(1.5, rep(1.5, length(non.mi.scan.list)+2)), 
               col=c(main.colors[1:length(non.mi.scan.list)], "orange", mi.colors[1:2]), bty="n", cex=my.legend.cex)
      }
      if(!use.lod){
        legend("topright", legend=c(names(non.mi.scan.list), paste0("MI - ", nrow(mi.mat)), "MI: Mean P-value (LOD scale)", "MI: Mean P-value (P scale)", "MI: Rubin MI LRT P-value"), 
               lty=rep(1, length(non.mi.scan.list)+4), lwd=c(1.5, rep(1.5, length(non.mi.scan.list)+3)), 
               col=c(main.colors[1:length(non.mi.scan.list)], "orange", mi.colors[1:3]), bty="n", cex=my.legend.cex)
      }
    }
    if(!is.null(bs.max[[1]]) & is.null(bs.max[[2]])){
      legend("topleft", legend=c("ROP 95% Bootstrap Thresholds - 100 samples", "ROP 90% Bootstrap Thresholds - 100 samples"), 
             lty=c(2, 2), lwd=c(1, 1), 
             col=c("red", "blue"), bty="n", cex=my.legend.cex)
    }
    if(!is.null(bs.max[[1]]) & !is.null(bs.max[[2]])){
      legend("topleft", 
             legend=c("ROP 95% Bootstrap Thresholds - 100 samples", "ROP 90% Bootstrap Thresholds - 100 samples", 
                      paste0("MI", nrow(mi.mat), " 95% Bootstrap Thresholds - 100 samples"), paste0("MI", nrow(mi.mat), " 90% Bootstrap Thresholds - 100 samples")), 
             lty=c(2, 2, 6, 6), lwd=c(1, 1, 1, 1), 
             col=c("red", "blue", "red", "blue"), bty="n", cex=my.legend.cex)
    }
  }
  if(!is.null(hard.thresholds)){
    for(i in 1:length(hard.thresholds)){
      abline(h=hard.thresholds[i], col=thresholds.col[i], lty=2)
    }
  }
}

#' @export
snp.genome.plotter.whole <- function(snp.scan, use.lod=TRUE, just.these.chr=NULL,
                                     scale="Mb",
                                     main.label=NULL,
                                     bs.max=NULL,
                                     cache.title="DO-QTL",
                                     approach.title="Exact", y.max.manual=NULL, title="", alt.col=NULL,
                                     hard.thresholds=NULL)
{
  main.object <- snp.scan
  if(use.lod){
    outcome <- main.object$LOD
    plot.this <- "LOD"
    this.ylab <- "LOD"
  }
  if(!use.lod){
    outcome <- -log10(main.object$p.value)
    plot.this <- "p.value"
    this.ylab <- expression("-log"[10]*"P")
  }
  
  # Allowing for special colors
  if(is.null(alt.col)){ use.col <- rep("black", length(outcome)) }
  if(!is.null(alt.col)){ use.col <- alt.col }
  chr <- main.object$chr
  has.X <- FALSE
  if(any(chr=="X")){
    has.X <- TRUE
    chr[chr=="X"] <- length(unique(chr))
  }
  pos <- ifelse(rep(scale=="Mb", length(outcome)), main.object$pos$Mb, main.object$pos$cM)
  
  if(!is.null(just.these.chr)){
    keep.chr <- chr %in% just.these.chr
    chr <- chr[keep.chr]
    outcome <- outcome[keep.chr]
    pos <- pos[keep.chr]
  }
  
  has.X <- FALSE
  if(any(chr=="X")){
    has.X <- TRUE
    chr[chr=="X"] <- length(unique(chr))
  }
  
  pre.chr <- as.factor(as.numeric(chr))
  order.i <- order(pre.chr, pos)
  
  outcome <- outcome[order.i]
  pre.chr <- pre.chr[order.i]
  pos <- pos[order.i]
  use.col <- use.col[order.i]
  
  min.pos <- tapply(pos, pre.chr, function(x) min(x))
  max.pos <- tapply(pos, pre.chr, function(x) max(x))
  chr.types <- levels(pre.chr)
  
  # Finding thresholds
  bs.threshold.95 <- bs.threshold.90 <- NULL
  if(!is.null(bs.max[[1]])){
    bs.threshold.95 <- get.gev.thresholds(bs.max[[1]])
    bs.threshold.90 <- get.gev.thresholds(bs.max[[1]], percentile=0.9)
  }
  
  mi.bs.threshold.95 <- mi.bs.threshold.90 <- NULL
  if(!is.null(bs.max[[2]])){
    mi.bs.threshold.95 <- get.gev.thresholds(bs.max[[2]])
    mi.bs.threshold.90 <- get.gev.thresholds(bs.max[[2]], percentile=0.9)
  }
  
  # Finding max y of plot window
  y.max <- ceiling(max(outcome, bs.threshold.95, mi.bs.threshold.95, hard.thresholds)) 
  if(!is.null(y.max.manual)){
    y.max <- y.max.manual
  }
  
  shift.left <- min(pos[chr==chr.types[1]])
  this.title <- c(title,
                  paste0(main.object$formula, " + SNP (", main.object$model.type, ")"))

  plot(1,
       xlim=c(shift.left, sum(max.pos)+(length(chr.types)-1)),
       ylim=c(-1, y.max),
       xaxt="n", yaxt="n", xlab="", ylab=this.ylab, main=this.title,
       frame.plot=F, type="n")
  axis(side=2, at=0:y.max, las=2)
  
  label.spots <- min.pos[1] + (max.pos[1] - min.pos[1])/2
  shift <- max.pos[1]
  if(length(chr.types) > 1){
    for(i in 2:length(chr.types)){
      this.pos <- pos[pre.chr==chr.types[i]] + shift
      if(i %% 2 == 0){
        polygon(x=c(min(this.pos), min(this.pos):max(this.pos), max(this.pos)), 
                y=c(0, rep(y.max, length(min(this.pos):max(this.pos))), 0), border=NA, col="gray88")
      }
      label.spots <- c(label.spots, min.pos[i] + shift + (max.pos[i] - min.pos[i])/2)

      shift <- shift + max.pos[i]
    }
  }
  # Adding in the points
  shift <- max.pos[1]
  
  if(length(chr.types) >= 1){
    points(pos[pre.chr==chr.types[1]], outcome[pre.chr==chr.types[1]], pch=20, cex=0.5, col=use.col[pre.chr==chr.types[1]])
    if(length(chr.types) > 1){
      for(i in 2:length(chr.types)){
        this.pos <- pos[pre.chr==chr.types[i]] + shift
        points(this.pos, outcome[pre.chr==chr.types[i]], type="p", pch=20, cex=0.5, col=use.col[pre.chr==chr.types[i]])
        
        shift <- shift + max.pos[i]
      }
    }
  }
  
  
  # Permutation threshold
  if(!is.null(bs.max[[1]])){
    abline(h=bs.threshold.95, lty=2, col="red")
    abline(h=bs.threshold.90, lty=2, col="blue")
    if(!is.null(bs.max[[2]])){
      abline(h=mi.bs.threshold.95, lty=6, col="red")
      abline(h=mi.bs.threshold.90, lty=6, col="blue")
    }
  }
  if(has.X){
    axis.label <- c(chr.types[-length(chr.types)], "X")
  }
  if(!has.X){
    axis.label <- chr.types
  }
  axis(side=1, tick=F, line=NA, at=label.spots, labels=axis.label, cex.axis=0.7, padj=-3.5)
}

#' @export
snp.genome.plotter.w.r2 <- function(snp.scan, r2.object, use.lod=TRUE,
                                    scale="Mb",
                                    main.label=NULL,
                                    bs.max=NULL,
                                    cache.title="DO-QTL",
                                    approach.title="Exact", y.max.manual=NULL, title="", alt.col=NULL, this.cex=1,
                                    hard.thresholds=NULL)
{
  main.object <- snp.scan
  if(use.lod){
    outcome <- main.object$LOD
    plot.this <- "LOD"
    this.ylab <- "LOD"
  }
  if(!use.lod){
    outcome <- -log10(main.object$p.value)
    plot.this <- "p.value"
    this.ylab <- expression("-log"[10]*"P")
  }
  
  # Allowing for special colors
  if(is.null(alt.col)){ use.col <- rep("black", length(outcome)) }
  if(!is.null(alt.col)){ use.col <- alt.col }
  
  pos <- ifelse(rep(scale=="Mb", length(outcome)), main.object$pos$Mb, main.object$pos$cM)
  
  point.locus <- r2.object$point.locus
  point.locus.outcome <- outcome[point.locus == main.object$loci]
  point.locus.pos <- pos[point.locus == main.object$loci]
  
  chr <- r2.object$chr
  outcome <- outcome[main.object$chr == chr]
  pos <- pos[main.object$chr == chr]
  
  order.i <- order(pos)
  
  outcome <- outcome[order.i]
  pos <- pos[order.i]
  use.col <- use.col[order.i]
  
  min.pos <- min(pos)
  max.pos <- max(pos)
  
  # Finding thresholds
  bs.threshold.95 <- bs.threshold.90 <- NULL
  if(!is.null(bs.max[[1]])){
    bs.threshold.95 <- get.gev.thresholds(bs.max[[1]])
    bs.threshold.90 <- get.gev.thresholds(bs.max[[1]], percentile=0.9)
  }
  
  mi.bs.threshold.95 <- mi.bs.threshold.90 <- NULL
  if(!is.null(bs.max[[2]])){
    mi.bs.threshold.95 <- get.gev.thresholds(bs.max[[2]])
    mi.bs.threshold.90 <- get.gev.thresholds(bs.max[[2]], percentile=0.9)
  }
  
  # Finding max y of plot window
  y.max <- ceiling(max(outcome, bs.threshold.95, mi.bs.threshold.95, hard.thresholds)) 
  if(!is.null(y.max.manual)){
    y.max <- y.max.manual
  }
  
  this.title <- c(title,
                  paste0(main.object$formula, " + SNP (", main.object$model.type, ")"))
  this.xlab <- paste0("Chr ", chr, " Position (", scale, ")")
  
  red2blue <- colorRampPalette(c("red", "blue"))
  these.colors <- rev(red2blue(1000))
  
  r2.col <- these.colors[ceiling(r2.object$r2*999.1)]
  
  plot(x=pos, y=outcome,
       xlim=c(min.pos, max.pos),
       ylim=c(0, y.max+1),
       pch=20, cex=this.cex, col=r2.col,
       yaxt="n", xlab=this.xlab, ylab=this.ylab, main=this.title,
       frame.plot=F)
  points(x=point.locus.pos, y=point.locus.outcome, 
         bg="red", pch=21, cex=1.5)
  axis(side=2, at=0:y.max, las=2)

  plotrix::color.legend(xl=floor(0.75*max.pos), yb=y.max, xr=max.pos, yt=y.max+0.5, legend=c(0, 0.25, 0.5, 0.75, 1), rect.col=these.colors, align="rb", gradient="x")  
  text(x=(max.pos - floor(0.75*max.pos))/2 + floor(0.75*max.pos),
       y=y.max+0.75,
       labels="r2 with peak SNP")
  
  # Permutation threshold
  if(!is.null(bs.max[[1]])){
    abline(h=bs.threshold.95, lty=2, col="red")
    abline(h=bs.threshold.90, lty=2, col="blue")
    if(!is.null(bs.max[[2]])){
      abline(h=mi.bs.threshold.95, lty=6, col="red")
      abline(h=mi.bs.threshold.90, lty=6, col="blue")
    }
  }
}

#' @export
single.chr.plotter.w.ci <- function(scan.object, qtl.ci.object, 
                                    ci.type, scan.type, 
                                    qtl.only=NULL, subsample.only=NULL, supress.ci=FALSE,
                                    scale="cM", these.col=c("#7BAFD4", "red"), manual.interval=NULL, manual.scale=NULL){
  
  outcome <- -log10(scan.object$p.value) 
  this.ylab <- expression("-log"[10]*"P")
  
  pos <- unlist(scan.object$pos[scale])
  
  order.i <- order(pos)
  pos <- pos[order.i]
  outcome <- outcome[order.i]
  
  all.loci <- dimnames(qtl.ci.object$full.results)[[3]]
  
  peak.pos <- qtl.ci.object$peak.pos[[scale]]
  loci <- qtl.ci.object$peak.loci
  if(!is.null(qtl.only)){
    peak.pos <- qtl.ci.object$peak.pos[[scale]][qtl.only,]
    loci <- qtl.ci.object$peak.loci[,qtl.only]
  }
  
  # Process per CI
  ci.names <- names(qtl.ci.object$ci)
  for(i in 1:length(qtl.ci.object$ci)){
    ci <- qtl.ci.object$ci[[ci.names[i]]][[scale]]
    
    if(!is.null(qtl.only)){ 
      ci <-  ci[qtl.only,] 
    }
    lb.dist <- pos - ci[1]
    low.locus <- all.loci[lb.dist <= 0][which.max(lb.dist[lb.dist <= 0])]
    low.locus.pos <- scan.object$pos[[scale]][which(all.loci == low.locus)]
    ub.dist <- pos - ci[2]
    high.locus <- all.loci[ub.dist >= 0][which.min(ub.dist[ub.dist >= 0])]
    high.locus.pos <- scan.object$pos[[scale]][which(all.loci == high.locus)]
    
    region <- scan.object$pos[[scale]] >= low.locus.pos & scan.object$pos[[scale]] <= high.locus.pos
    peak.locus <- all.loci[region][which.min(scan.object$p.value[region])]
    peak.locus.pos <- scan.object$pos[[scale]][region][which.min(scan.object$p.value[region])]
  }
  
  if(!is.null(manual.interval)){
    main.title <- c(paste0(scan.type, ": ", scan.object$formula, " + locus (", scan.object$model.type, ")"),
                    paste0("QTL interval type: ", ci.type),
                    paste0("Width: ", round(manual.interval[2] - manual.interval[1], 2), manual.scale),
                    paste0("peak locus: ", peak.locus, " (", round(manual.interval[3], 4), manual.scale, ")"),
                    paste0("(closest) lower locus: ", low.locus, " (", round(manual.interval[1], 4), manual.scale, ")"),
                    paste0("(closest) upper locus: ", high.locus, " (", round(manual.interval[2], 4), manual.scale, ")"))
  }
  if(!supress.ci & is.null(manual.interval)){
    main.title <- c(paste0(scan.type, ": ", scan.object$formula, " + locus (", scan.object$model.type, ")"),
                    paste0("QTL interval type: ", ci.type),
                    paste0("Width: ", round(high.locus.pos - low.locus.pos, 2), scale),
                    paste0("peak locus: ", peak.locus, " (", round(peak.locus.pos, 4), scale, ")"),
                    paste0("(closest) lower locus: ", low.locus, " (", round(low.locus.pos, 4), scale, ")"),
                    paste0("(closest) upper locus: ", high.locus, " (", round(high.locus.pos, 4), scale, ")"))
  }
  if(supress.ci){
    main.title <- c(paste0(scan.type, ": ", scan.object$formula, " + locus (", scan.object$model.type, ")"),
                    paste0("QTL interval type: ", ci.type))
  }
  #browser()
  this.xlab <- paste("Chr", qtl.ci.object$chr, paste0("(", scale, ")"))
  y.max <- max(outcome, -log10(qtl.ci.object$full.results))
  plot(1, 
       xlim=c(0, max(pos)), 
       ylim=c(0, y.max), 
       xlab=this.xlab, ylab=this.ylab, main=main.title,
       frame.plot=F, type="l", pch=20, cex=0.5, las=1, cex.main=0.8)
  if(!supress.ci){
    polygon(c(rep(low.locus.pos, 2), rep(high.locus.pos, 2)), c(0, rep(ceiling(max(outcome)), 2), 0), col="gray", border=NA)
  }
  
  full.results <- qtl.ci.object$full.results
  peaks <- qtl.ci.object$peak.pos[[scale]]
  if(!is.null(qtl.only)){
    full.results <- full.results[,qtl.only,,drop=FALSE]
    peaks <- peaks[qtl.only,,drop=FALSE]
  }
  if(!is.null(subsample.only)){
    full.results <- full.results[subsample.only,,,drop=FALSE]
    peaks <- peaks[,subsample.only,drop=FALSE]
  }
  for(i in 1:dim(full.results)[2]){
    for(j in 1:dim(full.results)[1]){
      lines(pos, -log10(full.results[j,i,]), lwd=0.5, col=scales::alpha(these.col[i], 0.5))
    }
  }
  rug(qtl.ci.object$peak.pos[[scale]], col=scales::alpha("black", 0.5))
  rug(qtl.ci.object$peak.pos[[scale]][,qtl.only], col=scales::alpha(these.col[i], 0.75))
  
  lines(pos, outcome, lwd=1.5)
}

inspect.ci.genome.plotter.whole <- function(ci.object, scan.type.label, which.ci=1, ...){
  this.scan <- list(p.value = ci.object$full.results[which.ci,],
                    chr=rep(ci.object$chr, length(ci.object$full.results[which.ci,])),
                    pos = ci.object$pos)
  this.scan.list <- list()
  this.scan.list[[scan.type.label]] <- this.scan
  genome.plotter.whole(non.mi.scan.list=this.scan.list, use.lod=FALSE, scale="cM", ...)
}


