#' Takes DO-QTL founder haplotype reconstruction output files and re-formats into a HAPPY genome cache
#'
#' This function produces a HAPPY-format genome cache from DO-QTL founder haplotype reconstruction output files.
#' The main output files of importance are the individual-level files with naming scheme [sample].genotype.probs.Rdata.
#'
#' @param DOQTL.recon.output.path The path to the directory containing DO-QTL founder haplotype output files. 
#' @param map.path The path to the map file. The map file contains data on the loci. It should be a tab-delimited 
#' file with columns labeled "marker", "chr", "bp", and "pos". "pos" represents map distance in cM. "bp" should 
#' be in bp, not Mb. 
#' @param HAPPY.output.path The path to a directory that will be created as the HAPPY-format genome cache.
#' @param allele.labels DEFAULT: NULL. Allows for specification of founder labels different from what is in the DO-QTL
#' output. The DEFAULT of NULL leads to using the labels from the DO-QTL output.
#' @param chr DEFAULT: c(1:19, "X"). Allows for specification of the chromosomes. DEFAULT is all the chromosomes from the mouse.
#' @export
#' @examples convert.DOQTL.to.HAPPY()
convert.DOQTL.to.HAPPY <- function(DOQTL.recon.output.path,
                                   map.path,
                                   HAPPY.output.path,
                                   allele.labels=NULL,
                                   chr=c(1:19, "X")){
  
  require(data.table)
  #----------------------------------
  # founder probs from DO-QTL
  #----------------------------------
  load(paste(DOQTL.recon.output.path, "founder.probs.Rdata", sep="/"))
  
  samples <- dimnames(model.probs)[[1]]
  simple.alleles <- dimnames(model.probs)[[2]]
  rm(model.probs)
  
  
  #----------------------------------
  # Putting together strain labels
  #----------------------------------
  if(is.null(allele.labels)){
    allele.labels <- simple.alleles
  }

  full.to.dosages <- straineff.mapping.matrix()
  
  diplotype.labels <- c(paste(allele.labels, allele.labels, sep="."), apply(full.to.dosages[-(1:8),], 1, 
                                                                            function(x) paste(allele.labels[sort(which(x==1), decreasing=FALSE)], collapse=".")))
  simple.het.labels <- c(paste(simple.alleles, simple.alleles, sep=""), apply(full.to.dosages[-(1:8),], 1, 
                                                                              function(x) paste(simple.alleles[sort(which(x==1), decreasing=FALSE)], collapse="")))
  
  
  #-------------------------------
  # Marker info
  #-------------------------------
  map <- read.table(map.path, header=TRUE, as.is=TRUE)
  
  
  #-------------------------------
  # Functions to output marker files
  #-------------------------------
  f <- function(chr, Marker.Name,
                aa, bb, cc, dd, ee, ff, gg, hh,
                ba, ca, cb, da, db, dc, ea, eb, ec, ed,
                fa, fb, fc, fd, fe, ga, gb, gc, gd, ge, gf,
                ha, hb, hc, hd, he, hf, hg,
                allele.labels,
                diplotype.labels,
                full.to.dosages.matrix){
    var_name <- Marker.Name[1]
    assign(var_name, matrix(data=c(aa, bb, cc, dd, ee, ff, gg, hh,
                                   ba, ca, cb, da, db, dc, ea, eb, ec, ed,
                                   fa, fb, fc, fd, fe, ga, gb, gc, gd, ge, gf,
                                   ha, hb, hc, hd, he, hf, hg), 
                            ncol=length(diplotype.labels),
                            dimnames=list(NULL, diplotype.labels)))
    
    temp <- get(var_name)
    colnames(temp) <- diplotype.labels
    temp.add <- temp %*% full.to.dosages.matrix
    colnames(temp.add) <- allele.labels
    
    dir.create(paste0(HAPPY.output.path, '/full/chr', chr[1], '/data/'),
               showWarnings=FALSE, recursive=TRUE)
    fn <- paste0(HAPPY.output.path, '/full/chr', chr[1], '/data/', var_name, '.RData')
    save(list=var_name, file=fn)
    
    assign(var_name, temp.add)
    dir.create(paste0(HAPPY.output.path, '/additive/chr', chr[1], '/data/'),
               showWarnings=FALSE, recursive=TRUE)
    fn <- paste0(HAPPY.output.path, '/additive/chr', chr[1], '/data/', var_name, '.RData')
    save(list = var_name, file = fn) 
  }
  # export file of marker names for each chr
  export_marker_name_file <- function(chr, Marker.Name) {
    markers <- as.character(Marker.Name)
    save(markers, file=paste0(HAPPY.output.path, '/additive/chr', chr[1], '/markers.RData'))
    save(markers, file=paste0(HAPPY.output.path, '/full/chr', chr[1], '/markers.RData'))
    dir.create(paste0(HAPPY.output.path, '/genotype/chr', chr[1]), showWarnings=FALSE, recursive=TRUE)
    save(markers, file=paste0(HAPPY.output.path, '/genotype/chr', chr[1], '/markers.RData'))
  }
  # export file of marker bp positions for each chr
  export_marker_position_file <- function(chr, pos) {
    bp <- as.character(pos)
    save(bp, file=paste0(HAPPY.output.path, '/additive/chr', chr[1], '/bp.RData'))
    save(bp, file=paste0(HAPPY.output.path, '/full/chr', chr[1], '/bp.RData'))
    save(bp, file=paste0(HAPPY.output.path, '/genotype/chr', chr[1], '/bp.RData'))
  }
  # export file of chromosome
  export_marker_chromosome_file <- function(chr) {
    chromosome <- as.character(chr)
    save(chromosome, file=paste0(HAPPY.output.path, '/additive/chr',chr[1], '/chromosome.RData'))
    save(chromosome, file=paste0(HAPPY.output.path, '/full/chr',chr[1], '/chromosome.RData'))
    save(chromosome, file=paste0(HAPPY.output.path, '/genotype/chr',chr[1], '/chromosome.RData'))
  }
  # export file of map distance (cM)
  export_marker_map_distance_file <- function(chr, pos) {
    map <- as.character(pos)
    save(map, file=paste0(HAPPY.output.path, '/additive/chr', chr[1], '/map.RData'))
    save(map, file=paste0(HAPPY.output.path, '/full/chr', chr[1], '/map.RData'))
    save(map, file=paste0(HAPPY.output.path, '/genotype/chr', chr[1], '/map.RData'))
  }
  
  #-----------------------------
  # combining data of individuals
  #-----------------------------
  for(i in 1:length(chr)){
    for(j in 1:length(samples[1:3])){
      cat(paste("Loading DOQTL output for individual", j, "for chr", chr[i]), "\n")
      load(paste0(DOQTL.recon.output.path, "/", samples[j], ".genotype.probs.Rdata"))
      marker <- rownames(prsmth)
      subject <- rep(samples[j], nrow(prsmth))
      one.sample.data <- data.frame(marker, subject,  prsmth)
      combined.data <- merge(x=map, y=one.sample.data, by.x="marker", by.y="marker")[,c(1:5,c(1,9,16,22,27,31,34,36,2,3,10,4,11,
                                                                                              17,5,12,18,23,6,13,19,24,28,7,14,
                                                                                              20,25,29,32,8,15,21,26,30,33,35)+5)]
      combined.data <- combined.data[combined.data$chr == chr[i],]
      #cat(paste(dim(combined.data), collapse=" "), "\n")
      #cat(paste(combined.data[1:10,]))
    
      if(!exists('all.subjects')){
        all.subjects <- combined.data
      } 
      else{ 
        all.subjects <- data.table::rbindlist(list(all.subjects, combined.data))
      }
    }
    #--------------------------------------------------------------------------
    # make each marker_name.Rdata
    # Subject order in marker_name.Rdata should match with SUBJECT.NAME in pheno
    #---------------------------------------------------------------------------
    var.names <- c(names(all.subjects)[1:5], simple.het.labels)
    data.table::setnames(all.subjects, names(all.subjects), var.names)
    data.table::setkey(all.subjects, NULL)
    data.table::setkey(all.subjects, chr, bp, marker, subject)
    
    # setnames(all.subjects, names(all.subjects), var.names)
    # setkey(all.subjects, NULL)
    # setkey(all.subjects, chr, bp, marker, subject)
    #all.subjects$marker.notkey <- all.subjects$marker
    #all.subjects$chr.notkey <- all.subjects$chr
    
    # all.subjects[, f(all.subjects$chr, all.subjects$marker,
    #                  all.subjects$AA, all.subjects$BB, all.subjects$CC, all.subjects$DD,
    #                  all.subjects$EE, all.subjects$FF, all.subjects$GG, all.subjects$HH,
    #                  all.subjects$AB, all.subjects$AC, all.subjects$BC, all.subjects$AD,
    #                  all.subjects$BD, all.subjects$CD, all.subjects$AE, all.subjects$BE,
    #                  all.subjects$CE, all.subjects$DE, all.subjects$AF, all.subjects$BF,
    #                  all.subjects$CF, all.subjects$DF, all.subjects$EF, all.subjects$AG,
    #                  all.subjects$BG, all.subjects$CG, all.subjects$DG, all.subjects$EG,
    #                  all.subjects$FG, all.subjects$AH, all.subjects$BH, all.subjects$CH,
    #                  all.subjects$DH, all.subjects$EH, all.subjects$FH, all.subjects$GH,
    #                  allele.labels, diplotype.labels, full.to.dosages), by=all.subjects$marker]
    all.subjects[, f(chr, marker,
                     AA, BB, CC, DD,
                     EE, FF, GG, HH,
                     AB, AC, BC, AD,
                     BD, CD, AE, BE,
                     CE, DE, AF, BF,
                     CF, DF, EF, AG,
                     BG, CG, DG, EG,
                     FG, AH, BH, CH,
                     DH, EH, FH, GH,
                     allele.labels, diplotype.labels, full.to.dosages), by="marker"]
    #cat(colnames(all.subjects), "\n")
    #cat(class(all.subjects), "\n")
    
    #cat(all.subjects[, AA, by="marker", with=TRUE])
    
    #-------------------------------------
    # make other necessary files
    #-------------------------------------
    one.subj = samples[1]
    markers.one.subj <- all.subjects[grepl(one.subj, subject), ]  
    
    markers.one.subj[, export_marker_name_file(all.subjects$chr.notkey, all.subjects$marker), by=all.subjects$chr]
    
    markers.one.subj[, export_marker_position_file(all.subjects$chr.notkey, all.subjects$bp), by=all.subjects$chr]
    
    markers.one.subj[, export_marker_chromosome_file(all.subjects$chr.notkey), by=all.subjects$chr]
    
    markers.one.subj[, export_marker_map_distance_file(all.subjects$chr.notkey, all.subjects$pos), by=all.subjects$chr]
    
    rm(all.subjects)
  }
  
  # export file of mice names and strains
  subjects <- samples
  strains <- allele.labels
  for(this.chr in chr){
    save(subjects, file = paste0(HAPPY.output.path, '/additive/chr', this.chr, '/subjects.RData'))
    save(subjects, file = paste0(HAPPY.output.path, '/full/chr', this.chr, '/subjects.RData'))
    save(subjects, file = paste0(HAPPY.output.path, '/genotype/chr', this.chr, '/subjects.RData'))
    
    save(strains, file = paste0(HAPPY.output.path, '/additive/chr', this.chr, '/strains.RData'))
    save(strains, file = paste0(HAPPY.output.path, '/full/chr', this.chr, '/strains.RData'))
    save(strains, file = paste0(HAPPY.output.path, '/genotype/chr', this.chr, '/strains.RData'))
  }
}






