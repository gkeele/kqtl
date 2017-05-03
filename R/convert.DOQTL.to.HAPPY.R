#' @export
convert.DOQTL.to.HAPPY <- function(DOQTL.recon.output.path,
                                   map.path,
                                   HAPPY.output.path,
                                   allele.labels=NULL,
                                   chr=c(1:20, "X")){
  
  #----------------------------------
  # founder probs from DO-QTL
  #----------------------------------
  load(paste(DOQTL.recon.output.path, "founder.probs.Rdata", sep="/"))
  
  samples <- dimnames(model.probs)[[1]]
  loci <- dimnames(model.probs)[[3]]
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
  map = read.table(map.path, header=TRUE, as.is=TRUE)
  
  
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
    for(j in 1:length(samples)){
      cat(paste("Loading DOQTL output for individual", j, "for chr", i), "\n")
      load(paste0(DOQTL.recon.output.path, "/", samples[i], ".genotype.probs.Rdata"))
      marker <- rownames(prsmth)
      subject <- rep(samples[i], nrow(prsmth))
      one.sample.data <- data.frame(marker, subject,  prsmth)
      combined.data <- merge(x=map, y=one.sample.data, by.x="marker", by.y="marker",
                             all.x=FALSE, all.y=TRUE)[,c(1:5,c(1,9,16,22,27,31,34,36,2,3,10,4,11,
                                                               17,5,12,18,23,6,13,19,24,28,7,14,
                                                               20,25,29,32,8,15,21,26,30,33,35)+5)]
      combined.data <- combined.data[combined.data$chr == chr[i],]
    
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
    #paste(is(all.subjects), collapse=" ")
    var.names <- c(names(all.subjects)[1:5], simple.het.labels)
    data.table::setnames(all.subjects, names(all.subjects), var.names)
    data.table::setkey(all.subjects, NULL)
    data.table::setkey(all.subjects, chr, bp, marker, subject)
    all.subjects$marker.notkey <- all.subjects$marker
    all.subjects$chr.notkey <- all.subjects$chr
    
    all.subjects[, f(all.subjects$chr, all.subjects$marker.notkey, 
                     all.subjects$AA, all.subjects$BB, all.subjects$CC, all.subjects$DD, 
                     all.subjects$EE, all.subjects$FF, all.subjects$GG, all.subjects$HH, 
                     all.subjects$AB, all.subjects$AC, all.subjects$BC, all.subjects$AD, 
                     all.subjects$BD, all.subjects$CD, all.subjects$AE, all.subjects$BE, 
                     all.subjects$CE, all.subjects$DE, all.subjects$AF, all.subjects$BF, 
                     all.subjects$CF, all.subjects$DF, all.subjects$EF, all.subjects$AG, 
                     all.subjects$BG, all.subjects$CG, all.subjects$DG, all.subjects$EG, 
                     all.subjects$FG, all.subjects$AH, all.subjects$BH, all.subjects$CH, 
                     all.subjects$DH, all.subjects$EH, all.subjects$FH, all.subjects$GH, 
                     allele.labels, diplotype.labels, full.to.dosages), by=all.subjects$marker]
    
    #-------------------------------------
    # make other necessary files
    #-------------------------------------
    one.subj = samples[1]
    markers.one.subj <- all.subjects[grepl(one.subj, subject), ]  
    
    markers.one.subj[, export_marker_name_file(all.subjects$chr.notkey, all.subjects$marker), by=all.subjects$chr]
    
    markers.one.subj[, export_marker_position_file(all.subjects$chr.notkey, all.subjects$bp), by=all.subjects$chr]
    
    markers.one.subj[, export_marker_chromosome_file(all.subjects$chr.notkey), by=all.subjects$chr]
    
    markers.one.subj[, export_marker_map_distance_file(all.subjects$chr.notkey, all.subjects$pos), by=all.subjects$chr]
  }
  
  # export file of mice names and strains
  #chr.names <- unique(all.subjects[,chr])
  subjects <- samples
  strains <- allele.labels
  for(this.chr in chr){
    save(subjects, file = paste0(final.output.path, 'additive/chr', this.chr, '/subjects.RData'))
    save(subjects, file = paste0(final.output.path, 'full/chr', this.chr, '/subjects.RData'))
    save(subjects, file = paste0(final.output.path, 'genotype/chr', this.chr, '/subjects.RData'))
    
    save(strains, file = paste0(final.output.path, 'additive/chr', this.chr, '/strains.RData'))
    save(strains, file = paste0(final.output.path, 'full/chr', this.chr, '/strains.RData'))
    save(strains, file = paste0(final.output.path, 'genotype/chr', this.chr, '/strains.RData'))
  }
}






