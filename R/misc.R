is.GR <- function(GR) {
	return(length(which(is(GR)=="GenomicRanges"))!=0)
}

add_motif_matched_to_ATACtion <- function(ace) {
  library(chromVAR)
  library(chromVARmotifs)
  library(motifmatchr)
  library(Matrix)
  library(BiocParallel)
  
  GR = SummarizedExperiment::rowRanges(ace)

  species = tolower(genome(GR))
  if(length(species) > 0)
	species = species[[1]]
  else {
	  print("Unknown genome");
	  return()
  }
  
  if(grepl('hg', species)) {
    data(human_pwms_v2)
    motifs = human_pwms_v2
    
    if(species == 'hg19' || species == 'grch37') {
      library(BSgenome.Hsapiens.UCSC.hg19)
      motif_ix <- matchMotifs(motifs, ace, genome = BSgenome.Hsapiens.UCSC.hg19)
    }
    else if(species == 'hg38' || species == 'grch38') {
      library(BSgenome.Hsapiens.UCSC.hg38)
      motif_ix <- matchMotifs(motifs, ace, genome = BSgenome.Hsapiens.UCSC.hg38)
    } else {
      R.utils::printf('Species %s not supported. Please add motif dataset manually\n', species)
      return(ace);
    }
  } 
  else if(grepl('mm', species)) {
    data(mouse_pwms_v2)
    motifs = human_pwms_v2
    
    if(species == 'mm10') {
      library(BSgenome.Mmusculus.UCSC.mm10)
      motif_ix <- matchMotifs(motifs, ace, genome = BSgenome.Mmusculus.UCSC.mm10)
    }
    else if(species == 'mm9') {
      library(BSgenome.Mmusculus.UCSC.mm9)
      motif_ix <- matchMotifs(motifs, ace, genome = BSgenome.Mmusculus.UCSC.mm9)
    } else {
      R.utils::printf('Species %s not supported. Please add motif dataset manually\n', species)
      return(ace);
    }
  } 
  else {
    R.utils::printf('Species %s not supported. Please add motif dataset manually\n', species)
    return(ace);
  }	
  
  library(stringr)
  motif_matches = as(assays(motif_ix)[["motifMatches"]], 'sparseMatrix')
  rownames(motif_matches) = rownames(motif_ix)
  colnames(motif_matches) = as.character(sapply(colnames(motif_ix), function(str) str_split(str, "_")[[1]][[3]])) #colnames(motif_ix)
  
  
  rowMaps(ace)[["motif_matches"]] = as(motif_matches, "dgCMatrix")
  
  return(ace)  
}




annotate_ATACTion_peaks <- function(ace) {
  library(ChIPseeker)
  
  GR = SummarizedExperiment::rowRanges(ace)
  species = tolower(genome(GR)[[1]])
  
  if(species == 'mm10') {
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    library(org.Mm.eg.db)
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    peakAnno <- annotatePeak(GR, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
  }
  else if (species == 'mm9') {
    library(TxDb.Mmusculus.UCSC.mm9.knownGene)
    library(org.Mm.eg.db)
    txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
    peakAnno <- annotatePeak(GR, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
  }
  else if(species == 'grch37' || species == 'hg19') {
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(org.Hs.eg.db)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    peakAnno <- annotatePeak(GR, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
  }
  else if(species == 'grch38' || species == 'hg38') {
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    peakAnno <- annotatePeak(GR, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
  }
  else {
	  R.utils::printf('Unknown species %s\n', species);
	  return()
  }
  peak.annotations.df <- as.data.frame(peakAnno)
  peak.annotations.df = peak.annotations.df[, -c(1:5)]
  if(!grepl('chr', peak.annotations.df$geneChr[[1]]))
    peak.annotations.df$geneChr = paste('chr', peak.annotations.df$geneChr, sep='')
  peak.annotations.df$gene = peak.annotations.df$SYMBOL
  peak.annotations.df$region.type = as.character(sapply(peak.annotations.df$annotation, function(x) strsplit(x, ' (', fixed = TRUE)[[1]][[1]]))
  
  
  idx = which(!(colnames(mcols(ace)) %in% colnames(peak.annotations.df)))
  if(length(idx) > 0) {
    mcols(ace) = cbind(mcols(ace)[, idx], peak.annotations.df)
  } else {
    mcols(ace) = peak.annotations.df
  }
  
  perc = round(100*table(peak.annotations.df$region.type) / length(peak.annotations.df$region.type), 2)
  print(perc)
  
  return(ace)  
}


annotate.GR.peaks <- function(GR) {
  library(ChIPseeker)
  species = tolower(genome(GR)[[1]])
  print(species)
  
  if(species == 'mm10') {
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    library(org.Mm.eg.db)
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    peakAnno <- annotatePeak(GR, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
  }
  else if (species == 'mm9') {
    library(TxDb.Mmusculus.UCSC.mm9.knownGene)
    library(org.Mm.eg.db)
    txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
    peakAnno <- annotatePeak(GR, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
  }
  else if(species == 'grch37' || species == 'hg19') {
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(org.Hs.eg.db)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    peakAnno <- annotatePeak(GR, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
  }
  else if(species == 'grch38' || species == 'hg38') {
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    peakAnno <- annotatePeak(GR, tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
  }
  else {
	  R.utils::printf('Unknown species %s\n', species);
	  return()
  }
  peak.annotations.df <- as.data.frame(peakAnno)
  peak.annotations.df = peak.annotations.df[, -c(1:5)]
  if(!grepl('chr', peak.annotations.df$geneChr[[1]]))
    peak.annotations.df$geneChr = paste('chr', peak.annotations.df$geneChr, sep='')
  peak.annotations.df$gene = peak.annotations.df$SYMBOL
  peak.annotations.df$region.type = as.character(sapply(peak.annotations.df$annotation, function(x) strsplit(x, ' (', fixed = TRUE)[[1]][[1]]))
  
  return(peak.annotations.df)  
}

liftOverGR <- function(GR, from, to) {
  from = tolower(from)
  to = tolower(to)
  
  library(rtracklayer)
  library(ATACtionDB)
  
  path = system.file(package="ATACtionDB", "extdata/overlift", sprintf('%sTo%s.over.chain', from, to))
  
  ch = import.chain(path)
  seqlevelsStyle(GR) = "UCSC"  # necessary
  liftedGR = liftOver(GR, ch)
  liftedGR = unlist(liftedGR)
  genome(liftedGR) = to
  return(liftedGR)
  
}

liftOverGR_list <- function(GR_list, from, to) {
  from = tolower(from)
  to = tolower(to)
  
  library(rtracklayer)
  library(ATACtionDB)
  
  path = system.file(package="ATACtionDB", "extdata/overlift", sprintf('%sTo%s.over.chain', from, to))
  
  ch = import.chain(path)
  
  liftedGR_list = lapply(GR_list, function(GR) {
    seqlevelsStyle(GR) = "UCSC"
    liftedGR = liftOver(GR, ch)
    liftedGR = unlist(liftedGR)
    genome(liftedGR) = to
    return(liftedGR)
  })
  
  return(liftedGR_list)
}


liftOverGR <- function(GR, from, to) {
  from = tolower(from)
  to = tolower(to)
  
  library(rtracklayer)
  library(ATACtionDB)
  
  path = system.file(package="ATACtionDB", "extdata/overlift", sprintf('%sTo%s.over.chain', from, to))
  
  ch = import.chain(path)
  seqlevelsStyle(GR) = "UCSC"  # necessary
  liftedGR = liftOver(GR, ch)
  liftedGR = unlist(liftedGR)
  genome(liftedGR) = to
  return(liftedGR)
  
}

liftOverGR_list <- function(GR_list, from, to) {
  from = tolower(from)
  to = tolower(to)
  
  library(rtracklayer)
  library(ATACtionDB)
  
  path = system.file(package="ATACtionDB", "extdata/overlift", sprintf('%sTo%s.over.chain', from, to))
  
  ch = import.chain(path)
  
  liftedGR_list = lapply(GR_list, function(GR) {
    seqlevelsStyle(GR) = "UCSC"
    liftedGR = liftOver(GR, ch)
    liftedGR = unlist(liftedGR)
    genome(liftedGR) = to
    return(liftedGR)
  })
  
  return(liftedGR_list)
}

impute_genome_wide_activity <- function(ace, chr, start, end, arch_subset = NULL, bp_resolution = 1, kernel_type = "matern_5_2", profile_slot = "unified_feature_profile") {
	library(FastGaSP)
	
	archGR = SummarizedExperiment::rowRanges(ace)
	species = genome(archGR)[[1]]
	Y = rowMaps(ace)[[profile_slot]]
	
	plotGR = makeGRangesFromDataFrame(data.frame(chr = chr, start = start, end = end))
	matches <- GenomicRanges::findOverlaps(archGR, plotGR, select = "all", maxgap = -1, minoverlap = 1)
	hits = queryHits(matches)
	if(length(hits) == 0) {
		warning(sprintf("No signal in the requested genomic regions %s_%d_%d", chr, start, end))
		return()
	}
	archGR = archGR[hits]
	Y = Y[hits, ]
	
	perm = order(start(archGR), decreasing = FALSE)
	archGR = archGR[perm]
	Y = Y[perm, ]
	
	if(!is.null(arch_subset)) {
		Y = Y[, arch_subset]
	}
	
	s = start(archGR)
	e = end(archGR)
	r = e - s
	x = s + round(r/2)
	
	bins = seq(start, end, by = bp_resolution)
	

	imputed_tracks = matrix(0, nrow=length(bins), ncol = ncol(Y))
	colnames(imputed_tracks) = colnames(Y)

	for(j in 1:ncol(Y)) {		
		y = Y[, j]
		
		fgasp.model = fgasp(input = x, output = y, have_noise = TRUE, kernel_type = kernel_type)
		est_all<-optim(c(log(1),log(.02)),log_lik,object=fgasp.model,method="L-BFGS-B",
	control = list(fnscale=-1))
		
		pred.model=predict(param=est_all$par,object=fgasp.model, testing_input=bins)
		imputed_tracks[, j] = pred.model@mean
		
		# lb=pred.model@mean+qnorm(0.025)*sqrt(pred.model@var)
		# ub=pred.model@mean+qnorm(0.975)*sqrt(pred.model@var)
		# 
		# plot(pred.model@testing_input,pred.model@mean,type='l',col='blue',
		# xlab='x',ylab='y')
	}	
	
	
	imputed_GR = GenomicRanges::makeGRangesFromDataFrame(data.frame(chr = chr, start = bins, end = bins + (bins[2]-bins[1]) + 1))
	mcols(imputed_GR) = imputed_tracks
	
	return(imputed_GR)
}
