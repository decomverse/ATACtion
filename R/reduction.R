reduce.peaks.inverse.ACTION <- function(sce, reduced_dim = 50, iters = 5, peak.subset = NULL) {
	if( !is.null(peak.subset) ) {
		sce = sce[peak.subset, ]
	}
		
	system.time( {reduction.out = reduceChromatinAccessibility(Matrix::t(sce@assays[["counts"]]), reduced_dim = reduced_dim, iters = iters)} )

	reduced.sce.inv <- SingleCellExperiment(assays = list(logcounts = reduction.out$S_r))
	colData(reduced.sce.inv)$GR = rowRanges(sce)
	reducedDim(reduced.sce.inv, "S_r") <- t(reduction.out$S_r)

	V = reduction.out$V
	colnames(V) = sapply(1:dim(V)[2], function(i) sprintf('PC%d', i))
	rownames(V) = colnames(sce)
	metadata(reduced.sce.inv)$V = V

	metadata(reduced.sce.inv)$ACTION.reduction.time = Sys.time()
	
	return(reduced.sce.inv)	
}


ATACtion.match.motifs <- function(sce) {
  library(chromVAR)
  library(chromVARmotifs)
  library(motifmatchr)
  library(Matrix)
  library(BiocParallel)
  
  GR = rowRanges(sce)

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
      motif_ix <- matchMotifs(motifs, sce, genome = BSgenome.Hsapiens.UCSC.hg19)
    }
    else if(species == 'hg38' || species == 'grch38') {
      library(BSgenome.Hsapiens.UCSC.hg38)
      motif_ix <- matchMotifs(motifs, sce, genome = BSgenome.Hsapiens.UCSC.hg38)
    } else {
      R.utils::printf('Species %s not supported. Please add motif dataset manually\n', species)
      return();
    }
  } 
  else if(grepl('mm', species)) {
    data(mouse_pwms_v2)
    motifs = human_pwms_v2
    
    if(species == 'mm10') {
      library(BSgenome.Mmusculus.UCSC.mm10)
      motif_ix <- matchMotifs(motifs, sce, genome = BSgenome.Mmusculus.UCSC.mm10)
    }
    else if(species == 'mm9') {
      library(BSgenome.Mmusculus.UCSC.mm9)
      motif_ix <- matchMotifs(motifs, sce, genome = BSgenome.Mmusculus.UCSC.mm9)
    } else {
      R.utils::printf('Species %s not supported. Please add motif dataset manually\n', species)
      return();
    }
  } 
  else {
    R.utils::printf('Species %s not supported. Please add motif dataset manually\n', species)
    return();
  }	
  
  matched.motifs = as(motif_ix@assays[["motifMatches"]], 'sparseMatrix')
  rownames(matched.motifs) = rownames(motif_ix)
  colnames(matched.motifs) = colnames(motif_ix)
  
  rowData(sce)[["matched.motifs"]] = matched.motifs
  
  return(sce)  
}


reduce.peaks.ACTION <- function(sce, reduced_dim = 50, iters = 5) {	
	if( !("cell.hashtag" %in% colnames(colData(sce))) ) {		
		print("tagging cells")
		time.tag = Sys.time()
		
		h = hashid_settings(salt = as.character(time.tag), min_length = 8)
		cell.hashtags = sapply(1:ncol(sce), function(i) ENCODE(i, h))
		sce$cell.hashtag = cell.hashtags

		colData(sce)$original.colnames = colnames(sce)
		colnames(sce) = sce$cell.hashtag
		
		metadata(sce)$tagging.time = time.tag		
	}
	
	reduction.out = reduceChromatinAccessibility(sce@assays[["counts"]], reduced_dim = reduced_dim, iters = iters)


	S_r = t(reduction.out$S_r)
	rownames(S_r) = colnames(sce)
	colnames(S_r) = sapply(1:ncol(S_r), function(i) sprintf('PC%d', i))	
	reducedDim(sce, "ACTION") <- S_r
	

    
    V = reduction.out$V
    colnames(V) = sapply(1:dim(V)[2], function(i) sprintf('PC%d', i))
    
    X = rowData(sce)
    PC.idx = -grep('^PC', colnames(X))
    if(length(PC.idx) > 0)
		X = X[, PC.idx]
	
    rowData(sce) = cbind(V, X)

	metadata(sce)$ACTION.reduction.time = Sys.time()
	return(sce)	
}

reduce.peaks.chromVAR <- function(sce, cores = 1) {
	library(chromVAR)
	library(chromVARmotifs)
	library(motifmatchr)
	library(Matrix)
	library(BiocParallel)
	if( !("cell.hashtag" %in% colnames(colData(sce))) ) {		
		print("tagging cells")
		time.tag = Sys.time()
		
		h = hashid_settings(salt = as.character(time.tag), min_length = 8)
		cell.hashtags = sapply(1:ncol(sce), function(i) ENCODE(i, h))
		sce$cell.hashtag = cell.hashtags

		colData(sce)$original.colnames = colnames(sce)
		colnames(sce) = sce$cell.hashtag
		
		metadata(sce)$tagging.time = time.tag		
	}
	
	if(cores == 1)
		register(SerialParam())
	else
		register(MulticoreParam(cores, progressbar = TRUE))	
		
	GR = rowRanges(sce)
	species = tolower(genome(GR))
	if(length(species) > 0)
		species = species[[1]]
	else {
	  print("Unknown genome");
	  return()
	}
  
	if( !("matched.motifs" %in% names(rowData(sce))) ) {
	  sce = ATACtion.match.motifs(sce)
	}  
  
	filtered.peaks = which(Matrix::rowSums(sce@assays[["counts"]]) == 0)
	if(length(filtered.peaks) > 0)
		sub.sce = sce[-filtered.peaks, ]
	else
		sub.sce = sce


	if(species == 'hg19' || species == 'grch37') {
		library(BSgenome.Hsapiens.UCSC.hg19)
		sub.sce <- addGCBias(sub.sce, genome = BSgenome.Hsapiens.UCSC.hg19)    
	}
	else if(species == 'hg38' || species == 'grch38') {
		library(BSgenome.Hsapiens.UCSC.hg38)
		sub.sce <- addGCBias(sub.sce, genome = BSgenome.Hsapiens.UCSC.hg38)    

	} 
	else if(species == 'mm10') {
		library(BSgenome.Mmusculus.UCSC.mm10)
		sub.sce <- addGCBias(sub.sce, genome = BSgenome.Mmusculus.UCSC.mm10)    
	}
	else if(species == 'mm9') {
		library(BSgenome.Mmusculus.UCSC.mm9)
		sub.sce <- addGCBias(sub.sce, genome = BSgenome.Mmusculus.UCSC.mm9)    
	} 
	else {
	  R.utils::printf('Species %s not supported. Please run chromVAR manually\n', species)
	  return();
	}

	bg <- getBackgroundPeaks(object = sub.sce)
  
	sce_chromVAR <- computeDeviations(object = sub.sce, annotations = rowData(sce)[["matched.motifs"]], background_peaks = bg)
	Z = sce_chromVAR@assays[['z']]

	filtered.rows = which(is.na(Matrix::rowSums(sce_chromVAR@assays[['z']])))
	if(length(filtered.rows) > 0)
	Z = Z[-filtered.rows, ]

	rownames(Z) = colnames(sce)
	sce@reducedDims[["chromVAR"]] = Z

	metadata(sce)$chromVAR.reduction.time = Sys.time()
	return(sce)
}


# Adopted from Cusanovich et al.
reduce.peaks.LSI <- function(sce, reduced_dim = 50, site_frequency_threshold = 0.0, logTF=FALSE, scale.factor=100000) {
	library(Matrix)
	library(irlba)	
	if( !("cell.hashtag" %in% colnames(colData(sce))) ) {		
		print("tagging cells")
		time.tag = Sys.time()
		
		h = hashid_settings(salt = as.character(time.tag), min_length = 8)
		cell.hashtags = sapply(1:ncol(sce), function(i) ENCODE(i, h))
		sce$cell.hashtag = cell.hashtags

		colData(sce)$original.colnames = colnames(sce)
		colnames(sce) = sce$cell.hashtag
		
		metadata(sce)$tagging.time = time.tag		
	}

	if(! ("bin_counts" %in% names(sce@assays)) ) {
		print("Adding binarized counts ... ")
		sce = add.binarized.counts(sce)
	}		
	atac_matrix = sce@assays[["bin_counts"]]

	if(site_frequency_threshold > 0) {
		rs = Matrix::rowSums(atac_matrix > 0)
		threshold = ncol(atac_matrix) * site_frequency_threshold
		atac_matrix = atac_matrix[rs >= threshold,]
	}


	#Calc TF-IDF
	print("Computing TF-IDF ...")
	
	npeaks <- Matrix::colSums(x = atac_matrix)
	tf <- Matrix::t(x = Matrix::t(x = atac_matrix) / npeaks)
	if(logTF){
		message("Epoch: running log term frequency ...");
        tf@x = log1p(tf@x * scale.factor);
	}	
	
		
	idf <- log(1+ ncol(x = atac_matrix) / Matrix::rowSums(x = atac_matrix))
	
	tfidf <- as(Diagonal(n = length(x = idf), x = as.vector(idf)), 'sparseMatrix') %*% tf
	tfidf[is.na(x = tfidf)] <- 0


	#Calc SVD then LSI
	print("Computing SVD using irlba...")
	svd <- irlba::irlba(tfidf, reduced_dim, reduced_dim)
	svdDiag <- matrix(0, nrow=reduced_dim, ncol=reduced_dim)
	diag(svdDiag) <- svd$d
	LSI <- t(svdDiag %*% t(svd$v))
	rownames(LSI) <- colnames(LSI)
	colnames(LSI) <- paste0("PC",seq_len(ncol(LSI)))


	LSI = as.matrix(LSI)
	rownames(LSI) = colnames(sce)
	colnames(LSI) = sapply(1:ncol(LSI), function(i) sprintf('LSI%d', i))	
	sce@reducedDims[["LSI"]] = LSI
	
	metadata(sce)$LSI.reduction.time = Sys.time()
	return(sce)
}
