run_ATACtion <- function(ace, assay_name = "bin_counts", flank.size = 10000, ...) {
	ace = as(ace, "ACTIONetExperiment")
	
	if(! (assay_name %in% names(assays(ace))) & ("counts" %in% names(assays(ace)))) {
		B = as(assays(ace)[["counts"]], 'sparseMatrix')	
		B@x = rep(1, length(B@x))
		assays(ace)[[assay_name]] = B		
	}
    if( sum(grepl("cisConnectome", names(rowMaps(ace)))) == 0 ) { # If there is no cisConenctome, add a base proximal one
		ace = add_proximal_peak_gene_interactions_to_ATACtion(ace, flank.size = flank.size)		
	}
	
	ace = run.ACTIONet(ace, assay_name = "bin_counts", ...) 	

	return(ace)
}


add_proximal_peak_gene_interactions_to_ATACtion <- function(ace, flank.size = 10000, promoter.GR = NULL, connectome_name = "proximal") {
  ace = as(ace, "ACTIONetExperiment")
  
  GR = SummarizedExperiment::rowRanges(ace)
  
  if(is.null(promoter.GR)) {
    genome.annotation = genome(GR)[[1]]
    if(genome.annotation %in% c('hg19', 'hg38', 'mm10', 'mm9')) {	  
	  path = system.file(package="ATACtionDB", "extdata/promoters", sprintf('%s.bed', genome.annotation))
      tbl = read.table(path, sep = "\t")
      colnames(tbl) = c("chr", "start", "end", "Gene")      
      gr = GenomicRanges::makeGRangesFromDataFrame(tbl, keep.extra.columns = TRUE)    
      gr <- sortSeqlevels(gr)
	  gr <- sort(gr)

      
      promoter.GR = gr
      promoter.GR = promoter.GR[!is.na(promoter.GR$Gene)]
    } else {
      R.utils::printf('%s genome.annotation is not supported. Please provide promoter.GR directly', genome.annotation)
      return(ace)
    }
  }

  matches <- GenomicRanges::findOverlaps(GR, promoter.GR, select = "all", maxgap = flank.size)

  ii = S4Vectors::queryHits(matches)
  
  gg = sort(unique(unique(promoter.GR$Gene)))
  jj = match(as.character(promoter.GR$Gene[S4Vectors::subjectHits(matches)]), gg)
  
  Ind = Matrix::sparseMatrix(i = ii, j = jj, x = 1, dims = c(length(GR), length(gg)))
  colnames(Ind) = gg

  Ind@x[Ind@x > 1] = 1
  rowMaps(ace)[[sprintf("%s_cisConnectome", connectome_name)]] = Ind
  
  return(ace)
}


add_physical_peak_gene_interactions_to_ATACtion <- function(ace, HiC.GI, promoter.GR=NULL, connectome_name = "physical") {
	ace = as(ace, "ACTIONetExperiment")

	GR = SummarizedExperiment::rowRanges(ace)
	
	if(is.null(promoter.GR)) {
		genome.annotation = genome(GR)[[1]]
		if(genome.annotation %in% c('hg19', 'hg38', 'mm10', 'mm9')) {
			path = system.file(package="ATACtionDB", "extdata/promoters", sprintf('%s.bed', genome.annotation))
			tbl = read.table(path, sep = "\t")
			colnames(tbl) = c("chr", "start", "end", "Gene")      
			gr = GenomicRanges::makeGRangesFromDataFrame(tbl, keep.extra.columns = TRUE)    
			gr <- sortSeqlevels(gr)
			gr <- sort(gr)


			promoter.GR = gr
			promoter.GR = promoter.GR[!is.na(promoter.GR$Gene)]
		} else {
			R.utils::printf('%s genome.annotation is not supported. Please provide promoter.GR directly', genome.annotation)
			return(ace)
		}
	}
	
	
	E_mask1 <- GenomicRanges::findOverlaps(GR, anchorOne(GI), select = "all", maxgap = -1, minoverlap = 1)
	P_mask1 <- GenomicRanges::findOverlaps(promoter.GR, anchorTwo(GI), select = "all", maxgap = -1, minoverlap = 1)
	
	Emat1 = sparseMatrix(i = S4Vectors::queryHits(E_mask1), j = S4Vectors::subjectHits(E_mask1), x = GI$scores[S4Vectors::subjectHits(E_mask1)], dims = c(length(GR), length(GI)))
	Pmat1 = sparseMatrix(i = S4Vectors::queryHits(P_mask1), j = S4Vectors::subjectHits(P_mask1), x = GI$scores[S4Vectors::subjectHits(P_mask1)], dims = c(length(promoter.GR), length(GI)))
	E2P1 = Emat1 %*% Matrix::t(Pmat1)
	
	
	E_mask2 <- GenomicRanges::findOverlaps(GR, anchorTwo(GI), select = "all", maxgap = -1, minoverlap = 1)
	P_mask2 <- GenomicRanges::findOverlaps(promoter.GR, anchorOne(GI), select = "all", maxgap = -1, minoverlap = 1)
	Emat2 = sparseMatrix(i = S4Vectors::queryHits(E_mask2), j = S4Vectors::subjectHits(E_mask2), x = GI$scores[S4Vectors::subjectHits(E_mask2)], dims = c(length(GR), length(GI)))
	Pmat2 = sparseMatrix(i = S4Vectors::queryHits(P_mask2), j = S4Vectors::subjectHits(P_mask2), x = GI$scores[S4Vectors::subjectHits(P_mask2)], dims = c(length(promoter.GR), length(GI)))
	E2P2 = Emat2 %*% Matrix::t(Pmat2)
	
	E2P = E2P1 + E2P2
	E2P@x = rep(1, length(E2P@x))
	
    rowMaps(ace)[[sprintf("%s_cisConnectome", connectome_name)]] = E2P
	
	return(ace)
}


