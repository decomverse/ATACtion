run_ATACtion <- function(ace, batch = NULL, k_min = 2, k_max = 30, assay_name = "bin_counts", 
    reduction_slot = "ACTION", net_slot_out = "ACTIONet", min_cells_per_arch = 2, 
    max_iter_ACTION = 50, min_specificity_z_thresh = -3, network_density = 1, 
    mutual_edges_only = TRUE, layout_compactness = 50, layout_epochs = 1000, 
    layout_algorithm = 0, layout_in_parallel = TRUE, unification_violation_threshold = 0, 
    footprint_alpha = 0.85, thread_no = 0, full_trace = FALSE, 
    seed = 0) {
	
	ace = as(ace, "ACTIONetExperiment")
	
	if(! (data_slot %in% names(assays(ace))) & ("counts" %in% names(assays(ace)))) {
		B = as(assays(ace)[["counts"]], 'sparseMatrix')	
		B@x = rep(1, length(B@x))
		assays(ace)[[data_slot]] = B		
	}
    if( sum(grepl("cisConnectome", names(rowMaps(ace)))) == 0 ) {
		ace = add_proximal_peak_gene_interactions_to_ATACtion(ace, flank.size = flank.size)		
	}
	
	ace = run.ACTIONet(ace, batch = batch, k_min = k_min, k_max = k_max, assay_name = assay_name, 
    reduction_slot = reduction_slot, net_slot_out = net_slot_out, min_cells_per_arch = min_cells_per_arch, 
    max_iter_ACTION = max_iter_ACTION, min_specificity_z_thresh = min_specificity_z_thresh, network_density = network_density, 
    mutual_edges_only = mutual_edges_only, layout_compactness = layout_compactness, layout_epochs = layout_epochs, 
    layout_algorithm = layout_algorithm, layout_in_parallel = layout_in_parallel, unification_violation_threshold = unification_violation_threshold, 
    footprint_alpha = footprint_alpha, thread_no = thread_no, full_trace = full_trace, 
    seed = seed) 	

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


