run_ATACtion <- function(ace, k_max = 30, min.cells.per.arch = 2, min_specificity_z_threshold = -3,
    network_density = 1, mutual_edges_only = TRUE, layout_compactness = 50, layout_epochs = 500,
    layout.in.parallel = FALSE, thread_no = 0, data_slot = "bin_counts", reduction_slot = "ACTION",
    unification_min_edge_weight = 0.5, unification_min_coreness = 2, unification_resolution = 1, unification_min_repeat = 0, unification_alpha = 0.05, unification_beta = 0.5,    
    footprint_alpha = 0.85, max_iter_ACTION = 50, full.trace = FALSE) {   
		
	
	if(! (data_slot %in% names(assays(ace))) & ("counts" %in% names(assays(ace)))) {
		B = as(assays(ace)[["counts"]], 'sparseMatrix')	
		B@x = rep(1, length(B@x))
		assays(ace)[[data_slot]] = B		
	}
    if( sum(grepl("cisConnectome", names(rowMaps(ace)))) == 0 ) {
		ace = add_proximal_peak_gene_interactions_to_ATACtion(ace)		
	}
	
	ace = run.ACTIONet(ace, k_max = k_max, min.cells.per.arch = min.cells.per.arch , min_specificity_z_threshold = min_specificity_z_threshold, 
    network_density = network_density, mutual_edges_only = mutual_edges_only, layout_compactness = layout_compactness, 
    layout_epochs = layout_epochs, layout.in.parallel = layout.in.parallel, thread_no = thread_no, 
    data_slot = data_slot, reduction_slot = reduction_slot,
    unification_min_edge_weight = unification_min_edge_weight, unification_min_coreness = unification_min_coreness, unification_resolution = unification_resolution, unification_min_repeat = unification_min_repeat, unification_alpha = unification_alpha, unification_beta = unification_beta,     
    footprint_alpha = footprint_alpha, max_iter_ACTION = max_iter_ACTION, full.trace = full.trace) 	

	return(ace)
}


add_proximal_peak_gene_interactions_to_ATACtion <- function(ace, flank.size = 10000, promoter.GR = NULL, connectome_name = "proximal") {
  GR = rowRanges(ace)
  
  if(is.null(promoter.GR)) {
    genome.annotation = genome(GR)[[1]]
    if(genome.annotation %in% c('hg19', 'hg38', 'mm10', 'mm9')) {
      path = system.file(package="ATACtionDB", "data", sprintf('refGene_%s.rda', genome.annotation))
      load(path)  
    
      promoter.GR = refGene
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
	GR = rowRanges(ace)
	
	if(is.null(promoter.GR)) {
		genome.annotation = genome(GR)[[1]]
		if(genome.annotation %in% c('hg19', 'hg38', 'mm10', 'mm9')) {
			path = system.file(package="ATACtionDB", "data", sprintf('refGene_%s.rda', genome.annotation))
			load(path)  
			
			promoter.GR = refGene
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


