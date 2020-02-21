initialize.cisConnectome <- function(sce, flank.size = 10000, promoter.GR = NULL) {
  GR = rowRanges(sce)
  
  if(is.null(promoter.GR)) {
    genome.annotation = genome(GR)[[1]]
    if(genome.annotation %in% c('hg19', 'hg38', 'mm10', 'mm9')) {
      path = system.file(package="atACTIONdb", "data", sprintf('refGene_%s.rda', genome.annotation))
      load(path)  
    
      promoter.GR = refGene
      promoter.GR = promoter.GR[!is.na(promoter.GR$Gene)]
    } else {
      R.utils::printf('%s genome.annotation is not supported. Please provide promoter.GR directly', genome.annotation)
      return(sce)
    }
  }

  matches <- GenomicRanges::findOverlaps(GR, promoter.GR, select = "all", maxgap = flank.size)

  ii = S4Vectors::queryHits(matches)
  
  gg = sort(unique(unique(promoter.GR$Gene)))
  jj = match(as.character(promoter.GR$Gene[S4Vectors::subjectHits(matches)]), gg)
  
  Ind = Matrix::sparseMatrix(i = ii, j = jj, x = 1, dims = c(length(GR), length(gg)))
  colnames(Ind) = gg

  Ind@x[Ind@x > 1] = 1
  metadata(sce)[["cisConnectome"]] = Ind
  
  return(sce)
}

	
run.ATAC.ACTIONet <- function(sce, k_max = 20, layout.compactness = 50, thread_no = 8, epsilon = 3, LC = 1, arch.specificity.z = -1, core.z = 3, 
    sce.data.attr = "bin_counts", sym_method = "AND", scale.initial.coordinates = TRUE, reduction_slot = "ACTION", k_min = 2, n_epochs = 100, compute.core = F, compute.signature = F, specificity.mode = "dense", unification.min.cor = 0.95) {
	
	if(! ("bin_counts" %in% names(sce@assays)) )
		sce = add.binarized.counts(sce)
    
	
	ACTIONet.out = run.ACTIONet(sce, k_max = k_max, layout.compactness = layout.compactness, thread_no = thread_no, epsilon = epsilon, LC = LC, arch.specificity.z = arch.specificity.z, core.z = core.z, 
    sce.data.attr = sce.data.attr, sym_method = sym_method, scale.initial.coordinates = scale.initial.coordinates, reduction_slot = reduction_slot, k_min = k_min, n_epochs = n_epochs, compute.core = compute.core, compute.signature = compute.signature, specificity.mode = specificity.mode)

	return(ACTIONet.out)
}

