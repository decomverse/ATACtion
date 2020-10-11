plot_ATACtion_TF_view <- function(ace, TF.counts = 3, TF_slot = "TFs", renormalize = F) {
	TF.enrichment.mat = log1p(get_ATACtion_enrichment(ace, enrichment_name = TF_slot))
	#rownames(TF.enrichment.mat) = as.character(sapply(rownames(TF.enrichment.mat), function(str) str_split(str, "_")[[1]][[3]])) #colnames(motif_ix)
	
	plot.ACTIONet.feature.view(ace, feature.enrichment.table = TF.enrichment.mat, title = "TF View", top.features = TF.counts, renormalize = renormalize)
}

plot_ATACtion_gene_view <- function(ace, gene.counts = 3, gene_slot = "genes", renormalize = F) {
	gene.enrichment.mat = log1p(get_ATACtion_enrichment(ace, enrichment_name = gene_slot))
	plot.ACTIONet.feature.view(ace, feature.enrichment.table = gene.enrichment.mat, title = "Gene View", top.features = gene.counts, renormalize = renormalize)
}

plot_archetype_tracks <- function(input, chr, start, end, arch_names = NULL, arch_subset = NULL, profile_slot = "unified_feature_profile") {
	if(class(input) == "GRanges") {
		archGR = input
	} else {
		archGR = SummarizedExperiment::rowRanges(input)		
		X = rowMaps(ace)[[profile_slot]]		
	
		if(!is.null(arch_names)) {
			colnames(X) = paste(colnames(X), arch_names, sep = "-")
		}
		if(!is.null(arch_subset)) {
			X = X[, arch_subset]
		}
		mcols(archGR) = X
	}
	reference_genome = genome(archGR)[[1]]
	
	plotGR = makeGRangesFromDataFrame(data.frame(chr = chr, start = start, end = end))
	matches <- GenomicRanges::findOverlaps(archGR, plotGR, select = "all", maxgap = -1, minoverlap = 1)
	hits = queryHits(matches)
	if(length(hits) == 0) {
		warning(sprintf("No signal in the requested genomic regions %s_%d_%d", chr, start, end))
		return()
	}
	archGR = archGR[hits]
	

	library(Gviz)
	gtrack <- GenomeAxisTrack()
	itrack <- IdeogramTrack(genome = reference_genome, chromosome = chr)


  if(reference_genome == 'mm10') {
    data(mm10_genes)
	gene_model = mm10_genes
  } else if (reference_genome == 'mm9') {
    data(mm9_genes)
	gene_model = mm9_genes
  } else if(reference_genome == 'grch37' || reference_genome == 'hg19') {
    data(hg19_genes)
	gene_model = hg19_genes
  } else if(reference_genome == 'grch38' || reference_genome == 'hg38') {
    data(hg38_genes)
	gene_model = hg38_genes
  } else {
	  R.utils::printf('Unknown reference_genome %s\n', reference_genome);
	  return()
  }
	grtrack <- Gviz::GeneRegionTrack(gene_model, chromosome = chr, start = start, end = end, geneSymbols=TRUE, showId=TRUE, name = "Gene Model", shape="smallArrow")


	dTrack <- DataTrack(range = archGR, name = "Archetypes", showSampleNames = TRUE, cex.sampleNames = 0.6, type = c("horiz"), start = from, end = to, chromosome = chr, genome = reference_genome, groups = colnames(mcols(archGR)))


	plotTracks(list(dTrack, itrack, gtrack, grtrack), from = from, to = to)
}
