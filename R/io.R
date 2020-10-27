export_ATACtion_to_cisTopic <- function(ace, project.name = "ACTIONet") {
	library(cisTopic)
	mat = as(assays(ace)[["counts"]], 'dgCMatrix')

	GR = SummarizedExperiment::rowRanges(ace)

	chr = as.character(GR@seqnames)
	start = GR@ranges@start
	end = GR@ranges@start + GR@ranges@width - 1
	coor.str = sapply(1:length(GR), function(i) sprintf('%s:%d-%d', chr[i], start[i], end[i]))
	rownames(mat) = coor.str

	cisTopicObject <- createcisTopicObject(mat, project.name=project.name)

	meta.data = as.data.frame(colData(ace))
	rownames(meta.data) = colnames(ace)

	cisTopicObject@cell.data = cbind(cisTopicObject@cell.data, meta.data)

	return(cisTopicObject)
}



export_ATACtion_to_SnapATAC <-function(ace) {
	library(SnapATAC)
	mat = as(t(assays(ace)[["counts"]]), 'dgCMatrix')
	barcodes = colnames(assays(ace)[["counts"]])
	peaks = SummarizedExperiment::rowRanges(ace)

	x.sp = createSnapFromPmat(mat, barcodes, peaks)

	x.sp@metaData = as.data.frame(colData(ace))
	return(x.sp)
}

import_ATACtion_from_10X <- function(input_path, matrix_file = "matrix.mtx", bed_file = "peaks.bed", sample_annotation_file = NULL, coordinate = NULL, bed.header = FALSE, sample_annotations.header = TRUE, sep = '\t') {
	library(Matrix)
	library(SingleCellExperiment)
	library(GenomicRanges)

	if(is.null(coordinate)) {
		print("Please provide genome coordinate system (i.e., hg19)")
		return()
	}

	peaks = readMM(paste(input_path, matrix_file, sep='/'))


	# BED = as.data.frame(stringr::str_split_fixed(BED.str[, 1],  ":|-|_", 3), stringsAsFactors = FALSE)
	BED = read.table(paste(input_path, bed_file, sep='/'), sep, header = bed.header, as.is = TRUE)
	colnames(BED) = c('chr', 'start', 'end')


	GR = makeGRangesFromDataFrame(BED)
	peak_names = paste(GR@seqnames, GR@ranges@start, GR@ranges@start+GR@ranges@width, sep = '_')
	rownames(peaks) = peak_names

	genome(GR) = coordinate

	if(!is.null(sample_annotation_file)) {
		sample_annotations = read.table(paste(input_path, 'sample_annotations.txt', sep='/'), header = sample_annotations.header, as.is = TRUE, sep = sep)
		ace <- ACTIONetExperiment(assays = list(counts = peaks), colData = sample_annotations, rowRanges = GR)
	} else {
		ace <- ACTIONetExperiment(assays = list(counts = peaks), rowRanges = GR)
	}

	return(ace)
}


import_ATACtion_from_SnapATAC <- function(x.sp) {
	ace.sp = SingleCellExperiment(assays = list(counts = Matrix::t(x.sp@pmat)), rowRanges = GR)

	meta = x.sp@metaData
	if(length(meta) > 0) {
		colData(ace.sp) = DataFrame(meta)
	}
	G = x.sp@graph@mat
	if(length(G) > 0) {
		colNets(ace.sp)$knn = G@mat
	}
	if(length(x.sp@cluster)) {
		ace.sp$clusters = x.sp@cluster
	}
	if(length(x.sp@tsne) > 0) {
		colMaps(x.sp)$tsne = x.sp@tsne
	}
	if(length(x.sp@umap) > 0) {
		colMaps(x.sp)$umap = x.sp@umap
	}

	return(x.sp)
}

import_ATACtion_from_ArchR <- function(proj, genome.name = "hg19") {
	require(ACTIONet)
	require(BiocParallel)

	proj_GeneScoreMatrix = getMatrixFromProject(ArchRProj = proj, useMatrix = "GeneScoreMatrix")
	proj_TileMatrix = getMatrixFromProject(ArchRProj = proj, useMatrix = "TileMatrix", binarize = TRUE)

	ace = as(proj_TileMatrix, "ACTIONetExperiment")
	names(assays(ace)) = c("bin_counts")
	colMaps(ace)$GeneScore = Matrix::t(assays(proj_GeneScoreMatrix)$GeneScoreMatrix)

	x = as.numeric(rowData(ace)$start)
	width = x[2]-x[1]
	BED = data.frame(chr = as.character(rowData(ace)$seqnames), start = x, end = x + (width-1))
	GR = GenomicRanges::makeGRangesFromDataFrame(BED)
	SummarizedExperiment::rowRanges(ace) = GR

	rownames(ace) = paste(BED$chr, BED$start, BED$end, sep = "_")

	genome(SummarizedExperiment::rowRanges(ace)) = genome.name
	return(ace)
}


impute_genomewide_ATAC_signals <- function(ace, profile.slot = "unified_feature_specificity", reference_genome = "hg38") {
	  if(reference_genome == 'mm10') {
		library(BSgenome.Mmusculus.UCSC.mm10)
		genome <- BSgenome.Mmusculus.UCSC.mm10
	  }
	  else if (reference_genome == 'mm9') {
		library(BSgenome.Mmusculus.UCSC.mm9)
		genome <- BSgenome.Mmusculus.UCSC.mm9
	  }
	  else if(reference_genome == 'grch37' || reference_genome == 'hg19') {
		library(BSgenome.Hsapiens.UCSC.hg19)
		genome <- BSgenome.Hsapiens.UCSC.hg19
	  }
	  else if(reference_genome == 'grch38' || reference_genome == 'hg38') {
		library(BSgenome.Hsapiens.UCSC.hg38)
		genome <- BSgenome.Hsapiens.UCSC.hg38
	  }
	  else {
		  R.utils::printf('Unknown reference_genome %s\n', reference_genome);
		  return()
	  }	
	profile = rowMaps(ace)[[profile.slot]]
	GR = rowRanges(ace)
	mcols(GR) = profile
	GR.chr = split(GR, seqnames(GR))

	chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
	chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
	
	chr_size = width(chromSizes)
	names(chr_size) = as.character(seqnames(chromSizes))
	
	bins = bin_genome(reference_genome, chunk_size)
	bins.chr = split(bins, seqnames(bins))

	common.chr = intersect(names(frags.chr), names(bins.chr))
	bins.chr = bins.chr[common.chr]	
	frags.chr = frags.chr[common.chr]
	
	imputed_signal = lapply(common.chr, function(chr) {
		bin.GR = bins.chr[[chr]]
		cur.GR = GR.chr[[chr]]
		
		Y = as.matrix(mcols(cur.GR))
		x = start(cur.GR) + (end(cur.GR) - start(cur.GR))/2
		
		W = apply(y, 2, function(y) {
			sm = smooth.spline(x, y) 
			plot(sm$x, sm$y)
			yy = predict(sm, x = 1:chr_size[chr])
			w = pmax(0, yy$y)
			return(w)			
		})
		BED = data.frame(chr = rep(chr, chr_size[[chr]]), start = 1, end = chr_size[[chr]])
		gg = makeGRangesFromDataFrame(BED)
		genome(gg) = reference_genome
		mcols(gg) = W
		
		return(gg)
	})
	
	names(imputed_signal) = common.chr
	
	return(imputed_signal)
}
