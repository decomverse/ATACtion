sce2SNAP <-function(sce) {
	mat = as(t(sce@assays[["counts"]]), 'dgCMatrix')
	barcodes = colnames(sce@assays[["counts"]])
	peaks = rowRanges(sce)  

	x.sp = createSnapFromPmat(mat, barcodes, peaks)
	x.sp = makeBinary(x.sp, mat="pmat");

	return(x.sp)
}

import.peaks.10X <- function(input_path, matrix_file = "matrix.mtx", bed_file = "peaks.bed", sample_annotation_file = NULL, coordinate = NULL, bed.header = FALSE, sample_annotations.header = TRUE, sep = '\t') {
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
		sce <- SingleCellExperiment(assays = list(counts = peaks), colData = sample_annotations, rowRanges = GR)   	
	} else {
		sce <- SingleCellExperiment(assays = list(counts = peaks), rowRanges = GR)   	
	}
	
	return(sce)
}


add.TFmat.10X <- function(sce, input_path, matrix_file = "matrix.mtx", motif_file = "motifs.tsv", motif.col = 2) {
	library(Matrix)
	library(SingleCellExperiment)
	
	scores = readMM(paste(input_path, matrix_file, sep='/'));
	motifs = read.table(paste(input_path, motif_file, sep='/'), as.is = TRUE)
	
	rownames(scores) = motifs[, motif.col]
	
	A = as(scores, 'dgTMatrix')
	cs = Matrix::colSums(A)    
	B = Matrix::sparseMatrix(i = A@i+1, j = A@j+1, x = log(1 + median(cs)*(A@x / cs[A@j + 1])), dims = dim(A))
	rownames(B) = rownames(A)
	
	sce@reducedDims[["TFs"]] = as.matrix(Matrix::t(B))
			
	
	return(sce)
}

load.gwascat.GR <- function(current = TRUE, coordinate = 'hg38') {	
	library(gwascat)
	
	if(current == TRUE)
		suppressWarnings({ebicat38 = gwascat::makeCurrentGwascat()})
	else
		data("ebicat38", package = "gwascat", envir = environment())
	
	if(coordinate == 'hg19') { # Liftover
		ebicat19 = liftOverGR(ebicat38, "hg38", "hg19")
		return(as(ebicat19, "GRanges"))
	} else {
		return(as(ebicat38, "GRanges"))
	}
}

