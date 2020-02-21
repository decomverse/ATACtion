is.sce <- function(sce) {
	return(length(which(is(sce)=="SingleCellExperiment"))!=0)
}

is.GR <- function(GR) {
	return(length(which(is(GR)=="GenomicRanges"))!=0)
}

add.binarized.counts <- function(sce) {
    require(Matrix)
    B = as(sce@assays[["counts"]], 'sparseMatrix')	
    B@x = rep(1, length(B@x))
    sce@assays[["bin_counts"]] = B
    
    return(sce)
}
