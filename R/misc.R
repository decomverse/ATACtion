annotatePeaks.sce <- function(sce, species) {
  library(ChIPseeker)
  species = tolower(species)
  
  if(species == 'mm10') {
    library(TxDb.Mmusculus.UCSC.mm10.knownGene)
    library(org.Mm.eg.db)
    txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
    peakAnno <- annotatePeak(rowRanges(sce), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
  }
  else if (species == 'mm9') {
    library(TxDb.Mmusculus.UCSC.mm9.knownGene)
    library(org.Mm.eg.db)
    txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
    peakAnno <- annotatePeak(rowRanges(sce), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Mm.eg.db")
  }
  else if(species == 'grch37' || species == 'hg19') {
    library(TxDb.Hsapiens.UCSC.hg19.knownGene)
    library(org.Hs.eg.db)
    txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
    peakAnno <- annotatePeak(rowRanges(sce), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
  }
  else if(species == 'grch38' || species == 'hg38') {
    library(TxDb.Hsapiens.UCSC.hg38.knownGene)
    library(org.Hs.eg.db)
    txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
    peakAnno <- annotatePeak(rowRanges(sce), tssRegion=c(-3000, 3000), TxDb=txdb, annoDb="org.Hs.eg.db")
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
  
  
  idx = which(!(colnames(mcols(sce)) %in% colnames(peak.annotations.df)))
  if(length(idx) > 0) {
    mcols(sce) = cbind(mcols(sce)[, idx], peak.annotations.df)
  } else {
    mcols(sce) = peak.annotations.df
  }
  
  perc = round(100*table(peak.annotations.df$region.type) / length(peak.annotations.df$region.type), 2)
  print(perc)
  
  return(sce)  
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


plot.ATACtion.motif.view <- function(ACTIONet.out, motif.enrichment, top.motif.count = 5, split = TRUE) {
	if(dim(motif.enrichment)[2] != length(ACTIONet.out$core.out$core.archs)) {
		motif.enrichment = motif.enrichment[, ACTIONet.out$core.out$core.archs]
	}

	if(split == TRUE) {
		parts = stringr::str_split_fixed(rownames(motif.enrichment),  "_", 5)
		motif.name = parts[, 3]
	} else {
		motif.name = rownames(motif.enrichment)
	}
	
	top.motifs = matrix(apply(motif.enrichment, 2, function(x) motif.name[order(x, decreasing = TRUE)[1:top.motif.count]]), ncol = dim(motif.enrichment)[2])
	unique.top.motifs = sort(unique(as.character(top.motifs)))


	arch.RGB = col2rgb(ACTIONet.out$arch.vis.out$colors[ACTIONet.out$core.out$core.archs]) / 256;
	arch.Lab = grDevices::convertColor(color= t(arch.RGB), from = 'sRGB', to = 'Lab')

	motif.color.Lab = t(sapply(unique.top.motifs, function(motif) {
	  v = as.numeric(apply(top.motifs, 2, function(mm) motif %in% mm))
	  v = v / sum(v)
	  motif.Lab = t(v) %*% arch.Lab
	  return(motif.Lab)
	}))
	motif.colors = rgb(grDevices::convertColor(color=motif.color.Lab, from = 'Lab', to = 'sRGB'))


	arch.coor = ACTIONet.out$arch.vis.out$coordinates[ACTIONet.out$core.out$core.archs, ]
	motif.coors = t(sapply(unique.top.motifs, function(motif) {
	  v = as.numeric(apply(top.motifs, 2, function(mm) motif %in% mm))
	  v = v / sum(v)
	  motif.coor = t(v) %*% arch.coor
	  return(motif.coor)
	}))


	motifs.df = data.frame(motif = unique.top.motifs, x = motif.coors[, 1], y = motif.coors[, 2])


	require(ggrepel)
	require(ggplot2)
	p <- ggplot(motifs.df, aes(x, y, label = motif, color=motif)) + scale_colour_manual(values=motif.colors) + geom_point(show.legend = FALSE) + geom_label_repel(show.legend = FALSE) + theme_void()

	plot(p)	

}


liftOverGR <- function(GR, from, to) {
  from = tolower(from)
  to = tolower(to)
  
  library(rtracklayer)
  library(ATACtiondb)
  
  path = system.file(package="ATACtiondb", "extdata/overlift", sprintf('%sTo%s.over.chain', from, to))
  
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
  library(ATACtiondb)
  
  path = system.file(package="ATACtiondb", "extdata/overlift", sprintf('%sTo%s.over.chain', from, to))
  
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
  library(ATACtiondb)
  
  path = system.file(package="ATACtiondb", "extdata/overlift", sprintf('%sTo%s.over.chain', from, to))
  
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
  library(ATACtiondb)
  
  path = system.file(package="ATACtiondb", "extdata/overlift", sprintf('%sTo%s.over.chain', from, to))
  
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


correlate.ATACtionet.with.Bulk <- function(ACTIONet.out, Bulk.sce, arch.assay = "signature") {
	require(GenomicRanges)
	
	Bulk.GR = rowRanges(Bulk.sce)
	
	matches <- GenomicRanges::findOverlaps(Bulk.GR, ACTIONet.out$archetype.accessibility, select = "all", maxgap = -1, minoverlap = 1)

	archetype.accessibility = ACTIONet.out$archetype.accessibility@assays[[arch.assay]]
	archetype.accessibility = orthoProject(archetype.accessibility, Matrix::rowMeans(archetype.accessibility))
	X = archetype.accessibility[S4Vectors::subjectHits(matches), ]
	X.orth = orthoProject(X, Matrix::rowMeans(X))
	
	Y = Bulk.sce@assays[["logcounts"]][S4Vectors::queryHits(matches), ]
	Y.orth = orthoProject(Y, Matrix::rowMeans(Y))

	CC = cor(as.matrix(X.orth), as.matrix(Y.orth))
	
	return(CC)
}

