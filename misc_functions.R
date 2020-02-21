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

annotate.scATAC.cells.using.markers <- function(sce, marker.genes, proximity = 10000) {
	library(Matrix)
	
	peak.annotations = rowData(sce)
	celltype.mask = sapply(marker.genes, function(genes) {
	  mask = (peak.annotations$gene %in% genes) & (abs(peak.annotations$distanceToTSS) < proximity)
	  return(mask)
	})

	mask = rowSums(celltype.mask) > 0

	pop.size = sum(mask)
	sample.sizes = Matrix::colSums(sce@assays[["bin_counts"]][mask, ])
		
	success.sizes = Matrix::t(sce@assays[["bin_counts"]][mask, ]) %*% celltype.mask[mask, ]
	pos.size = Matrix::colSums(celltype.mask)

	celltype.scores = sapply(1:length(marker.genes), function(i) {
	  logPvals = HGT_tail(pop.size, pos.size[i], sample.sizes, success.sizes[, i])
	})
	colnames(celltype.scores) = names(marker.genes)

	Labels.conf = apply(celltype.scores, 1, max)
	Labels = colnames(celltype.scores)[apply(celltype.scores, 1, which.max)]
	Labels = factor(Labels, levels = sort(unique(Labels)))
		
	annotations = list(Labels = Labels, confidence = Labels.conf)
	
	return(annotations)	
}

peakset.Enrichment <- function(peaks.signature.profile, peakset.ind.mat) {
	# Normalize scores to avoid heavy-tail side-effect(s)
	peaks.signature.profile[peaks.signature.profile < 0] = 0
	scores = scale(peaks.signature.profile, center = FALSE, scale = Matrix::colSums(peaks.signature.profile))

	# Oberved statistic
	Obs = t(scores) %*% peakset.ind.mat

	# Expectation
	p_c = Matrix::colMeans(peakset.ind.mat)
	rho = p_c / mean(p_c)

	p_r = Matrix::rowMeans(peakset.ind.mat)
	Exp = t(scores) %*% p_r %*% t(rho)

	Nu = t(scores^2) %*% p_r %*% t(rho^2)

	Lambda = Obs - Exp;

	a = apply(scores, 2, max) %*% t(rho)

	# Chernoff tail bound
	log.pval = (Lambda^2) / (2 * (Nu + a*Lambda/3))
	log.pval[Lambda  < 0] = 0
	log.pval[is.na(log.pval)] = 0

	return(log.pval)
}

annotate.scATAC.archetypes.using.markers <- function(sce, archetype.accessibility.sce, markers, proximity = 20000) {
	peak.annotations = rowData(sce)
	celltype.mask = sapply(marker.genes, function(genes) {
	  mask = (peak.annotations$gene %in% genes) & (abs(peak.annotations$distanceToTSS) < proximity)
	  return(mask)
	})
	
	celltype.Enrichment = peakset.Enrichment(archetype.accessibility.sce@assays[["significance"]], celltype.mask)
	
	Labels = colnames(celltype.Enrichment)[apply(celltype.Enrichment, 1, which.max)]
	Labels = factor(Labels, levels = sort(unique(Labels)))
	
	out.list = list(Labels = Labels, Enrichment = celltype.Enrichment)
	
	return(out.list)
}


atACTIONet.geneset.Enrichment <- function(ACTIONet.out, ind.mat, max.dist2TSS = 2.5e3) {
	require(Matrix)

  accessibility.sce = ACTIONet.out$archetype.accessibility
    
  proximal.mask = (abs(rowData(accessibility.sce)$distanceToTSS) <= max.dist2TSS)
  proximal.accessibility.sce = accessibility.sce[proximal.mask, ]
  proximal.accessibility.significance = proximal.accessibility.sce@assays[["significance"]]
  
  proximal.ATAC.genes = rowData(proximal.accessibility.sce)$gene
  annotation.genes = rownames(ind.mat)
  
  common.genes = intersect(annotation.genes, proximal.ATAC.genes)
  rows = match(common.genes, annotation.genes)
  ind.mat = ind.mat[rows, ]
  
  ind.mat = as(ind.mat, 'dgTMatrix')
  

  ii = (ind.mat@i + 1)
  jj = (ind.mat@j + 1)
  gg = common.genes[ii] 
  gg.ii = match(gg, proximal.ATAC.genes)
  
  ATAC.ind.mat = sparseMatrix(i = gg.ii, j = jj, x = 1, dims = c(length(proximal.ATAC.genes), dim(ind.mat)[2]))
  
  colnames(ATAC.ind.mat) = colnames(ind.mat)
  rownames(ATAC.ind.mat) = proximal.ATAC.genes
  

	A = log(1+proximal.accessibility.significance)
	X = ATAC.ind.mat

	p_c = Matrix::colMeans(X)


	Obs = as.matrix(Matrix::t(X) %*% A)
	Exp = as.matrix(p_c %*% Matrix::t(Matrix::colSums(A)))
	Nu = as.matrix(p_c %*% Matrix::t(Matrix::colSums(A^2)))

	Lambda = Obs - Exp
	
	a = apply(A, 2, max)
	ones = array(1, dim = dim(Nu)[1])
  
	logPvals = Lambda^2 / (2*(Nu + (ones %*% t(a))*Lambda/3))
	
	rownames(logPvals) = colnames(X)
	colnames(logPvals) = colnames(A)

	logPvals[Lambda  < 0] = 0
	logPvals[is.na(logPvals)] = 0
	
	logPvals = logPvals / log(10)
	
		
	return(logPvals)
}


atACTIONet.motif.Enrichment <- function(ACTIONet.out, core_no = 1) {

	A = log(1+ACTIONet.out$archetype.accessibility@assays[["significance"]])	
	if( !("matched.motifs" %in% names(ACTIONet.out)) ) {
		print("Motif matching matrix doesn't exist in the ACTIONet.out. Running motif matching first");
		ACTIONet.out = atACTIONet.match.motifs(ACTIONet.out, core_no)
	}
	X = as(ACTIONet.out$matched.motifs, 'dgCMatrix')

	p_c = Matrix::colMeans(X)


	Obs = as.matrix(Matrix::t(X) %*% A)
	Exp = as.matrix(p_c %*% t(Matrix::colSums(A)))
	Nu = as.matrix(p_c %*% Matrix::t(Matrix::colSums(A^2)))

	Lambda = Obs - Exp
	
	a = apply(A, 2, max)
  ones = array(1, dim = dim(Nu)[1])
  
	logPvals = Lambda^2 / (2*(Nu + (ones %*% t(a))*Lambda/3))
	
	rownames(logPvals) = colnames(X)
	colnames(logPvals) = colnames(A)

	logPvals[Lambda  < 0] = 0
	logPvals[is.na(logPvals)] = 0
	
	logPvals = logPvals / log(10)
	
		
	return(logPvals)
}

plot.atACTIONet.motif.view <- function(ACTIONet.out, motif.enrichment, split = TRUE) {
	if(split == TRUE) {
		parts = stringr::str_split_fixed(rownames(motif.enrichment),  "_", 5)
		motif.name = parts[, 3]
	} else {
		motif.name = rownames(motif.enrichment)
	}
	top.motifs = apply(motif.enrichment, 2, function(x) motif.name[order(x, decreasing = TRUE)[1:5]])
	
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
