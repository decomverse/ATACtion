compute_ATACtion_archSE <- function(ace, arch_slot = "unified_feature_specificity") {
	GR = SummarizedExperiment::rowRanges(ace)
	archs = rowMaps(ace)[[arch_slot]]
	SE = SummarizedExperiment(assays = list(peaks = archs), SummarizedExperiment::rowRanges = GR)
	
	return(SE)
}

run_ATACtion_cicerify <- function(SE, assay_slot = "peaks") {
	library(cicero)
  gn = genome(SummarizedExperiment::rowRanges(SE))[[1]]
  if(is.na(seqlengths(SummarizedExperiment::rowRanges(SE))[[1]])) {
      if(is.na(gn)) {
        R.utils::printf("Error: No genome has been set for the SE object.")
      } else {
        fname = system.file("extdata", "chrLen", sprintf('%s_chrLen.txt', gn) , package = "ATACtionDB")
        if(! file.exists(fname) ) {
          R.utils::printf("Error: Genome: %s is not supported. Please pass chromosome length as part of the seqlength() of SE object.")
        } else {
          chr.table = read.csv(fname, sep = '\t')
        }
      }
  } else {
	  SL = seqlengths(SummarizedExperiment::rowRanges(SE))
	  chr.table = data.frame(chr = names(SL), Len = as.numeric(SL), stringsAsFactors = FALSE)
  }

	peak.labels = rownames(SE)
	arch.labels = colnames(SE)
	archs = as(assays(SE)[[assay_slot]], 'dgTMatrix')
	
	df = data.frame(Peak = peak.labels[archs@i+1], Cell = arch.labels[archs@j+1], Count = as.numeric(archs@x), stringsAsFactors = FALSE)
	
	arch.CDS = make_atac_cds(df, binarize = FALSE)
	perm = match(paste("A", 1:ncol(arch.CDS), sep = ""), colnames(arch.CDS))
	arch.CDS = arch.CDS[, perm]
	
	## Annotate with gene promoters for future reference (needed in build_gene_activity_matrix())
	path = system.file(package="ATACtionDB", "data", sprintf('refGene_%s.rda', gn))
	load(path)  
	
	
	promoter.bed = data.frame(chr = as.character(refGene@seqnames), start = as.numeric(refGene@ranges@start), end = as.numeric(refGene@ranges@start + refGene@ranges@width - 1), gene = refGene$Gene, stringsAsFactors = FALSE)
	arch.CDS = annotate_cds_by_site(arch.CDS, promoter.bed)
	
	system.time( {conns <- run_cicero(arch.CDS, chr.table)} )
	
	out = list(conns = conns, arch.CDS = arch.CDS)
	return(out)  
}

compute_cicero_ccans <- function(conns) {
	system.time( {CCAN_assigns <- generate_ccans(conns)} )

	return(CCAN_assigns)	
}


compute_cicero_gene_expression <- function(arch.CDS, conns) {
	require(cicero)
	
	
	unnorm_ga <- build_gene_activity_matrix(arch.CDS, conns)

	num_genes <- pData(arch.CDS)$num_genes_expressed
	names(num_genes) <- row.names(pData(arch.CDS))
	cicero_gene_activities <- as.matrix(normalize_gene_activities(unnorm_ga, num_genes))

	return(cicero_gene_activities)
}  
  
cicero.to.GenomicInteractions <- function(ace, conns, coaccess_cutoff = 0.25) {
  GR = SummarizedExperiment::rowRanges(ace)
  values(GR) = c()
  
  # peak.labels = paste(GR@seqnames, GR@ranges@start, GR@ranges@start+GR@ranges@width-1, sep = '_')  
  peak.labels = paste(GR@seqnames, GR@ranges@start, GR@ranges@start+GR@ranges@width, sep = '_')
  
  ii = match(conns$Peak1, peak.labels)
  jj = match(conns$Peak2, peak.labels)
  vv = conns$coaccess
  mask = !is.na(vv) & vv > coaccess_cutoff
  

  anchor.one = GR[ii[mask]]
  anchor.two = GR[jj[mask]]
  
  GR.interactions = GenomicInteractions::GenomicInteractions(anchor.one, anchor.two, coaccess=vv[mask])
  
  return(GR.interactions)  
}


GenomicInteractions.to.cicero <- function(ace, GI) {
  GR = SummarizedExperiment::rowRanges(ace)
  values(GR) = c()
  
  # peak.labels = paste(GR@seqnames, GR@ranges@start, GR@ranges@start+GR@ranges@width-1, sep = '_')  
  peak.labels = paste(GR@seqnames, GR@ranges@start, GR@ranges@start+GR@ranges@width, sep = '_')


  GI.conns = data.frame(Peak1 = peak.labels[GI@anchor1], Peak2 = peak.labels[GI@anchor2])
  GI.conns = cbind(GI.conns, as.data.frame(GI@elementMetadata))
  
  
  ii = match(conns$Peak1, peak.labels)
  jj = match(conns$Peak2, peak.labels)
  vv = conns$coaccess
  mask = !is.na(vv) & vv > coaccess_cutoff
  

  anchor.one = GR[ii[mask]]
  anchor.two = GR[jj[mask]]
  
  GR.interactions = GenomicInteractions::GenomicInteractions(anchor.one, anchor.two, coaccess=vv[mask])
  
  return(GR.interactions)  
}


conn.to.peakConnectivity.mat <- function(ace, conns, coaccess_cutoff = 0.25) {
	GR = SummarizedExperiment::rowRanges(ace)
	values(GR) = c()

	# peak.labels = paste(GR@seqnames, GR@ranges@start, GR@ranges@start+GR@ranges@width-1, sep = '_')  
	peak.labels = paste(GR@seqnames, GR@ranges@start, GR@ranges@start+GR@ranges@width, sep = '_')

	ii = match(conns$Peak1, peak.labels)
	jj = match(conns$Peak2, peak.labels)
	vv = conns$coaccess
	mask = !is.na(vv) & vv > coaccess_cutoff
	connectivity.map = Matrix::sparseMatrix(i = ii[mask], j = jj[mask], x = vv[mask], dims = c(length(GR), length(GR)))  	
	
	return(connectivity.map)
}

GI.to.peakConnectivity.mat <- function(ace, GI, maxgap = 100) {
	matches <- GenomicRanges::findOverlaps(SummarizedExperiment::rowRanges(ace), GI@regions, select = "all", maxgap)
	Ind = Matrix::sparseMatrix(i = S4Vectors::queryHits(matches), j = S4Vectors::subjectHits(matches), x = 1, dims =  c(nrow(ace), length(GI@regions)))

	A = Matrix::sparseMatrix(i = GI@anchor1, j = GI@anchor2, x = 1, dims = c(length(GI@regions), length(GI@regions)))
	A = A + t(A)

	connectivity.map = Ind %*% A %*% t(Ind)
	connectivity.map@x[connectivity.map@x > 1] = 1
	diag(connectivity.map) = 0

	return(connectivity.map)
}	


extend.cisConnectome.usingConns <-function(ace, conns, coaccess_cutoff=0.25) {
  if( !("cisConnectome" %in% names(metadata(ace))) ) {
	  print("cisConnectome slot does not exist in ace metadata")
	  return()
  }

  connectivity.map = conn.to.peakConnectivity.mat(ace, conns, coaccess_cutoff)
  cisConnectome.mat = rowMaps(ace)[["cisConnectome"]]
  
  extended.cisConnectome = connectivity.map %*% cisConnectome.mat
  
  return(extended.cisConnectome)
}  

extend.cisConnectome.usingGI <-function(ace, GI, max_gap = 100) {
  if( !("cisConnectome" %in% names(metadata(ace))) ) {
	  print("cisConnectome slot does not exist in ace metadata")
	  return()
  }

  connectivity.map = GI.to.peakConnectivity.mat(ace, GI, max_gap)
  cisConnectome.mat = rowMaps(ace)[["cisConnectome"]]
  
  extended.cisConnectome = connectivity.map %*% cisConnectome.mat
  
  return(extended.cisConnectome)
}  

