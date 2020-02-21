archCicerify <- function(ACTIONet.out, sce, accessibility.slot = "signature") {
  require(cicero)

  archs = ACTIONet.out$archetype.accessibility@assays[[accessibility.slot]]
  GR = rowRanges(sce)
  gn = genome(rowRanges(sce))[[1]]
  
  if(is.na(seqlengths(rowRanges(sce))[[1]])) {
      if(is.na(gn)) {
        R.utils::printf("Error: No genome has been set for the sce object.")
      } else {
        fname = system.file("extdata", "chrLen", sprintf('%s_chrLen.txt', gn) , package = "ATACtiondb")
        if(! file.exists(fname) ) {
          R.utils::printf("Error: Genome: %s is not supported. Please pass chromosome length as part of the seqlength() of sce object.")
        } else {
          chr.table = read.csv(fname, sep = '\t')
        }
      }
  } else {
	  SL = seqlengths(rowRanges(sce))
	  chr.table = data.frame(chr = names(SL), Len = as.numeric(SL), stringsAsFactors = FALSE)
  }


  
  peak.labels = paste(GR@seqnames, GR@ranges@start, GR@ranges@start+GR@ranges@width, sep = '_')
  arch.labels = sapply(1:ncol(archs), function(x) sprintf('Arch%03d', x))
  
  archs = as(archs, 'dgTMatrix')
  df = data.frame(Peak = peak.labels[archs@i+1], Cell = arch.labels[archs@j+1], Count = as.numeric(archs@x), stringsAsFactors = FALSE)
  
  arch.CDS = make_atac_cds(df, binarize = FALSE)
  arch.CDS@assayData[['exprs']] = as.matrix(arch.CDS@assayData[['exprs']])

	## Annotate with gene promoters for future reference (needed in build_gene_activity_matrix())
  ds.name = sprintf('refGene_%s', gn)

  path = system.file(package="ATACtiondb", "data", sprintf('refGene_%s.rda', gn))
  load(path)  

  promoter.bed = data.frame(chr = as.character(refGene@seqnames), start = as.numeric(refGene@ranges@start), end = as.numeric(refGene@ranges@start + refGene@ranges@width - 1), gene = refGene$Gene, stringsAsFactors = FALSE)
  arch.CDS = annotate_cds_by_site(arch.CDS, promoter.bed)

  
  system.time( {conns <- run_cicero(arch.CDS, chr.table)} )
  
  
  CCAN_assigns <- generate_ccans(conns)

  
  out = list(connections = conns, modules = CCAN_assigns, arch.CDS = arch.CDS)
  
  return(out)  
}


inferExpr.cicero <- function(arch.CDS, conns) {
	require(cicero)
	
	unnorm_ga <- build_gene_activity_matrix(arch.CDS, conns)

	num_genes <- pData(arch.CDS)$num_genes_expressed
	names(num_genes) <- row.names(pData(arch.CDS))
	cicero_gene_activities <- as.matrix(normalize_gene_activities(unnorm_ga, num_genes))

	return(cicero_gene_activities)
}

  

cicero.to.GenomicInteractions <- function(sce, conns, coaccess_cutoff = 0.25) {
  GR = rowRanges(sce)
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


GenomicInteractions.to.cicero <- function(sce, GI) {
  GR = rowRanges(sce)
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


conn.to.peakConnectivity.mat <- function(sce, conns, coaccess_cutoff = 0.25) {
	GR = rowRanges(sce)
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

GI.to.peakConnectivity.mat <- function(sce, GI, maxgap = 100) {
	matches <- GenomicRanges::findOverlaps(rowRanges(sce), GI@regions, select = "all", maxgap)
	Ind = Matrix::sparseMatrix(i = S4Vectors::queryHits(matches), j = S4Vectors::subjectHits(matches), x = 1, dims =  c(nrow(sce), length(GI@regions)))

	A = Matrix::sparseMatrix(i = GI@anchor1, j = GI@anchor2, x = 1, dims = c(length(GI@regions), length(GI@regions)))
	A = A + t(A)

	connectivity.map = Ind %*% A %*% t(Ind)
	connectivity.map@x[connectivity.map@x > 1] = 1
	diag(connectivity.map) = 0

	return(connectivity.map)
}	


extend.cisConnectome.usingConns <-function(sce, conns, coaccess_cutoff=0.25) {
  if( !("cisConnectome" %in% names(metadata(sce))) ) {
	  print("cisConnectome slot does not exist in sce metadata")
	  return()
  }

  connectivity.map = conn.to.peakConnectivity.mat(sce, conns, coaccess_cutoff)
  cisConnectome.mat = metadata(sce)[["cisConnectome"]]
  
  extended.cisConnectome = connectivity.map %*% cisConnectome.mat
  
  return(extended.cisConnectome)
}  

extend.cisConnectome.usingGI <-function(sce, GI, max_gap = 100) {
  if( !("cisConnectome" %in% names(metadata(sce))) ) {
	  print("cisConnectome slot does not exist in sce metadata")
	  return()
  }

  connectivity.map = GI.to.peakConnectivity.mat(sce, GI, max_gap)
  cisConnectome.mat = metadata(sce)[["cisConnectome"]]
  
  extended.cisConnectome = connectivity.map %*% cisConnectome.mat
  
  return(extended.cisConnectome)
}  

