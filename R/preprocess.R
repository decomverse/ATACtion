
# Adapted from Granja, et al. (2019)
# Original at https://github.com/GreenleafLab/MPAL-Single-Cell-2019

insertionProfileSingles_helper <- function(feature, fragments, by = "RG", getInsertions = TRUE, fix = "center", flank = 2000, norm = 100, smooth = 51, range = 100, batchSize = 100){
  #Convert To Insertion Sites
  if(getInsertions){
    insertions <- c(
      GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), RG = mcols(fragments)[,by]),
      GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), RG = mcols(fragments)[,by])
    )
    by <- "RG"
  }else{
    insertions <- fragments
  }
  remove(fragments)
  gc()

  #center the feature
  center <- unique(resize(feature, width = 1, fix = fix, ignore.strand = FALSE))

  #get overlaps between the feature and insertions only up to flank bp
  overlap <- DataFrame(findOverlaps(query = center, subject = insertions, maxgap = flank, ignore.strand = TRUE))
  overlap$strand <- strand(center)[overlap[,1]]
  overlap$name <- mcols(insertions)[overlap[,2],by]
  overlap <- transform(overlap, id=match(name, unique(name)))
  ids <- length(unique(overlap$name))

  #distance
  overlap$dist <- NA
  minus <- which(as.logical(overlap$strand == "-"))
  other <- which(as.logical(overlap$strand != "-"))
  overlap$dist[minus] <- start(center[overlap[minus,1]]) - start(insertions[overlap[minus,2]])
  overlap$dist[other] <- start(insertions[overlap[other,2]]) - start(center[overlap[other,1]])

  #Insertion Mat
  profile_mat <- tabulate2dCpp(x1 = overlap$id, y1 = overlap$dist, xmin = 1, xmax = ids, ymin = -flank, ymax = flank)
  colnames(profile_mat) <- unique(overlap$name)
  profile <- ACTIONet::fast_row_sums(profile_mat)

  #normalize
  profile_mat_norm <- apply(profile_mat, 2, function(x) x/max(mean(x[c(1:norm,(flank*2-norm+1):(flank*2+1))]), 0.5)) #Handles low depth cells
  profile_norm <- profile/mean(profile[c(1:norm,(flank*2-norm+1):(flank*2+1))])

  #smooth
  profile_mat_norm_smooth <- apply(profile_mat_norm, 2, function(x) zoo::rollmean(x, smooth, fill = 1))
  profile_norm_smooth <- zoo::rollmean(profile_norm, smooth, fill = 1)

  #enrichment
  max_finite <- function(x){
    suppressWarnings(max(x[is.finite(x)], na.rm=TRUE))
  }
  e_mat <- apply(profile_mat_norm_smooth, 2, function(x) max_finite(x[(flank-range):(flank+range)]))
  names(e_mat) <- colnames(profile_mat_norm_smooth)
  e <- max_finite(profile_norm_smooth[(flank-range):(flank+range)])

  #Summary
  df_mat <- data.frame(
    enrichment = e_mat,
    insertions = as.vector(table(mcols(insertions)[,by])[names(e_mat)]),
    insertionsWindow = as.vector(table(overlap$name)[names(e_mat)])
  )
  df_sum <- data.frame(bp = (-flank):flank, profile = profile, norm_profile = profile_norm, smooth_norm_profile = profile_norm_smooth, enrichment = e)
  rownames(df_sum) <-  NULL

  return(list(df = df_sum, dfall = df_mat, profileMat = profile_mat_norm, profileMatSmooth = profile_mat_norm_smooth))
}

insertionProfileSingles <- function(feature, fragments, sample_name = NULL, by = "RG", getInsertions = TRUE, fix = "center", flank = 2000, norm = 100, smooth = 51, range = 100, batchSize = 100){

  uniqueTags <- as.character(unique(mcols(fragments)[,by]))
  splitTags <- split(uniqueTags, ceiling(seq_along(uniqueTags)/batchSize))

  pb <- txtProgressBar(min = 0, max = 100, initial = 0, style = 3)
  batchTSS <- lapply(seq_along(splitTags), function(x){
    setTxtProgressBar(pb, round(x * 100/length(splitTags), 0))
    profilex <- insertionProfileSingles_helper(
      feature=feature,
      fragments=fragments[which(mcols(fragments)[,by] %in% splitTags[[x]])],
      by = by,
      getInsertions = getInsertions,
      fix = fix,
      flank = flank,
      norm = norm,
      smooth = smooth,
      range = range
    )

    return(profilex)
  })
  df <- lapply(batchTSS, function(x) x$df) %>% Reduce("rbind",.)
  dfall <- lapply(batchTSS, function(x) x$dfall) %>% Reduce("rbind",.)
  profileMat <- lapply(batchTSS, function(x) x$profileMat) %>% Reduce("cbind",.)
  profileMatSmooth <- lapply(batchTSS, function(x) x$profileMatSmooth) %>% Reduce("cbind",.)
  return(list(df = df, dfall = dfall, profileMat = profileMat, profileMatSmooth = profileMatSmooth))
}

countInsertions <- function(query, fragments, by = "RG"){
	require(GenomicRanges)
  inserts <- c(
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), RG = mcols(fragments)[,by]),
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), RG = mcols(fragments)[,by])
  )
  suppressWarnings({ overlapDF <- DataFrame(GenomicRanges::findOverlaps(query, inserts, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any")) })
  overlapDF$name <- mcols(inserts)[overlapDF[, 2], "RG"]
  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))

  inPeaks <- table(overlapDF$name)
  total <- table(mcols(inserts)[, "RG"])
  total <- total[names(inPeaks)]
  frip <- inPeaks / total

  sparseM <- Matrix::sparseMatrix(
    i = overlapTDF[, 1],
    j = overlapTDF[, 4],
    x = rep(1, nrow(overlapTDF)),
    dims = c(length(query), length(unique(overlapDF$name))))
  colnames(sparseM) <- unique(overlapDF$name)
  total <- total[colnames(sparseM)]
  frip <- frip[colnames(sparseM)]
  out <- list(counts = sparseM, frip = frip, total = total)
  invisible(gc())
  return(out)
}


import.and.filter.frags <- function(input_path, reference_genome = "hg38", sample_name = NULL, min_frags_per_cell = 1000, min_TSS_per_cell = 8, min_cells_per_frag = 100, keep_filtered = FALSE, save_frags = FALSE, plot_TSS = FALSE, save_qc = FALSE, save_dir = NULL, header = F){
	require(viridis)
	  if(reference_genome == 'mm10') {
		library(TxDb.Mmusculus.UCSC.mm10.knownGene)
		txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
	  }
	  else if (reference_genome == 'mm9') {
		library(TxDb.Mmusculus.UCSC.mm9.knownGene)
		txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
	  }
	  else if(reference_genome == 'grch37' || reference_genome == 'hg19') {
		library(TxDb.Hsapiens.UCSC.hg19.knownGene)
		txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
	  }
	  else if(reference_genome == 'grch38' || reference_genome == 'hg38') {
		library(TxDb.Hsapiens.UCSC.hg38.knownGene)
		txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
	  }
	  else {
		  R.utils::printf('Unknown reference_genome %s\n', reference_genome);
		  return()
	  }	
	
  require(ggplot2)
  if(is.null(save_dir)){
    save_dir = getwd()
  }
  if(save_frags | plot_TSS | save_qc){
    if(is.null(sample_name)){
      stop('sample_name needed to save output')
    }
    if(plot_TSS | save_qc){
      md_dir = file.path(save_dir, 'qc_metadata')
      dir.create(md_dir, recursive = T, showWarnings = F)
    }
  }

  minFrags <- min_cells_per_frag
  filterFrags <- min_frags_per_cell
  filterTSS <- min_TSS_per_cell
  file_fragments <- input_path
  name <- sample_name

  ## Read Fragment Files
  message("Reading in fragment files...")
  fragments <- data.frame(readr::read_tsv(file_fragments, col_names = header, col_types = c(chr = "c", start = "d", end = "d", RG = "c", N ="d"), progress = TRUE))

  fragments <- GRanges(
    seqnames = fragments[,1],
    IRanges(fragments[,2]+1, fragments[,3]),
    RG = fragments[,4],
    N = fragments[,5]
  )

  message("Filtering Lowly Represented Cells...")
  tabRG <- table(fragments$RG)
  keep <- names(tabRG)[which(tabRG >= minFrags)]
  fragments <- fragments[fragments$RG %in% keep,]
  fragments <- sort(sortSeqlevels(fragments))

  genome(fragments) = reference_genome

  ## Compute TSS Profile
  feature <- txdb %>% transcripts(.) %>% resize(., width = 1, fix = "start") %>% unique
  suppressWarnings({ tssProfile <- insertionProfileSingles(feature = feature, fragments = fragments,
                                        getInsertions = TRUE, batchSize = 1000) })

  tssSingles <- tssProfile$dfall
  tssSingles$uniqueFrags <- 0
  tssSingles[names(tabRG),"uniqueFrags"] <- tabRG
  tssSingles$cellCall <- 0
  tssSingles$cellCall[tssSingles$uniqueFrags >= filterFrags & tssSingles$enrichment >= filterTSS] <- 1
  tssSingles <- tssSingles[complete.cases(tssSingles),]
  tssProfile$dfall <- tssSingles
  fragments$qc_pass <- 0
  mcols(fragments)[,"qc_pass"][mcols(fragments)[,"RG"] %in% rownames(tssSingles)[tssSingles$cellCall == 1]] <- 1

  ## Plot Stats
  if(plot_TSS){
    tssSingles <- tssSingles[complete.cases(tssSingles),]
    nPass  <- sum(tssSingles$cellCall==1)
    nTotal <- sum(tssSingles$uniqueFrags >= filterFrags)

    plot_dir = file.path(md_dir, paste0(name, '_tss.pdf') )
    p <- ggplot(tssSingles[tssSingles$uniqueFrags > 500,], aes(x = log10(uniqueFrags), y = enrichment)) +
      # geom_bin2d(bins = 100) +
      geom_hex(bins = 100) +
      theme_bw() + scale_fill_viridis() +
      xlab("log10 Unique Fragments") +
      ylab("TSS Enrichment") +
      geom_hline(yintercept = filterTSS, lty = "dashed") +
      geom_vline(xintercept = log10(filterFrags), lty = "dashed") +
      ggtitle(sprintf("Pass Rate : %s of %s (%s)", nPass, nTotal, round(100*nPass/nTotal,2)))
    ggsave(plot_dir, plot = p, device = 'pdf')
  }


  ## Save Output
  if(save_qc){
    tss_dir = file.path(md_dir, paste0(name, '_filter-cells.txt') )
    write.table(tssSingles, tss_dir)
  }
  ## Filter fragments
  if(!keep_filtered){
    fragments <- fragments[fragments$qc_pass==1]
  }

  if(save_frags){
    dir.create(save_dir, recursive = T, showWarnings = F)
    frag_dir = file.path(save_dir, paste0(name, '_fragments.rds') )
    readr::write_rds(fragments, frag_dir)
  }
  out <- list(fragments = fragments, tssProfile = tssProfile[c('df', 'dfall')])
  invisible(gc())
  return(out)
}

frag.counts.to.ace <- function(mat, features, binarize = TRUE) {

	mat = as(mat, "sparseMatrix")
	ace = ACTIONetExperiment(assays = list(counts = mat))
	if(binarize == TRUE) {
        mat@x = rep(1, length(mat@x))
        assays(ace)[["bin_counts"]] = mat	
	}
	rowRanges(ace) = features
	
	rnames = paste(as.character(seqnames(features)), start(features), end(features), sep = "_")
	rownames(ace) = rnames
	
	return(ace)
}

frags.to.ace <- function(fragments, features, binarize = TRUE, reference.genome = "hg38", by = "RG"){
  genome(features) = reference.genome

  counts <- countInsertions(features, fragments, by = by)[[1]]
  
  ace = frag.counts.to.ace(counts, features, binarize = binarize)
  return(ace)
}

bin_genome <- function(reference_genome = "hg38", bin_size = 200) {
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
		
	
	chromSizes <- GRanges(names(seqlengths(genome)), IRanges(1, seqlengths(genome)))
	chromSizes <- GenomeInfoDb::keepStandardChromosomes(chromSizes, pruning.mode = "coarse")
	GR <- unlist(tile(chromSizes, width = bin_size)) 
	return(GR)	
}


construct_optimal_ATAC_ace <- function(frags, flank = 2500, bin_size = 250, enrichment.threshold = 1, frag_by = "RG") {
	reference_genome = tolower(genome(frags)[[1]])
	
	  if(reference_genome == 'mm10') {
		library(TxDb.Mmusculus.UCSC.mm10.knownGene)
		txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene
	  }
	  else if (reference_genome == 'mm9') {
		library(TxDb.Mmusculus.UCSC.mm9.knownGene)
		txdb <- TxDb.Mmusculus.UCSC.mm9.knownGene
	  }
	  else if(reference_genome == 'grch37' || reference_genome == 'hg19') {
		library(TxDb.Hsapiens.UCSC.hg19.knownGene)
		txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
	  }
	  else if(reference_genome == 'grch38' || reference_genome == 'hg38') {
		library(TxDb.Hsapiens.UCSC.hg38.knownGene)
		txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
	  }
	  else {
		  R.utils::printf('Unknown reference_genome %s\n', reference_genome);
		  return()
	  }
	
	# Focus on the proximity of genes first
	message("Aggregating fragments within gene proximities ... ")
	
	suppressMessages( {features.GR <- txdb %>% genes(.) %>% unique})
	start(features.GR) = start(features.GR) - flank
	system.time( {insert.out = countInsertions(query = features.GR, frags, by = frag_by)} )

	# Run ACTIONet
	message("Run ACTIONet ... ")

	gene.ace = ACTIONetExperiment(assays = list(counts = insert.out$counts))
	rowRanges(gene.ace) = features.GR

	gene.ace = reduce.ace(gene.ace)
	gene.ace = run.ACTIONet(gene.ace, data_slot = "logcounts", reduction_slot = "ACTION")


	# Used fixed-binning and then select specific bins
	message("Aggregating fragments within fixed bins ... ")

	bins.GR = bin_genome(reference_genome, bin_size)
	system.time( {insert.out.bins = countInsertions(query = bins.GR, frags, by = frag_by)} )


	# Compute specificity of bins based on archetypes from gene-inferred ACTIONet
	message("Computing the most informative bins/peaks ... ")

	S = insert.out.bins$counts
	S@x = rep(1, length(S@x))
	H_unified = as.matrix(Matrix::t(colMaps(gene.ace)[["H_unified"]]))
	idx = match(colnames(S), colnames(gene.ace))
	H_unified = H_unified[, idx]
	system.time( {specificity.out = compute_archetype_feature_specificity(S, H_unified) })

	spec = specificity.out$upper_significance
	selected.bins = sort(unique(unlist(apply(spec, 2, function(x) which(x > enrichment.threshold)))))


	peak.counts = insert.out.bins$counts[selected.bins, ]
	ace = ACTIONetExperiment(assays = list(counts = peak.counts))
	rowRanges(ace) = bins.GR[selected.bins]
	genome(rowRanges(ace)) = reference_genome
	ace = ace[, match(colnames(gene.ace), colnames(ace))]

	colMaps(ace)$geneAce_H = Matrix::t(H_unified)
	rowMaps(ace)$geneAce_peak_specificity = spec[selected.bins, ]
		
	# Doublet detection
	message("Computing cell QC metrics ... ")
	ace$centrality_score = gene.ace$node_centrality
	

	scores = rowMaps(ace)$geneAce_peak_specificity
	scores = apply(scores, 2, function(x)  exp(scale(x)))
	
	associations = counts(ace)
	associations@x = rep(1, length(associations@x))
	system.time( {enrichment.out = assess_enrichment(scores, associations, 1)} )


	X = apply(enrichment.out$logPvals, 2, function(x) x / sum(x))
	h = apply(X, 1, ineq::entropy)
			
	ace$doublet_score = scale(-h)
	
	out = list(atac.ace = ace, gene.ace = gene.ace, selected_peaks = bins.GR[selected.bins])

	return(out)
}
