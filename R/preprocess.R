
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
  profile <- rowSums(profile_mat)

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

  inserts <- c(
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(start(fragments), start(fragments)), RG = mcols(fragments)[,by]),
    GRanges(seqnames = seqnames(fragments), ranges = IRanges(end(fragments), end(fragments)), RG = mcols(fragments)[,by])
  )
  overlapDF <- DataFrame(GenomicRanges::findOverlaps(query, inserts, ignore.strand = TRUE, maxgap=-1L, minoverlap=0L, type = "any"))
  overlapDF$name <- mcols(inserts)[overlapDF[, 2], by]
  overlapTDF <- transform(overlapDF, id = match(name, unique(name)))

  inPeaks <- table(overlapDF$name)
  total <- table(mcols(inserts)[, by])
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


import.and.filter.frags <- function(input_path, txdb, sample_name = NULL, min_frags_per_cell = NULL, min_TSS_per_cell = NULL, min_cells_per_frag = NULL, keep_filtered = FALSE, save_frags = FALSE, plot_TSS = FALSE, save_qc = FALSE, save_dir = NULL){
  require(ggplot2)
  if(is.null(save_dir)){
    save_dir = getwd()
  }
  if(save_frags | plot_TSS | save_qc){
    if(is.null(sample_name)){
      stop('sample_name needed to save requested output')
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
  fragments <- data.frame(readr::read_tsv(file_fragments, col_names=FALSE))

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

  ## Compute TSS Profile
  feature <- txdb %>% transcripts(.) %>% resize(., width = 1, fix = "start") %>% unique
  tssProfile <- insertionProfileSingles(feature = feature, fragments = fragments,
                                        getInsertions = TRUE, batchSize = 1000)
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
    saveRDS(fragments, frag_dir)
  }
  out <- list(fragments = fragments, tssProfile = tssProfile[c('df', 'dfall')])
  invisible(gc())
  return(out)
}

frag.counts.to.ace <- function(mat, features, binarize = TRUE, nFeatures = NULL){

  rownames(mat) <- seq_len(nrow(mat))
  names(windows) <- rownames(mat)

  # message("Making ACE Object...")
  sce <- SingleCellExperiment(
    assays = SimpleList(counts = mat),
    rowRanges = windows
  )
  rownames(sce) <- paste(seqnames(sce),start(sce),end(sce), sep = "_")

  if(binarize){
    # message(paste0("Binarizing matrix..."))
    sce <- add.binarized.counts(sce)
  }

  if(!is.null(nFeatures)){
    # message(paste0("Getting top ", nFeatures, " features..."))
    if(binarize){
      sce <- sce[head(order(Matrix::rowSums(assays(sce_pre)[["bin_counts"]]), decreasing = TRUE), nFeatures),]
    } else{
      sce <- sce[head(order(Matrix::rowSums(assays(sce_pre)[["counts"]]), decreasing = TRUE), nFeatures),]
    }
  }

frags.to.ace <- function(fragments, features, by = "RG", binarize = TRUE, nFeatures = NULL){

  counts <- countInsertions(features, fragments, by = by)[[1]]
  ace = frag.counts.to.ace(counts, features, binarize = binarize, nFeatures = nFeatures)
  return(ace)
}