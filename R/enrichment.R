impute_gene_expression_using_ATACtion <- function(ace, genes, enrichment_name = "genes", slot_name = "unified") {
	if(! ("enrichments" %in% names(metadata(ace))) ) {
		warning(sprintf("Couldn't find enrichment in the metadata(ace).\n"))
		return()
	}
	if( !(slot_name %in% names(metadata(ace)$enrichments)) ) {
		warning(sprintf("Couldn't find %s in enrichment list.\n", slot_name))
		return()
	} 
	
	enrichments = metadata(ace)$enrichments[[slot_name]]
	if( !(enrichment_name %in% names(enrichments)) ) {
		warning(sprintf("Couldn't find %s in enrichment list of %s.\n", enrichment_name, slot_name))
		return()
	} 	
	arch_gene_expression = enrichments[[enrichment_name]]
	
	common.genes = sort(unique(intersect(genes, rownames(arch_gene_expression))))
	sub_arch_gene_expression = arch_gene_expression[common.genes, ]
	
	H = colMaps(ace)[[sprintf("H_%s", slot_name)]]
	gene.expression = sub_arch_gene_expression %*% Matrix::t(H)
	
	rownames(gene.expression) = common.genes
	colnames(gene.expression) = colnames(ace)
	
	return(gene.expression)
}

impute_TF_activity_using_ATACtion <- function(ace, TFs, enrichment_name = "TFs", slot_name = "unified") {
	if(! ("enrichments" %in% names(metadata(ace))) ) {
		warning(sprintf("Couldn't find enrichment in the metadata(ace).\n"))
		return()
	}
	if( !(slot_name %in% names(metadata(ace)$enrichments)) ) {
		warning(sprintf("Couldn't find %s in enrichment list.\n", slot_name))
		return()
	} 
	
	enrichments = metadata(ace)$enrichments[[slot_name]]
	if( !(enrichment_name %in% names(enrichments)) ) {
		warning(sprintf("Couldn't find %s in enrichment list of %s.\n", enrichment_name, slot_name))
		return()
	} 	
	arch_TF_activity = enrichments[[enrichment_name]]
	
	common.TFs = sort(unique(intersect(TFs, rownames(arch_TF_activity))))
	sub_arch_TF_activity = arch_TF_activity[common.TFs, ]
	
	H = colMaps(ace)[[sprintf("H_%s", slot_name)]]
	TF.activity = sub_arch_TF_activity %*% Matrix::t(H)
	
	rownames(TF.activity) = common.TFs
	colnames(TF.activity) = colnames(ace)
	
	return(TF.activity)
}

visualize_marker_genes_over_ATACtion <- function(ace, genes, enrichment_name = "genes", slot_name = "unified") {
	imputed_gene_expression = impute_gene_expression_using_ATACtion(ace, genes = genes, enrichment_name = enrichment_name, slot_name = slot_name)

	imputed_genes = rownames(imputed_gene_expression)
	
	sapply(imputed_genes, function(gene) {
		plot.ACTIONet.gradient(ace, imputed_gene_expression[gene, ], title = gene, transparency.attr = ace$node_centrality, node.size = 0.01)
	})
}


get_ATACtion_enrichment <- function(ace, enrichment_name = "genes", slot_name = "unified") {
	if(! ("enrichments" %in% names(metadata(ace))) ) {
		warning(sprintf("Couldn't find enrichment in the metadata(ace).\n"))
		return()
	}
	if( !(slot_name %in% names(metadata(ace)$enrichments)) ) {
		warning(sprintf("Couldn't find %s in enrichment list.\n", slot_name))
		return()
	} 
	
	enrichments = metadata(ace)$enrichments[[slot_name]]
	if( !(enrichment_name %in% names(enrichments)) ) {
		warning(sprintf("Couldn't find %s in enrichment list of %s.\n", enrichment_name, slot_name))
		return()
	} 	
	enrichment = enrichments[[enrichment_name]]
	
	return(enrichment)
}


compute_gene_enrichment_from_ATACtion <- function(ace, min_peaks = 5, slot_name = "unified", gene_slot_name = "genes", thread_no = 0, cisConnectome = "proximal_cisConnectome", min.score = 0, min.association = 0) {
	if(! (cisConnectome %in% names(rowMaps(ace))) ) {
		warning(sprintf("Couldn't find %s in the rowMaps(ace). Please call Add_proximal_peak_gene_interactions_to_ATACtion() or Add_physical_peak_gene_interactions_to_ATACtion() prior to calling gene_enrichment_from_ATACtion()", cisConnectome))
		return(ace)
	}	

	associations = as(rowMaps(ace)[[cisConnectome]], 'dgCMatrix')



	scores = rowMaps(ace)[[sprintf("%s_feature_specificity", slot_name)]]
	if(max(scores) > 100) {
		scores = log1p(scores)
	}
	
	x = apply(scores, 1, max)
	y = Matrix::colSums(Matrix::t(associations))

	mask = (x > min.score) & (y > min.association)
	associations = associations[mask, ]
	scores = scores[mask, ]
	
	enrichment.out = assess_enrichment(scores, associations, thread_no)
	enrichment.mat = enrichment.out$logPvals
	rownames(enrichment.mat) = colnames(associations)

	if(! ("enrichments" %in% names(metadata(ace))) ) {
		metadata(ace)[["enrichments"]] = list()
	}
	if( !(slot_name %in% names(metadata(ace)$enrichments)) ) {
		L = list(genes = enrichment.mat)
		names(L) = gene_slot_name
		metadata(ace)$enrichments[[slot_name]] = L
	} else {
		L = metadata(ace)$enrichments[[slot_name]]
		L[[gene_slot_name]] = enrichment.mat
	}
	metadata(ace)$enrichments[[slot_name]] = L
	
	return(ace)
}




compute_geneset_enrichment_from_ATACtion <- function(ace, genesets, min_peaks = 5, slot_name = "unified", geneset_slot_name = "genesets", thread_no = 0, cisConnectome = "proximal_cisConnectome", min.score = 0, min.association = 0) {
	if(! (cisConnectome %in% names(rowMaps(ace))) ) {
		warning(sprintf("Couldn't find %s in the rowMaps(ace). Please call Add_proximal_peak_gene_interactions_to_ATACtion() or Add_physical_peak_gene_interactions_to_ATACtion() prior to calling gene_enrichment_from_ATACtion()", cisConnectome))
		return(ace)
	}	

	peak.gene.associations = as(rowMaps(ace)[[cisConnectome]], 'dgCMatrix')
	mask = (Matrix::colSums(peak.gene.associations) > min_peaks)
	peak.gene.associations = peak.gene.associations[, mask]
	
	
	gene.to.geneset.associations = do.call(cbind, lapply(genesets, function(geneset) {
		common.genes = intersect(geneset, colnames(peak.gene.associations))
		ii = match(common.genes, colnames(peak.gene.associations))
		v = as(sparseVector(x = 1, i = ii, length = ncol(peak.gene.associations)), 'sparseMatrix')
		return(v)
	}))

	associations = peak.gene.associations %*% gene.to.geneset.associations
	associations@x = rep(1, length(associations@x))


	scores = rowMaps(ace)[[sprintf("%s_feature_specificity", slot_name)]]
	if(max(scores) > 100) {
		scores = log1p(scores)
	}
	
	
	x = apply(scores, 1, max)
	y = Matrix::colSums(Matrix::t(associations))

	mask = (x > min.score) & (y > min.association)
	associations = associations[mask, ]
	scores = scores[mask, ]
		
	
	enrichment.out = assess_enrichment(scores, associations, thread_no)
	enrichment.mat = enrichment.out$logPvals
	rownames(enrichment.mat) = names(genesets)

	
	if(! ("enrichments" %in% names(metadata(ace))) ) {
		metadata(ace)[["enrichments"]] = list()
	}
	if( !(slot_name %in% names(metadata(ace)$enrichments)) ) {
		L = list(genesets = enrichment.mat)
		names(L) = geneset_slot_name
	} else {
		L = metadata(ace)$enrichments[[slot_name]]
		L[[geneset_slot_name]] = enrichment.mat
	}
	metadata(ace)$enrichments[[slot_name]] = L
	
	return(ace)
}


compute_TF_enrichment_from_ATACtion <- function(ace, min_peaks = 5, slot_name = "unified", TF_slot_name = "TFs", thread_no = 0, cisConnectome = "proximal_cisConnectome", min.score = 0, min.association = 0) {
	if(! ("motif_matches" %in% names(rowMaps(ace)))) {
		warning("motif_matches is not in rowMaps(ace). Please run add_motif_matched_to_ATACtion() first.")
		return(ace)
	}
	
	associations = as(rowMaps(ace)[["motif_matches"]], 'dgCMatrix')
	mask = (Matrix::colSums(associations) > min_peaks)
	associations = associations[, mask]


	scores = rowMaps(ace)[[sprintf("%s_feature_specificity", slot_name)]]
	if(max(scores) > 100) {
		scores = log1p(scores)
	}

	
	x = apply(scores, 1, max)
	y = Matrix::colSums(Matrix::t(associations))

	mask = (x > min.score) & (y > min.association)
	associations = associations[mask, ]
	scores = scores[mask, ]
	
	enrichment.out = assess_enrichment(scores, associations, thread_no)
	enrichment.mat = enrichment.out$logPvals
	rownames(enrichment.mat) = colnames(associations)

	if(! ("enrichments" %in% names(metadata(ace))) ) {
		metadata(ace)[["enrichments"]] = list()
	}
	if( !(slot_name %in% names(metadata(ace)$enrichments)) ) {
		L = list(TFs = enrichment.mat)
		names(L) = TF_slot_name
	} else {
		L = metadata(ace)$enrichments[[slot_name]]
		L[[TF_slot_name]] = enrichment.mat
	}
	metadata(ace)$enrichments[[slot_name]] = L
	
	return(ace)
}


compute_GRList_enrichment_from_ATACtion <- function(ace, GRlist, GR_slot_name = "GRs", slot_name = "unified", thread_no = 0, min.score = 0, min.association = 0) {
	rowGR = rowRanges(ace)
	DF = do.call(rbind, sapply(1:length(GRList), function(i) {
		GR = GRList[[i]]
		matches <- GenomicRanges::findOverlaps(rowGR, GR, select = "first", maxgap = -1, minoverlap = 1)
		ii = which(!is.na(matches))
		df = cbind(row = ii, col = rep(i, length(ii)))
		#v = as(sparseVector(x = 1, i = ii, length = length(rowGR)), 'sparseMatrix')
		
		return(df)
	}))
	associations = sparseMatrix(i = DF[, 1], j = DF[, 2], x = 1, dims = c(nrow(ace), length(GRList)))
	colnames(associations) = names(GRList)


	scores = rowMaps(ace)[[sprintf("%s_feature_specificity", slot_name)]]
	if(max(scores) > 100) {
		scores = log1p(scores)
	}

	
	x = apply(scores, 1, max)
	y = Matrix::colSums(Matrix::t(associations))

	mask = (x > min.score) & (y > min.association)
	associations = associations[mask, ]
	scores = scores[mask, ]
	
	
	enrichment.out = assess_enrichment(scores, associations, thread_no)
	enrichment.mat = enrichment.out$logPvals
	rownames(enrichment.mat) = colnames(associations)

	if(! ("enrichments" %in% names(metadata(ace))) ) {
		metadata(ace)[["enrichments"]] = list()
	}
	if( !(slot_name %in% names(metadata(ace)$enrichments)) ) {
		L = list(GRs = enrichment.mat)
		names(L) = GR_slot_name
	} else {
		L = metadata(ace)$enrichments[[slot_name]]
		L[[GR_slot_name]] = enrichment.mat		
	}
	metadata(ace)$enrichments[[slot_name]] = L
	
	return(ace)
}



annotate_ATACtion_cells_using_markers <- function(ace, markers, slot_name = "unified", postprocess = TRUE, LFR.threshold = 1.5) {
	ace = compute_geneset_enrichment_from_ATACtion(ace, markers, geneset_slot_name = "marker_enrichment", slot_name = slot_name, min.score = -1, min.association = -1)
	
	enrichment.mat = (Matrix::t(get_ATACtion_enrichment(ace, "marker_enrichment")))
	W = (enrichment.mat) #scale(enrichment.mat, center = T, scale = T)
	W[W < -log10(0.01/length(W))] = 0
	W = scale(W, center = F, scale = T)

    H.slot = sprintf("H_%s", slot_name)
    cell.scores.mat = colMaps(ace)[[H.slot]]
    cell.enrichment.mat = cell.scores.mat %*% W
    
    cell.annotations = colnames(cell.enrichment.mat)[apply(cell.enrichment.mat, 
        1, which.max)]
    Labels = colnames(cell.enrichment.mat)[apply(cell.enrichment.mat, 
        1, which.max)]
    Labels.confidence = apply(cell.enrichment.mat, 1, max)
    res = list(Labels = Labels, Labels.confidence = Labels.confidence, 
        Enrichment = cell.enrichment.mat)
    
    if(postprocess == T) {
    	res$Labels.updated = correct.cell.annotations(ace, res$Labels, LFR.threshold = LFR.threshold)	
    }
    return(res)			
}


project_ATACtion_enrichment_to_cells <- function(ace, enrichment.mat, slot_name = "unified", postprocess = TRUE) {
	W = Matrix::t(enrichment.mat) #scale(enrichment.mat, center = T, scale = T)
	W[W < -log10(0.01/length(W))] = 0
	W = scale(W, center = F, scale = T)
	W[is.na(W)] = 0
	

    H.slot = sprintf("H_%s", slot_name)
    cell.scores.mat = colMaps(ace)[["H_unified"]]
    cell.enrichment.mat = cell.scores.mat %*% W
    
    return(cell.enrichment.mat)			
}

