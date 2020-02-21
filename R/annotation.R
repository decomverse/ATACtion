annotate.ATAC.clusters.using.markers <- function(ACTIONet.out, annotation.cluster, sce, marker.genes, flank.size = 10000, reduction.slot = "ACTION", sce.data.attr = "bin_counts") {
	cl.idx = which(names(ACTIONet.out$annotations) == annotation.cluster)
	if(length(cl.idx) == 0) {
		R.utils::printf('Error in correct.cell.labels: annotation.cluster "%s" not found\n', annotation.cluster)
		return(ACTIONet.out)
	}		
	clusters = ACTIONet.out$annotations[[cl.idx]]$Labels    
	clusters = preprocess.labels(ACTIONet.out, clusters)


	if( is.null(ACTIONet.out$annotations[[cl.idx]]$DE.profile) ) {
		ACTIONet.out = compute.annotations.feature.specificity(ACTIONet.out, sce, annotation.cluster, sce.data.attr = sce.data.attr)
	}
	A = log1p(as.matrix(ACTIONet.out$annotations[[cl.idx]]$DE.profile@assays[["significance"]]))


	if(! ("cisConnectome" %in% names(metadata(sce))) ) {
	  sce = initialize.cisConnectome(sce, flank = flank.size)
	}
	cisConnectome = metadata(sce)[["cisConnectome"]]
	
	markerInd = sapply(marker.genes, function(genes) {
	  genes = intersect(genes, colnames(cisConnectome))
	  v = as.numeric(Matrix::rowSums(cisConnectome[, genes]) > 0)
	})


	celltype.Enrichment = Chernoff.enrichment(A, markerInd)

    Annot = names(marker.genes)    
    clusterLabels = Annot[apply(celltype.Enrichment, 2, which.max)]
   
    cellLabels = match(clusterLabels[clusters], Annot)
    names(cellLabels) = clusterLabels[clusters]
    
    fullLabels = sapply(sort(unique(clusters)), function(i) {
        return(sprintf("Cluster %d (%s)", i, clusterLabels[[i]]))
    })

    cellFullLabels = match(fullLabels[clusters], fullLabels)
    names(cellFullLabels) = fullLabels[clusters]
    
    res = list(Labels = clusterLabels, cellLabels = cellLabels, fullLabels = fullLabels, cellFullLabels = cellFullLabels, Enrichment = celltype.Enrichment)
    
    return(res)
    
    #ACTIONet.out$annotations[[cl.idx]]$markerEnrichment = res    
    #return(ACTIONet.out)
	# clusterCelltypes = factor(rownames(celltype.Enrichment)[apply(celltype.Enrichment, 2, which.max)])
	# clusterLabels = sapply(1:length(clusterCelltypes), function(i) sprintf('%s (%s)', i, clusterCelltypes[[i]]))
	# clusterLabels = factor(clusterLabels, levels = clusterLabels)
	# 
	# clusterInferredLabels = clusterLabels[as.numeric(ACTIONet.out$annotations[[annot.idx]]$Labels)]
	# clusterInferredCelltypes = clusterCelltypes[as.numeric(ACTIONet.out$annotations[[annot.idx]]$Labels)]
	# 
	# out.list = list(clusterLabels = clusterLabels, clusterCelltypes = clusterCelltypes, clusterInferredLabels = clusterInferredLabels, clusterInferredCelltypes = clusterInferredCelltypes, Enrichment = celltype.Enrichment)
	# 
	# return(out.list)	
}


annotate.ATAC.cells.using.markers <- function(ACTIONet.out, sce, marker.genes, annotation.name, confidence.threshold = 10, post.update = T, flank.size = 10000) {
	arch.annot.out = annotate.ATAC.archetypes.using.markers(ACTIONet.out, sce, marker.genes, flank.size = flank.size, core = T)
	confidence = apply(arch.annot.out$Enrichment, 1, max)
	arch.labels = as.character(arch.annot.out$Labels)
	arch.labels[confidence < confidence.threshold] = "?"

	cell.annotations = arch.labels[ACTIONet.out$unification.out$assignments.core]

	ACTIONet.out = add.cell.annotations(ACTIONet.out, cell.annotations = cell.annotations, annotation.name = annotation.name)

    if(post.update == T) {
		ACTIONet.out = correct.cell.annotations(ACTIONet.out, annotation.in = annotation.name, annotation.out = annotation.name, adjust.levels = TRUE, highlight = F)
	}
	
	return(ACTIONet.out)
}


annotate.ATAC.archetypes.using.markers <- function(ACTIONet.out, sce, marker.genes, flank.size = 10000, core = T) {
	if(! ("cisConnectome" %in% names(metadata(sce))) ) {
	  sce = initialize.cisConnectome(sce, flank = flank.size)
	}
	cisConnectome = metadata(sce)[["cisConnectome"]]
	
	markerInd = sapply(marker.genes, function(genes) {
	  genes = intersect(genes, colnames(cisConnectome))
		v = as.numeric(Matrix::rowSums(cisConnectome[, genes]) > 0)
	})

	if(core == T) {
		if (("unification.out" %in% names(ACTIONet.out))) {
			print("Using unification.out$DE.core (merged archetypes)")
			A = as.matrix(log1p(ACTIONet.out$unification.out$DE.core@assays[["significance"]]))
			
		} else {
			print("unification.out is not in ACTIONet.out. Please run unify.cell.states() first.")
			return()
		}
	} else {
		if (("archetype.differential.signature" %in% names(ACTIONet.out))) {
			print("Using archetype.differential.signature (all archetypes)")
			A = as.matrix(log1p(ACTIONet.out$archetype.differential.signature@assays[["significance"]]))
		} else {
			print("archetype.differential.signature is not in ACTIONet.out. Please run compute.archetype.feature.specificity() first.")
			return()
		}
	}      
	
	celltype.Enrichment = Chernoff.enrichment(A, markerInd)
	
	rownames(celltype.Enrichment) = names(marker.genes)
	
	idx = apply(celltype.Enrichment, 2, which.max)	
	Labels = rownames(celltype.Enrichment)[idx]	
	Labels = factor(Labels, levels = names(marker.genes))
	
	out.list = list(Labels = Labels, Enrichment = t(celltype.Enrichment))
	
	return(out.list)	
}

