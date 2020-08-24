compute_RNA_archetype_to_ATAC_archetype_alignment <- function(RNA_ace, ATAC_ace, RNA_prefix = "unified", ATAC_enrichment_slot = "unified", ATAC_enrichment_name = "genes", deflate = FALSE, reduced_dim = 50, specificity_filter_threshold = 1) {
		
	reference_profile = rowMaps(RNA_ace)[[sprintf("%s_feature_specificity", RNA_prefix)]]
	query_profile = get_ATACtion_enrichment(ATAC_ace, enrichment_name = ATAC_enrichment_name, slot_name = ATAC_enrichment_slot)
	
	g1 = rownames(reference_profile)[apply(reference_profile, 1, max) > specificity_filter_threshold]
	g2 = rownames(query_profile)[apply(query_profile, 1, max) > specificity_filter_threshold]
	common.genes = intersect(g1, g2)
		
	reference_profile = reference_profile[common.genes, ]
	query_profile = query_profile[common.genes, ]

	alignment = compute_pairwise_alignment(reference_profile = reference_profile, query_profile = query_profile, reduced_dim = reduced_dim, deflate = deflate)
	
	return(alignment)
}

compute_RNA_cell_to_ATAC_cell_alignment <- function(RNA_ace, ATAC_ace, archetype_alignment = NULL, specificity_filter_threshold = 1, alignment_threshold = 0.1, footprint_threshold = 0.1, RNA_prefix = "unified", ATAC_enrichment_slot = "unified", ATAC_enrichment_name = "genes", deflate = FALSE, reduced_dim = 50) {
	if(is.null(archetype_alignment)) {
		archetype_alignment = compute_RNA_archetype_to_ATAC_archetype_alignment(RNA_ace = RNA_ace, ATAC_ace = ATAC_ace, RNA_prefix = RNA_prefix, ATAC_enrichment_slot = ATAC_enrichment_slot, ATAC_enrichment_name = ATAC_enrichment_name, deflate = deflate, reduced_dim = reduced_dim, specificity_filter_threshold = specificity_filter_threshold)
	}
	
	cell_to_cell_alignment = compute_cell_alignments_from_archetype_alignments(reference_ace = RNA_ace, query_ace = ATAC_ace, alignment = archetype_alignment, alignment_threshold = alignment_threshold, reference_slot_name = reference_slot_name, query_slot_name = query_slot_name, alignment_threshold = alignment_threshold, footprint_threshold = footprint_threshold)

	return(cell_to_cell_alignment)	
}

compute_ATAC_archetype_to_RNA_archetype_alignment <- function(RNA_ace, ATAC_ace, RNA_prefix = "unified", ATAC_enrichment_slot = "unified", ATAC_enrichment_name = "genes", deflate = FALSE, reduced_dim = 50, specificity_filter_threshold = 1) {
		
	reference_profile = get_ATACtion_enrichment(ATAC_ace, enrichment_name = ATAC_enrichment_name, slot_name = ATAC_enrichment_slot)
	query_profile = rowMaps(RNA_ace)[[sprintf("%s_feature_specificity", RNA_prefix)]]
	
	g1 = rownames(reference_profile)[apply(reference_profile, 1, max) > specificity_filter_threshold]
	g2 = rownames(query_profile)[apply(query_profile, 1, max) > specificity_filter_threshold]
	common.genes = intersect(g1, g2)
		
	reference_profile = reference_profile[common.genes, ]
	query_profile = query_profile[common.genes, ]

	alignment = compute_pairwise_alignment(reference_profile = reference_profile, query_profile = query_profile, reduced_dim = reduced_dim, deflate = deflate)
	
	return(alignment)
}

compute_ATAC_cell_to_RNA_cell_alignment <- function(RNA_ace, ATAC_ace, archetype_alignment = NULL, specificity_filter_threshold = 1, alignment_threshold = 0.1, footprint_threshold = 0.1, RNA_prefix = "unified", ATAC_enrichment_slot = "unified", ATAC_enrichment_name = "genes", deflate = FALSE, reduced_dim = 50) {
	if(is.null(archetype_alignment)) {
		archetype_alignment = compute_ATAC_archetype_to_RNA_archetype_alignment(RNA_ace = RNA_ace, ATAC_ace = ATAC_ace, RNA_prefix = RNA_prefix, ATAC_enrichment_slot = ATAC_enrichment_slot, ATAC_enrichment_name = ATAC_enrichment_name, deflate = deflate, reduced_dim = reduced_dim, specificity_filter_threshold = specificity_filter_threshold)
	}
	
	cell_to_cell_alignment = compute_cell_alignments_from_archetype_alignments(reference_ace = ATAC_ace, query_ace = RNA_ace, alignment = archetype_alignment, alignment_threshold = alignment_threshold, reference_slot_name = reference_slot_name, query_slot_name = query_slot_name, alignment_threshold = alignment_threshold, footprint_threshold = footprint_threshold)

	return(cell_to_cell_alignment)	
}



compute_bulkRNA_to_ATAC_archetype_alignment <- function(bulk_RNA_profile, ATAC_ace,  ATAC_enrichment_slot = "unified", ATAC_enrichment_name = "genes", deflate = FALSE, reduced_dim = 50, specificity_filter_threshold = 1) {
		
	reference_profile = bulk_RNA_profile#orthoProject(bulk_RNA_profile, Matrix::rowMeans(bulk_RNA_profile))
	query_profile = get_ATACtion_enrichment(ATAC_ace, enrichment_name = ATAC_enrichment_name, slot_name = ATAC_enrichment_slot)
		
	gg = rownames(query_profile)[apply(query_profile, 1, max) > specificity_filter_threshold]
	common.genes = intersect(rownames(reference_profile), gg)
	
	reference_profile = reference_profile[common.genes, ]
	query_profile = query_profile[common.genes, ]

	alignment = compute_pairwise_alignment(reference_profile = reference_profile, query_profile = query_profile, reduced_dim = reduced_dim, deflate = deflate)
	
	return(alignment)
}
