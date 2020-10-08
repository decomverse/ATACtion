compute_ATAC_archetype_to_ATAC_archetype_alignment <- function(reference_ace, query_ace, reference_slot_name = "unified", query_slot_name = "unified", deflate = FALSE, reduced_dim = 50, specificity_filter_threshold = 0) {
	
	reference_profile = rowMaps(reference_ace)[[sprintf("%s_feature_specificity", reference_slot_name)]]
	if(max(reference_profile) > 100) {
		reference_profile = log1p(reference_profile)
	}

	query_profile = rowMaps(query_ace)[[sprintf("%s_feature_specificity", query_slot_name)]]	
	if(max(query_profile) > 100) {
		query_profile = log1p(query_profile)
	}
	
	
	GR.reference = SummarizedExperiment::rowRanges(reference_ace)	
	GR.query = SummarizedExperiment::rowRanges(query_ace)
	

	mask.reference = apply(reference_profile, 1, max) > specificity_filter_threshold
	GR.reference = GR.reference[mask.reference]
	reference_profile = reference_profile[mask.reference, ]
	
	mask.query = apply(query_profile, 1, max) > specificity_filter_threshold	
	GR.query = GR.query[mask.query]
	query_profile = query_profile[mask.query, ]

	
	matches <- GenomicRanges::findOverlaps(GR.ace, GR.bulk, select = "first", maxgap = -1, minoverlap = 1)
	
	query_profile = query_profile[!is.na(matches), ]
	reference_profile = reference_profile[matches[!is.na(matches)], ]
	
	alignment = compute_pairwise_alignment(reference_profile = reference_profile, query_profile = query_profile, reduced_dim = reduced_dim, deflate = deflate)
	
	return(alignment)
}


compute_ATAC_cell_to_ATAC_cell_alignment <- function(reference_ace, query_ace, archetype_alignment = NULL, specificity_filter_threshold = 1, alignment_threshold = 0.1, footprint_threshold = 0.1, reference_slot_name = "unified", query_slot_name = "unified", deflate = FALSE, reduced_dim = 50) {	
	if(is.null(archetype_alignment)) {
		archetype_alignment = compute_ATAC_archetype_to_ATAC_archetype_alignment(reference_ace = reference_ace, query_ace = query_ace, reference_slot_name = reference_slot_name, query_slot_name = query_slot_name, deflate = deflate, reduced_dim = reduced_dim, specificity_filter_threshold = specificity_filter_threshold)
	}

	cell_to_cell_alignment = compute_cell_alignments_from_archetype_alignments(reference_ace = reference_ace, query_ace = query_ace, alignment = archetype_alignment, alignment_threshold = alignment_threshold, reference_slot_name = reference_slot_name, query_slot_name = query_slot_name, alignment_threshold = alignment_threshold, footprint_threshold = footprint_threshold)

	return(cell_to_cell_alignment)		
}


compute_bulkATAC_to_ATAC_archetype_alignment <- function(bulk_se, bulk_assay = "logcounts", ace,  slot_name = "unified", deflate = FALSE, reduced_dim = 50, specificity_filter_threshold = 1) {
		
	reference_profile = assays(bulk_se)[[bulk_assay]]
	query_profile = rowMaps(ace)[[sprintf("%s_feature_specificity", slot_name)]]
	if(max(query_profile) > 100) {
		query_profile = log1p(query_profile)
	}
	GR.ace = SummarizedExperiment::rowRanges(ace)	
	GR.bulk = SummarizedExperiment::rowRanges(bulk_se)
		
	x = apply(query_profile, 1, max)	
	mask = x > specificity_filter_threshold
	query_profile = query_profile[mask, ]
	GR.ace = GR.ace[mask]		
	
	matches <- GenomicRanges::findOverlaps(GR.ace, GR.bulk, select = "first", maxgap = -1, minoverlap = 1)
	
	query_profile = query_profile[!is.na(matches), ]
	reference_profile = reference_profile[matches[!is.na(matches)], ]
	
	alignment = compute_pairwise_alignment(reference_profile = reference_profile, query_profile = query_profile, reduced_dim = reduced_dim, deflate = deflate)
	
	return(alignment)
}
