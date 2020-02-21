plotEnrichmentHeatmap <- function(Enrichment.profile, rowAnnot, top.terms = 5, rowCPal = "d3", colCPal = "d3", row_title = NULL, col_title = "Cell states", row.normalize = TRUE) {

  #selected.rows= sort(unique(as.numeric(apply(Enrichment.profile, 2, function(x) {
  #  return(order(x, decreasing = T)[1:top.terms])
  #}))))
  

  max.scores = apply(Enrichment.profile, 1, max)
  nnz = round(sum(max.scores)^2 / sum(max.scores^2))

  selected.rows = order(max.scores, decreasing = T)[1:top.terms]


  
  
  require(ComplexHeatmap)
  
  if( !is.null(colnames(Enrichment.profile)) ) {
    L = sort(unique(colnames(Enrichment.profile)))
    
    if(length(colCPal) > 1) {
      if(is.null(names(colCPal))) {
		colPal = colCPal[1:length(L)]
        names(colCPal) = L        
      }
      else {
		colPal = colCPal[L]		 
	  }
    } else {
      colPal = ggpubr::get_palette(colCPal, length(unique(L)))
      names(colPal) = L
    }
    
    
    ha_cols = HeatmapAnnotation(df = list(Celltype = as.character(colnames(Enrichment.profile))), col = list(Celltype = colPal), annotation_legend_param = list(Celltype=list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 5))), which = "column")
  } else {
    ha_cols = NULL
  }
  
  if( !is.null(rowAnnot) ) {
    L = sort(unique(rowAnnot[selected.rows]))
    if(is.list(rowCPal)) {        
      if(is.null(names(rowPal))) {
		rowPal = rowPal[1:length(L)]
        names(colCPal) = L        
      }
      else {
		rowPal = rowPal[L]		 
	  }                        
    } else {
      rowPal = ggpubr::get_palette(rowCPal, length(L))
      names(rowPal) = L
    }
  
    
    ha_rows = HeatmapAnnotation(df = list(Group = rowAnnot[selected.rows]), col = list(Group = rowPal), annotation_legend_param = list(Group=list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 5))), which = "row", annotation_name_side = "top")
  } else {
    ha_rows = NULL
  }
  
  
  
  gradPal = grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = 9, name = "RdYlBu")))(100)

	if(row.normalize == TRUE) {
		Z = t(scale(t(Enrichment.profile)))
		Z[is.na(Z)] = 0
		Z[Z > 3] = 3
	} else {
		Z = Enrichment.profile
	}
  
  Heatmap(Z[selected.rows, ], col = gradPal, row_names_gp = gpar(fontsize = 8), column_names_gp = gpar(fontsize = 0), left_annotation = ha_rows, top_annotation = ha_cols, name = "Enrichment", column_title = col_title, row_title = row_title)

}  





plot.ACTIONet.motif.view <- function(ACTIONet.out, Enrichment, top.motifs = 3, arch.colors = NA) {

  if( nrow(Enrichment) == nrow(ACTIONet.out$reconstruct.out$H_stacked) )
    Enrichment = t(Enrichment)
  
  
  
  
  	Enrichment.core = exp(2*scale(t(scale(t(Enrichment[, ACTIONet.out$core.out$core.archs])))))
  	Enrichment.core[is.na(Enrichment.core)] = 0
  	rs = Matrix::rowSums(Enrichment.core)
  	rs[rs == 0] = 1
  	Enrichment.norm = t(scale(t(Enrichment.core), center = FALSE, scale = rs))
  
  
  # selected.rows = sort(unique(unlist(apply(Enrichment, 2, function(x) {
  # 		nnz = round(sum(abs(x^2))^2 / sum(x^4))
  # 		
  # 		rows = order(x, decreasing = TRUE)[1:nnz]
  # }))))
  
  selected.rows= sort(unique(as.numeric(apply(Enrichment, 2, function(x) {
    return(order(x, decreasing = T)[1:top.motifs])
  }))))
  
  Enrichment.norm = Enrichment.norm[selected.rows, ]
  
  colnames(Enrichment.norm) = archLabels.org[ACTIONet.out$core.out$core.archs]
  
  arch.coors = ACTIONet.out$arch.vis.out$coordinates
  motif.coors = Enrichment.norm %*% arch.coors[ACTIONet.out$core.out$core.archs, ]
    
  if(!is.na(arch.colors)) {
    motif.colors = as.character(arch.colors[ACTIONet.out$core.out$core.archs][apply(Enrichment.norm, 1, which.max)])
  } else {
    arch.colors = ACTIONet.out$arch.vis.out$colors
    arch.colors.RGB = col2rgb(arch.colors[ACTIONet.out$core.out$core.archs])/256
    
    
    arch.colors.Lab = grDevices::convertColor(color = t(arch.colors.RGB), from = "sRGB",
        to = "Lab")
    
    motif.colors.Lab = Enrichment.norm %*% arch.colors.Lab
    motif.colors = rgb(grDevices::convertColor(color = motif.colors.Lab, from = "Lab", to = "sRGB"))
  }
  motif.colors = colorspace::darken(motif.colors, 0.1)
  
  motifs.df = data.frame(motif = rownames(Enrichment.norm), x = motif.coors[, 1], y = motif.coors[, 2])
  
  names(motif.colors) = rownames(Enrichment.norm)
  
  
  require(ggrepel)
  require(ggplot2)
  p <- ggplot(motifs.df, aes(x, y, label = motif, color = motif)) + 
      scale_colour_manual(values = motif.colors) + geom_point(show.legend = FALSE) + 
      geom_label_repel(show.legend = FALSE, force = 5) + theme_void()
  plot(p)

}
