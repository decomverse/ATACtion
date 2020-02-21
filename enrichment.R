identify.associated.REs <- function(peakset.ind.mat, pval.threshold = 0.001, frac.threshold = 0.1) {
  perm = order(a, decreasing = TRUE)
  a = a[perm]
  X = as(peakset.ind.mat[perm, ], 'dgCMatrix')
  
  p_c = Matrix::colMeans(X)
  
  p_r = Matrix::rowMeans(X)
  rho = mean(p_r)
  
  X.Null = p_r %*% t(p_c) / rho
  X.Null.sq = p_r %*% (t(p_c) / rho) ^ 2
  
  
  Obs = apply(X, 2, function(x) x*a)
  Exp = apply(X.Null, 2, function(x) x*a)
  Nu = apply(X.Null, 2, function(x) x*(a^2))
  
  a.max = a[1]
  
  associated.REs = unlist(sapply(1:ncol(X), function(j) {
    lambda = cumsum(Obs[, j]) - cumsum(Exp[, j])
    
    nu = Nu[, j]
    
    logPvals = (lambda^2) / (2*(nu + a.max*lambda/3))
    
    logPvals[lambda  < 0] = 0
    logPvals[is.na(logPvals)] = 0
    
    max.val = max(logPvals)
    max.idx = which.max(logPvals)
  
    if(max.val > log(pval.threshold) + log(length(X))) {
      idx = which(X[1:max.idx, j] == 1)
      return(perm[idx])
    } else {
      return ()
    }
  }))

  #counts = table(associated.REs)
  #associated.REs = as.numeric(names(counts)[counts > round(ncol(X)*frac.threshold)])
  associated.REs = sort(unique(associated.REs))
  return(associated.REs)
}

peakset.Enrichment <- function(A, peakset.ind.mat) {
	if(max(A) > 20) {	
		A = log(1+A)	
	}
	
	if(nrow(A) != nrow(peakset.ind.mat)) {
		message("peakset.Enrichment:: Number of rows do not match");
		return()		
	}
	
	X = as(peakset.ind.mat, 'dgCMatrix')

	p_c = Matrix::colMeans(X)

	p_r = Matrix::rowMeans(X)
	rho = mean(p_r)
	
	Obs = as.matrix(Matrix::t(X) %*% A)
	Exp = as.matrix( p_c %*% ( Matrix::t(p_r) %*% A) ) / rho 
	Nu = as.matrix( (p_c^2) %*% (Matrix::t(p_r) %*% (A^2)) ) / (rho^2) 

	Lambda = Obs - Exp
	
	a = apply(A, 2, max)  
	logPvals = Lambda^2 / (2*(Nu + (p_c %*% t(a))*Lambda/3))
	
	rownames(logPvals) = colnames(X)
	colnames(logPvals) = colnames(A)

	logPvals[Lambda  < 0] = 0
	logPvals[is.na(logPvals)] = 0
	
	logPvals = logPvals / log(10)
			
	return(logPvals)
}

peakset.Enrichment.noRowScaling <- function(A, peakset.ind.mat) {
	if(max(A) > 20) {	
		A = log(1+A)	
	}
	
	if(nrow(A) != nrow(peakset.ind.mat)) {
		message("peakset.Enrichment:: Number of rows do not match");
		return()		
	}
	
	X = as(peakset.ind.mat, 'dgCMatrix')

	p_c = Matrix::colMeans(X)

	p_r = Matrix::rowMeans(X)
	rho = mean(p_r)
	
	X.Null = p_r %*% t(p_c) / rho
	X.Null.sq = p_r %*% (t(p_c) / rho) ^ 2
	
	Obs = as.matrix(Matrix::t(X) %*% A)
	Exp = as.matrix(Matrix::t(X.Null) %*% A)	
	Nu = as.matrix(Matrix::t(X.Null.sq) %*% (A^2))

	Lambda = Obs - Exp
	
	a = apply(A, 2, max)  
	logPvals = Lambda^2 / (2*(Nu + (p_c %*% t(a))*Lambda/3))
	
	
	rownames(logPvals) = colnames(X)
	colnames(logPvals) = colnames(A)

	logPvals[Lambda  < 0] = 0
	logPvals[is.na(logPvals)] = 0
	
	logPvals = logPvals / log(10)
			
	return(logPvals)
}



atACTION.RE.enrichment <- function(A, R, GRList) {
	require(GenomicRanges)

	if(is.GR(R)) {
		GR = R
	}
	else {
		GR = rowRanges(R)
	}
	
	if(is.matrix(GRList) | is.sparseMatrix(GRList)) {
		Ind.mat = as(GRList, 'sparseMatrix')
	} else {
		values(GR) = c()
		if(is.matrix(A) | is.sparseMatrix(A)) {
			A = A
		} else if(is.sce(A)) {
			A = ACTIONet.out@assays[[A.slot]]
		} else {
			A = ACTIONet.out$archetype.accessibility@assays[[A.slot]]
		}
		if(max(A) > 10 & min(A) > 0)
			A = log(1 + A)
			
		values(GR) = A
		
		src.genome = genome(GR)[[1]]
		dst.genome = genome(GRList[[1]])[[1]]
			
		if(src.genome != dst.genome) {
			print("liftOver GR")
			GR = liftOverGR(GR, src.genome, dst.genome)
		}
		
		print("Compute overlap with query sets")
		Ind.mat = sapply(GRList, function(query) {
			#return(as.numeric(GR %in% query))
			
			matches <- GenomicRanges::findOverlaps(query, GR, select = "all", maxgap = -1, minoverlap = 1)
			ind = as.numeric(Matrix::sparseVector(x = 1, i = sort(unique(S4Vectors::subjectHits(matches))), length = length(GR)))
		})
		colnames(Ind.mat) = paste(names(GRList))
	}
	
	scores = as.matrix(values(GR))
	
	#excluded.REs = identify.associated.REs(Ind.mat, Matrix::rowMeans(scores))	
	#selected.REs = setdiff(1:nrow(Ind.mat), excluded.REs)
	
	Enrichment.profile = peakset.Enrichment.noRowScaling(scores, Ind.mat)
	rownames(Enrichment.profile) = names(GRList)
	
	out.list = list(Enrichment.profile = Enrichment.profile, Ind.mat = Ind.mat)
	
	return(out.list)
}


atACTIONet.GWAS.enrichment <- function(ACTIONet.out, sce, varGR, split.var = "DISEASE/TRAIT", min.count = 30, flank = 100) {
  
  GR = rowRanges(sce)
  values(GR) = c()
  
  src.genome = genome(GR)[[1]]
  if(src.genome == 'GRCh38') {
    src.genome = 'hg38'
    genome(GR) = 'hg38'
  }
  dst.genome = genome(varGR)[[1]]
  if(dst.genome == 'GRCh38') {
    dst.genome = 'hg38'
    genome(varGR) = 'hg38'
  }
  
  seqlevelsStyle(varGR) <- "UCSC"
  seqlevelsStyle(GR) <- "UCSC"
  
  counts = table(values(varGR)[[split.var]])
  mask = values(varGR)[[split.var]] %in% names(counts)[counts >= min.count]
  
  varGR = varGR[mask]
  GRList = split(varGR, values(varGR)[[split.var]])


  A = ACTIONet.out$archetype.accessibility@assays[["signature"]]
  values(GR) = A
	if(src.genome != dst.genome) {
		print("liftOver GR")
		GR = liftOverGR(GR, src.genome, dst.genome)
	}

	print("Compute overlap with query sets")
	Ind.mat = sapply(GRList, function(query) {
	  start(ranges(query)) = start(ranges(query)) - flank
	  end(ranges(query)) = end(ranges(query)) + flank
	  v = as.numeric(GR %over% query)
		return(v)
	})
	
	colnames(Ind.mat) = paste(names(GRList))

	



	varEnrichment.profile = peakset.Enrichment(as.matrix(values(GR)), Ind.mat)

	rownames(varEnrichment.profile) = names(GRList)

	return(varEnrichment.profile)
}



atACTION.motif.enrichment <- function(ACTIONet.out, sce, log.transform = TRUE, data.slot = "signature") {

	if(log.transform == TRUE)
		A = log(1+ACTIONet.out$archetype.accessibility@assays[[data.slot]])	
	else
		A = ACTIONet.out$archetype.accessibility@assays[[data.slot]]
	
	if( !("matched.motifs" %in% names(rowData(sce))) ) {
		print("Motif matching matrix doesn't exist in the ACTIONet.out. Running motif matching first");
		sce = atACTIONet.match.motifs(sce)
	}  	
	
	X = as(rowData(sce)[["matched.motifs"]], 'dgCMatrix')

	p_c = Matrix::colMeans(X)

	p_r = Matrix::rowMeans(X)
	rho = mean(p_r)
	
	X.Null = p_r %*% t(p_c) / rho
	X.Null.sq = p_r %*% (t(p_c) / rho) ^ 2
	
	Obs = as.matrix(Matrix::t(X) %*% A)
	Exp = as.matrix(Matrix::t(X.Null) %*% A)	
	Nu = as.matrix(Matrix::t(X.Null.sq) %*% (A^2))

	Lambda = Obs - Exp
	
	a = apply(A, 2, max)  
	logPvals = Lambda^2 / (2*(Nu + (p_c %*% t(a))*Lambda/3))
	
	rownames(logPvals) = colnames(X)
	colnames(logPvals) = colnames(A)

	logPvals[Lambda  < 0] = 0
	logPvals[is.na(logPvals)] = 0
	
	logPvals = logPvals / log(10)
			
	return(logPvals)
}


atACTIONet.proximalGene.enrichment <- function(ACTIONet.out, sce, genes = NA) {
	
	gene.Ind = metadata(sce)[["cisConnectome"]]
	sig.profile = ACTIONet.out$archetype.accessibility@assays[["signature"]]

	if(!is.na(genes)) {
		gene.Ind = gene.Ind[, intersect(colnames(gene.Ind), genes)]
	}


	Enrichment.profile = peakset.Enrichment(sig.profile, gene.Ind)
	rownames(Enrichment.profile) = colnames(gene.Ind)
		
	return(Enrichment.profile)	
	
}

atACTIONet.geneSet.enrichment <- function(ACTIONet.out, sce, annotations) {
	
	proximal.map = metadata(sce)[["cisConnectome"]]
	#if(is.sce(ACTIONet.out)) 

	if(log.transform == TRUE)
		A = log(1+ACTIONet.out$archetype.accessibility@assays[[data.slot]])	
	else
		A = ACTIONet.out$archetype.accessibility@assays[[data.slot]]
		
	
	if(is.list(annotations)) {
		total.genes = sort(unique(unlist(annotations)))

		annotations = sapply(annotations, function(gs) {
			return(as.numeric(total.genes %in% gs))
		})
		rownames(annotations) = total.genes
		annotations = as(annotations, 'sparseMatrix')		
	}

	common.genes = intersect(colnames(proximal.map), rownames(annotations))

	X = proximal.map[, common.genes]
	Y = annotations[common.genes, ]

	ind.mat = X %*% Y
	ind.mat@x[ind.mat@x > 1] = 1
	
	counts = Matrix::rowSums(ind.mat)
	mask = counts > 0

	A = A[mask, ]
	Z = ind.mat[mask, ]

	p_c = Matrix::colMeans(Z)


	Obs = as.matrix(Matrix::t(Z) %*% A)
	Exp = as.matrix(p_c %*% Matrix::t(Matrix::colSums(A)))
	Nu = as.matrix(p_c %*% Matrix::t(Matrix::colSums(A^2)))

	Lambda = Obs - Exp

	a = apply(A, 2, max)
	ones = array(1, dim = dim(Nu)[1])

	logPvals = Lambda^2 / (2*(Nu + (ones %*% t(a))*Lambda/3))

	rownames(logPvals) = colnames(Z)
	colnames(logPvals) = colnames(A)

	logPvals[Lambda  < 0] = 0
	logPvals[is.na(logPvals)] = 0

	logPvals = logPvals / log(10)  

	return(logPvals)	
}

is.sparseMatrix <- function(aa) {
	return(length(which(is(aa)=="sparseMatrix"))!=0)
}

is.sce <- function(sce) {
	return(length(which(is(sce)=="SingleCellExperiment"))!=0)
}

is.GR <- function(GR) {
	return(length(which(is(GR)=="GenomicRanges"))!=0)
}

# From: Testing for Network and Spatial Autocorrelation (CRAN:: netdep
network.autocorrlation.Labels <- function(A, Labels) {
	if(is.igraph(A)) {
		A = as(get.adjacency(A, attr = "weight"), 'dgTMatrix')		
	} else if(is.list(A)) {
		A = as(A$build.out$ACTIONet, 'dgTMatrix')
	} else if(is.matrix(A) | is.sparseMatrix(A)) {
		A = as(A, 'dgTMatrix')				
	} else {
		print("Unknown format for A");
		print(class(A))
		return()
	}
	
	if(!is.factor(Labels)) {
		Labels = factor(Labels, levels = sort(unique(Labels)))
	}

	n = length(Labels)
	counts = table(Labels)
	pvec = as.numeric(counts / sum(counts))
	names(pvec) = levels(Labels)	
	k = length(pvec)

	# Compute statistics
	L_i = Labels[A@i+1];
	L_j = Labels[A@j+1];

	p_i = pvec[L_i]
	p_j = pvec[L_j]

	rawphi = sum(A@x*(2*(L_i == L_j)-1) / (p_i*p_j))

	# Compute moments
	s0 = sum(A + Matrix::t(A)) / 2
	s1 = sum((A + Matrix::t(A))^2) / 2
	s2 = sum( (Matrix::colSums(A) + Matrix::rowSums(A))^2 )

	m1.rawphi =  (s0/(n*(n-1)))*(n^2*k*(2-k) - n*sum(1/pvec))

	Q1 = sum(1/pvec)
	Q2 = sum(1/pvec^2)
	Q3 = sum(1/pvec^3)
	Q22 = sum((1/pvec)%*%t(1/pvec))
	E1 = (n^2*Q22 - n*Q3)/(n*(n-1))
	E2 = 4*n^3*Q1 - 4*n^3*k*Q1 + n^3*k^2*Q1 - 2*( 2*n^2*Q2 - n^2*k*Q2) + 2*n*Q3 - n^2*Q22
	E2 = E2/(n*(n-1)*(n-2))

	A1 = 4*n^4*k^2 - 4*n^4*k^3 + n^4*k^4 - (2*n^3*k*Q1 - n^3*k^2*Q1)
	A2 = 4*n^3*Q1 - 4*n^3*k*Q1 + n^3*k^2*Q1 - (2*n^2*Q2 - n^2*k*Q2)
	Apart = A1 - 2*A2

	B1 = 4*n^3*Q1 - 4*n^3*k*Q1 + n^3*k^2*Q1 - (2*n^2*Q2 - n^2*k*Q2)
	B2 = 2*n^2*Q2 - n^2*k*Q2 - n*Q3
	B3 = n^2*Q22 - n*Q3
	Bpart = B1 - B2 - B3

	C1 = 2*n^3*k*Q1 - n^3*k^2*Q1 - n^2*Q22
	C2 = 2*n^2*Q2 - n^2*k*Q2 - n*Q3
	Cpart = C1 - 2*C2

	E3 = (Apart - 2*Bpart - Cpart) / (n*(n-1)*(n-2)*(n-3))

	m2.rawphi = s1*E1 + (s2 - 2*s1)*E2 + (s0^2 - s2 + s1)*E3
	
		
	# Compute summary statistics
	mean.rawphi = m1.rawphi
	var.rawphi = m2.rawphi - mean.rawphi^2
	
	phi.z = (rawphi - mean.rawphi) / sqrt(var.rawphi)
	phi.logPval = -log10(pnorm(phi.z, lower.tail = FALSE))

	return(list(z = phi.z, phi = rawphi, mean.rawphi = mean.rawphi, var.rawphi = var.rawphi, logPval = phi.logPval))
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
