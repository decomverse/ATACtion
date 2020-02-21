# A is always the DE matrix
ATACtion.GRList.enrichment <- function(Diff.sce, GRList) {
	require(GenomicRanges)
	A = log1p(Diff.sce@assays[["significance"]])
	GR = rowRanges(Diff.sce)

	if(is.matrix(GRList) | is.sparseMatrix(GRList)) {
		Ind.mat = as(GRList, 'sparseMatrix')
		scores = as.matrix(A)		
	} else {
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
		
		scores = as.matrix(values(GR))
	}
	
	
	Enrichment.profile = Chernoff.enrichment(scores, Ind.mat)
	rownames(Enrichment.profile) = names(GRList)
	
	out.list = list(Enrichment.profile = Enrichment.profile, Ind.mat = Ind.mat)
	
	return(out.list)
}

ATACtion.archetype.enhancer.enrichment <- function(ACTIONet.out, enrichment.name = NULL, core = T, sce, reduction_slot = "ACTION", sce.data.attr = "bin_counts", enhancer_DB = "Roadmap") {
	if(core == T) {
		if( !("unification.out" %in% names(ACTIONet.out)) ) {
			print("Running unification")
			ACTIONet.out$unification.out = unify.cell.states(ACTIONet.out, sce, reduction_slot = "ACTION", sce.data.attr = sce.data.attr)    
		} 
		DE.profile = ACTIONet.out$unification.out$DE.core		
	} else {
		if( !("unification.out" %in% names(ACTIONet.out)) ) {
			print("Running archetype DE")
			ACTIONet.out$archetype.differential.signature = compute.archetype.feature.specificity(ACTIONet.out, sce, mode = "dense", sce.data.attr = sce.data.attr)
		}
		DE.profile = ACTIONet.out$archetype.differential.signature
	}
	values(DE.profile) = c()

	species = genome(DE.profile)[[1]][[1]]
	if(! (species %in% c("hg19", "hg38")) ) {
		R.utils::printf("Genome %s is not recognized\n", species)
		return(ACTIONet.out)
	}
	db.name = sprintf('EnhDB_%s_%s', enhancer_DB, species)


	if ( !exists(db.name) ) {
		cmd = sprintf('data(\"%s\")', db.name)
		eval(parse(text=cmd))
	}
	
	cmd = sprintf('Enrichment.out = ATACtion.GRList.enrichment(DE.profile, %s)', db.name)
	eval(parse(text=cmd))
	
	if( is.null(enrichment.name) ) {
		enrichment.name = db.name
	}
	if(core == T) {
		if( !("Enrichments" %in% names(ACTIONet.out$unification.out)) ) {
			ACTIONet.out$unification.out$Enrichments = list()
		}
		cmd = sprintf('ACTIONet.out$unification.out$Enrichments$\"%s\" = Enrichment.out', enrichment.name)
		eval(parse(text=cmd))
	} else {
		if( !("Enrichments" %in% names(ACTIONet.out)) ) {
			ACTIONet.out$Enrichments = list()
		}
		cmd = sprintf('ACTIONet.out$Enrichments$\"%s\" = Enrichment.out', enrichment.name)
		eval(parse(text=cmd))		
	}

	
	return(ACTIONet.out)
	
			
}

ATACtion.cluster.enhancer.enrichment <- function(ACTIONet.out, annotation.cluster, enrichment.name = NULL, sce, sce.data.attr = "bin_counts", enhancer_DB = "Roadmap") {
	cl.idx = which(names(ACTIONet.out$annotations) == annotation.cluster)
	if(length(cl.idx) == 0) {
		R.utils::printf('Error in ATACtion.cluster.enhancer.enrichment: annotation.cluster "%s" not found\n', annotation.cluster)
		return(ACTIONet.out)
	}		
	if(is.null(ACTIONet.out$annotations[[cl.idx]]$DE.profile)) {
		ACTIONet.out = compute.annotations.feature.specificity(ACTIONet.out, sce, annotation.name = annotation.cluster, sce.data.attr = sce.data.attr)
	}
	DE.profile = ACTIONet.out$annotations[[cl.idx]]$DE.profile
	values(DE.profile) = c()

	species = genome(DE.profile)[[1]][[1]]
	if(! (species %in% c("hg19", "hg38")) ) {
		R.utils::printf("Genome %s is not recognized\n", species)
		return(ACTIONet.out)
	}
	db.name = sprintf('EnhDB_%s_%s', enhancer_DB, species)


	if ( !exists(db.name) ) {
		cmd = sprintf('data(\"%s\")', db.name)
		eval(parse(text=cmd))
	}
	
	cmd = sprintf('Enrichment.out = ATACtion.GRList.enrichment(DE.profile, %s)', db.name)
	eval(parse(text=cmd))
	
	if( is.null(enrichment.name) ) {
		enrichment.name = db.name
	}
	if( !("Enrichments" %in% names(ACTIONet.out$annotations[[cl.idx]])) ) {
		ACTIONet.out$annotations[[cl.idx]]$Enrichments = list()
	}
	cmd = sprintf('ACTIONet.out$annotations[[cl.idx]]$Enrichments$%s = Enrichment.out', enrichment.name)
	eval(parse(text=cmd))

	return(ACTIONet.out)
}


ATACtion.archetype.geneset.enrichment <- function(ACTIONet.out, sce, genesets, enrichment.name = NULL, core = T) {
	if(core == T) {
		if( !("unification.out" %in% names(ACTIONet.out)) ) {
			print("Running unification")
			ACTIONet.out$unification.out = unify.cell.states(ACTIONet.out, sce, reduction_slot = "ACTION", sce.data.attr = sce.data.attr)    
		} 
		DE.profile = ACTIONet.out$unification.out$DE.core		
	} else {
		if( !("unification.out" %in% names(ACTIONet.out)) ) {
			print("Running archetype DE")
			ACTIONet.out$archetype.differential.signature = compute.archetype.feature.specificity(ACTIONet.out, sce, mode = "dense", sce.data.attr = sce.data.attr)
		}
		DE.profile = ACTIONet.out$archetype.differential.signature
	}
	values(DE.profile) = c()

	if(is.sce(sce)) {
		cisConnectome = metadata(sce)[["cisConnectome"]]
	} else {
		cisConnectome = sce
	}
	
    require(Matrix)
    if (is.matrix(genesets) | is.sparseMatrix(genesets)) {    
		genesets = apply(genesets, 2, function(x) intersect(colnames(cisConnectome), rownames(genesets)[x > 0]))
	}
	annotations = as(sapply(genesets, function(gs) as.numeric(colnames(cisConnectome) %in% gs)), "sparseMatrix")
	rownames(annotations) = colnames(cisConnectome)

	X = cisConnectome
	Y = annotations

	ind.mat = as(X %*% Y, 'sparseMatrix')
	ind.mat@x = rep(1, length(ind.mat@x))
	
	scores = log1p(DE.profile@assays[["significance"]])
	
	
	counts = Matrix::rowSums(ind.mat)
	mask = counts > 0
    Enrichment.profile = Chernoff.enrichment.noRowScaling(scores[mask, ], ind.mat[mask, ])
    
	Enrichment.out = list(Enrichment.profile = Enrichment.profile, Ind.mat = ind.mat)
	
		
	if( is.null(enrichment.name) ) {
		enrichment.name = sprintf("Pathways (%s)", as.character(Sys.time()))
	}
	if(core == T) {
		if( !("Enrichments" %in% names(ACTIONet.out$unification.out)) ) {
			ACTIONet.out$unification.out$Enrichments = list()
		}
		cmd = sprintf('ACTIONet.out$unification.out$Enrichments$\"%s\" = Enrichment.out', enrichment.name)
		eval(parse(text=cmd))
	} else {
		if( !("Enrichments" %in% names(ACTIONet.out)) ) {
			ACTIONet.out$Enrichments = list()
		}
		cmd = sprintf('ACTIONet.out$Enrichments$\"%s\" = Enrichment.out', enrichment.name)
		eval(parse(text=cmd))		
	}
	
	return(ACTIONet.out)
}


ATACtion.cluster.impute.gene.expression <- function(ACTIONet.out, sce, genesets, annotation.cluster, enrichment.name = NULL, sce.data.attr = "bin_counts") {

	cl.idx = which(names(ACTIONet.out$annotations) == annotation.cluster)
	if(length(cl.idx) == 0) {
		R.utils::printf('Error in ATACtion.cluster.impute.gene.expression: annotation.cluster "%s" not found\n', annotation.cluster)
		return(ACTIONet.out)
	}		
	if(is.null(ACTIONet.out$annotations[[cl.idx]]$DE.profile)) {
		ACTIONet.out = compute.annotations.feature.specificity(ACTIONet.out, sce, annotation.name = annotation.cluster, sce.data.attr = sce.data.attr)
	}
	DE.profile = ACTIONet.out$annotations[[cl.idx]]$DE.profile

	values(DE.profile) = c()

	if(is.sce(sce)) {
		cisConnectome = metadata(sce)[["cisConnectome"]]
	} else {
		cisConnectome = sce
	}
	
    require(Matrix)
    if (is.matrix(genesets) | is.sparseMatrix(genesets)) {
		genesets = apply(genesets, 2, function(x) intersect(colnames(cisConnectome), rownames(genesets)[x > 0]))
    }
	annotations = as(sapply(genesets, function(gs) as.numeric(colnames(cisConnectome) %in% gs)), "sparseMatrix")
rownames(annotations) = colnames(cisConnectome)

	X = cisConnectome
	Y = annotations
	
	
	ind.mat = as(X %*% Y, 'sparseMatrix')
	ind.mat@x = rep(1, length(ind.mat@x))

	scores = log1p(DE.profile@assays[["significance"]])
	
	counts = Matrix::rowSums(ind.mat)
	mask = counts > 0
    Enrichment.profile = Chernoff.enrichment.noRowScaling(scores[mask, ], ind.mat[mask, ])
    
	Enrichment.out = list(Enrichment.profile = Enrichment.profile, Ind.mat = ind.mat)
	
	if( is.null(enrichment.name) ) {
		enrichment.name = db.name
	}
	if( !("Enrichments" %in% names(ACTIONet.out$annotations[[cl.idx]])) ) {
		ACTIONet.out$annotations[[cl.idx]]$Enrichments = list()
	}
	cmd = sprintf('ACTIONet.out$annotations[[cl.idx]]$Enrichments$%s = Enrichment.out', enrichment.name)
	eval(parse(text=cmd))

	return(ACTIONet.out)
}

ATACtion.archetype.motif.enrichment <- function(ACTIONet.out, sce, matched.motifs, enrichment.name = "motif.enrichment", core = T) {
	if(core == T) {
		if( !("unification.out" %in% names(ACTIONet.out)) ) {
			print("Running unification")
			ACTIONet.out$unification.out = unify.cell.states(ACTIONet.out, sce, reduction_slot = "ACTION", sce.data.attr = sce.data.attr)    
		} 
		DE.profile = ACTIONet.out$unification.out$DE.core		
	} else {
		if( !("unification.out" %in% names(ACTIONet.out)) ) {
			print("Running archetype DE")
			ACTIONet.out$archetype.differential.signature = compute.archetype.feature.specificity(ACTIONet.out, sce, mode = "dense", sce.data.attr = sce.data.attr)
		}
		DE.profile = ACTIONet.out$archetype.differential.signature
	}
	values(DE.profile) = c()
	scores = log1p(DE.profile@assays[["significance"]])

	
	if( ncol(matched.motifs) == nrow(scores) ) {
		matched.motifs = Matrix::t(matched.motifs)
	}

	motif.enrichment  = Chernoff.enrichment.noRowScaling(scores, matched.motifs)
    
		
	if( is.null(enrichment.name) ) {
		enrichment.name = sprintf("Pathways (%s)", as.character(Sys.time()))
	}
	if(core == T) {
		if( !("Enrichments" %in% names(ACTIONet.out$unification.out)) ) {
			ACTIONet.out$unification.out$Enrichments = list()
		}
		cmd = sprintf('ACTIONet.out$unification.out$Enrichments$\"%s\" = motif.enrichment', enrichment.name)
		eval(parse(text=cmd))
	} else {
		if( !("Enrichments" %in% names(ACTIONet.out)) ) {
			ACTIONet.out$Enrichments = list()
		}
		cmd = sprintf('ACTIONet.out$Enrichments$\"%s\" = Enrichment.out', enrichment.name)
		eval(parse(text=cmd))		
	}
	
	return(ACTIONet.out)
}


ATACtion.cluster.motif.enrichment <- function(ACTIONet.out,  annotation.cluster, matched.motifs, enrichment.name = "motif.enrichment", sce.data.attr = "bin_counts") {

	cl.idx = which(names(ACTIONet.out$annotations) == annotation.cluster)
	if(length(cl.idx) == 0) {
		R.utils::printf('Error in ATACtion.cluster.motif.enrichment: annotation.cluster "%s" not found\n', annotation.cluster)
		return(ACTIONet.out)
	}		
	if(is.null(ACTIONet.out$annotations[[cl.idx]]$DE.profile)) {
		ACTIONet.out = compute.annotations.feature.specificity(ACTIONet.out, sce, annotation.name = annotation.cluster, sce.data.attr = sce.data.attr)
	}
	scores = log1p(ACTIONet.out$annotations[[cl.idx]]$DE.profile@assays[["significance"]])


	if( ncol(matched.motifs) == nrow(scores) ) {
		matched.motifs = Matrix::t(matched.motifs)
	}

	motif.enrichment = Chernoff.enrichment(scores, matched.motifs)
	
	if( is.null(enrichment.name) ) {
		enrichment.name = db.name
	}
	if( !("Enrichments" %in% names(ACTIONet.out$annotations[[cl.idx]])) ) {
		ACTIONet.out$annotations[[cl.idx]]$Enrichments = list()
	}
	cmd = sprintf('ACTIONet.out$annotations[[cl.idx]]$Enrichments$%s = motif.enrichment ', enrichment.name)
	eval(parse(text=cmd))

	return(ACTIONet.out)
}

ATACtion.archetype.gene.expression.imputation <- function(ACTIONet.out, sce, enrichment.name = "ImputedExpression", core = T) {
	if(core == T) {
		if( !("unification.out" %in% names(ACTIONet.out)) ) {
			print("Running unification")
			ACTIONet.out$unification.out = unify.cell.states(ACTIONet.out, sce, reduction_slot = "ACTION", sce.data.attr = sce.data.attr)    
		} 
		DE.profile = ACTIONet.out$unification.out$DE.core		
	} else {
		if( !("unification.out" %in% names(ACTIONet.out)) ) {
			print("Running archetype DE")
			ACTIONet.out$archetype.differential.signature = compute.archetype.feature.specificity(ACTIONet.out, sce, mode = "dense", sce.data.attr = sce.data.attr)
		}
		DE.profile = ACTIONet.out$archetype.differential.signature
	}
	values(DE.profile) = c()

	if(is.sce(sce)) {
		cisConnectome = metadata(sce)[["cisConnectome"]]
	} else {
		cisConnectome = sce
	}
	linked.genes.count = Matrix::rowSums(cisConnectome)
	mask = linked.genes.count > 0
	
	scores = log1p(DE.profile@assays[["significance"]])
	
    RNA.profile = Chernoff.enrichment.noRowScaling(scores[mask, ], cisConnectome[mask, ])
    


		
	if( is.null(enrichment.name) ) {
		enrichment.name = sprintf("Pathways (%s)", as.character(Sys.time()))
	}
	if(core == T) {
		if( !("Enrichments" %in% names(ACTIONet.out$unification.out)) ) {
			ACTIONet.out$unification.out$Enrichments = list()
		}
		cmd = sprintf('ACTIONet.out$unification.out$Enrichments$\"%s\" = RNA.profile', enrichment.name)
		eval(parse(text=cmd))
	} else {
		if( !("Enrichments" %in% names(ACTIONet.out)) ) {
			ACTIONet.out$Enrichments = list()
		}
		cmd = sprintf('ACTIONet.out$Enrichments$\"%s\" = Enrichment.out', enrichment.name)
		eval(parse(text=cmd))		
	}
	
	return(ACTIONet.out)
}


ATACtion.cluster.impute.gene.expression <- function(ACTIONet.out, sce, annotation.cluster, enrichment.name = "ImputedExpression", sce.data.attr = "bin_counts") {

	cl.idx = which(names(ACTIONet.out$annotations) == annotation.cluster)
	if(length(cl.idx) == 0) {
		R.utils::printf('Error in ATACtion.cluster.impute.gene.expression: annotation.cluster "%s" not found\n', annotation.cluster)
		return(ACTIONet.out)
	}		
	if(is.null(ACTIONet.out$annotations[[cl.idx]]$DE.profile)) {
		ACTIONet.out = compute.annotations.feature.specificity(ACTIONet.out, sce, annotation.name = annotation.cluster, sce.data.attr = sce.data.attr)
	}
	DE.profile = ACTIONet.out$annotations[[cl.idx]]$DE.profile

	values(DE.profile) = c()

	if(is.sce(sce)) {
		cisConnectome = metadata(sce)[["cisConnectome"]]
	} else {
		cisConnectome = sce
	}
	linked.genes.count = Matrix::rowSums(cisConnectome)
	mask = linked.genes.count > 0
	
	scores = log1p(DE.profile@assays[["significance"]])
	
    RNA.profile = Chernoff.enrichment.noRowScaling(scores[mask, ], cisConnectome[mask, ])

	if( is.null(enrichment.name) ) {
		enrichment.name = db.name
	}
	if( !("Enrichments" %in% names(ACTIONet.out$annotations[[cl.idx]])) ) {
		ACTIONet.out$annotations[[cl.idx]]$Enrichments = list()
	}
	cmd = sprintf('ACTIONet.out$annotations[[cl.idx]]$Enrichments$%s = RNA.profile', enrichment.name)
	eval(parse(text=cmd))

	return(ACTIONet.out)
}


identify.associated.REs <- function(scores, peakset.ind.mat, pval.threshold = 0.001, score.threshold = 3) {
  perm = order(scores, decreasing = TRUE)

  cap = sum(scores >= score.threshold)
  
  a = scores[perm]
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
    
    max.val = max(logPvals[1:cap])
    max.idx = which.max(logPvals[1:cap])
  
    if(max.val > log(pval.threshold) + log(length(X))) { # Bonferroni correction
		print(j)
		print(max.val)
		
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


ATACtionet.GWAS.enrichment <- function(A, R, GRList, flank = 250, A.slot = "signature", log.transform = TRUE) {
	if(is.matrix(A) | is.sparseMatrix(A)) {
		A = A
	} else if(is.sce(A)) {
		A = A@assays[[A.slot]]
	} else {
		A = A$archetype.accessibility@assays[[A.slot]]
	}
	if( (log.transform == TRUE) & (min(A) >= 0) )
		A = log(1 + A)

  
  if(is.GR(R)) {
	  GR = R
  } else {
	GR = rowRanges(R)
  }  
  values(GR) = c()
  
  src.genome = genome(GR)[[1]]
  if(src.genome == 'GRCh38') {
    src.genome = 'hg38'
    genome(GR) = 'hg38'
  }
  dst.genome = genome(GRList[[1]])[[1]]
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


getActiveEnhancers <- function(ACTIONet.out, sce, EnhEnrich.out, enrichment.z.threshold = 1.96, score.threshold = 3, core.only = TRUE) {
  Enrich.profile = EnhEnrich.out$Enrichment.profile
  scores = ACTIONet.out$archetype.accessibility@assays[["signature"]]

  
  if(core.only) {
    Enrich.profile = Enrich.profile[, ACTIONet.out$core.out$core.archs]
    scores = scores[, ACTIONet.out$core.out$core.archs]
  }
  
  GR = rowRanges(sce)
  values(GR) = c()
  
  peakset.ind.mat = EnhEnrich.out$Ind.mat
  
  
  Enh.modules.GList = GenomicRangesList(sapply(1:ncol(Enrich.profile), function(j) {
    x = scores[, j]
    rows = order(x, decreasing = TRUE)[1:sum(x >= score.threshold)]
    
    enrichment.z = scale(Enrich.profile[, j])
    cols = which(enrichment.z >= enrichment.z.threshold)
    
    if(length(cols) == 0)
    return()
    
    if(length(cols) == 1) {
      idx = which(peakset.ind.mat[rows, cols] > 0)
      Enh.GR = GR[rows[idx]]
    } else {
      idx = which(Matrix::rowSums(peakset.ind.mat[rows, cols]) > 0)
      Enh.GR = GR[rows[idx]]
    }  
  }))
  genome(Enh.modules.GList) = genome(GR)
  names(Enh.modules.GList) = ACTIONet.out$core.out$core.archs
  
  return(Enh.modules.GList)
}


