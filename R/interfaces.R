run_SnapATAC <- function(x.sp, reduced_dim = 50, min.cell = 0, thread_no = 8, knn = 15, seed = 0) {
	library(SnapATAC)
	x.sp = makeBinary(x.sp, mat="pmat");

	if(reduction_method == "LDA") {
		x.sp = runLDA(obj=x.sp,
			input.mat="pmat",
			method = "Probability",
			topic=reduced_dim, 
			min.cell=min.cell,
			num.cores=thread_no,
			seed.use=seed)
	} else if(reduction_method == "JDA") {
		x.sp = runJDA(
			obj=x.sp,
			input.mat="pmat",
			norm.method="normOVE",
			bin.cov.zscore.lower=-2,
			bin.cov.zscore.upper=2,
			pc.num=reduced_dim,
			max.var=5000,
			do.par=TRUE,
			ncell.chunk=1000,
			num.cores=thread_no,
			seed.use=seed,
			tmp.folder=tempdir()
			);		
	} else if(reduction_method == "LSA") {
		x.sp = runLSA(
			obj=x.sp, 
			input.mat="pmat", 
			pc.num=reduced_dim, 
			logTF=FALSE,
			min.cell=min.cells, 
			seed.use=seed
			);		
	} else if(reduction_method == "LDAlogTF") {
		x.sp = runLSA(
			obj=x.sp, 
			input.mat="pmat", 
			pc.num=reduced_dim, 
			logTF=TRUE,
			min.cell=min.cells, 
			seed.use=seed
			);		
	} else {
		stop(sprintf("Unknown reduction method %s", reduction_method))
		return()
	}


    x.sp = runKNN(
        obj=x.sp,
        pca.dims=1:reduced_dim,
        weight.by.sd=TRUE,
        k=knn
	)

	return(x.sp)
	
}

extract_SnapATAC_reduction <- function(x.sp) {
	reduction = x.sp@smat@dmat %*% diag(x.sp@smat@sdev)	
	
	return(reduction)
}

run_cisTopic <- function(cisTopicObject, topic.counts = seq(5, 50, by = 5), seed=0, thread_no = 8, burnin = 120, iterations = 150, addModels=FALSE, project.name = "cisTopic") {	
	cisTopicObject <- runModels(cisTopicObject, topic=topic.counts, seed=seed, nCores=thread_no, burnin = burnin, iterations = iterations, addModels=FALSE)
	
	return(cisTopicObject)
}

extract_cisTopic_factors <- function(cisTopicObject) {
	W = as.matrix(cisTopic:::.getScores(cisTopicObject))
	H <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
	
	out = list(W = W, H = H)
	
	return(out)
}
