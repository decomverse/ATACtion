#include "ATACtion.h"


namespace ATACtion {
	Projection reduceChromatinAccessibility(sp_mat &peaks, int reduced_dim = DEFAULT_PCA_DIM, int method = ACTIONplusPCA, int iter = 10) {
		printf("Reducing expression matrix\n");
		Projection projection;

		peaks = spones(peaks);
		
		switch(method) {
			case ACTIONplusPCA:	// Uses dense formulation
			{
				printf("\tReduce expression matrix using orthogonalization followed by PCA (k = %d) using sparse formulation ... \n", reduced_dim); fflush(stdout);
				projection = reducedKernel(peaks, reduced_dim, iter, 0);		
				printf("done\n"); fflush(stdout);
				
			}
			break;
				
			default:
				fprintf(stderr, "Unknown ATAC reduction method code %d\n", method);
		}

		return projection;
	}	
}
