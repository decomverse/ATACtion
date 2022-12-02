#ifndef ATACtion_H
#define ATACtion_H

#include <stdio.h>
#include <vector>
#include <string>
#include <cmath>
#include <sys/types.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <errno.h>
#include <sys/stat.h>
#include <getopt.h>
#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <unordered_map>
#include <ctime>


#define ARMA_DONT_USE_WRAPPER
#undef ARMA_BLAS_CAPITALS
#define ARMA_BLAS_UNDERSCORE

/*
#undef ARMA_BLAS_CAPITALS
#define ARMA_BLAS_UNDERSCORE
#define ARMA_64BIT_WORD
#define ARMA_BLAS_LONG_LONG
*/

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace arma;
using namespace std;

#define DEFAULT_PCA_DIM 50

// Expression reduction methods
#define ACTIONplusPCA 1

// Generic structure returned by all reduction methods
struct Projection {
	mat S_r;
	mat V;
	vec lambda;
	vec exp_var;
};

namespace ATACtion {
	Projection reducedKernel(sp_mat &profile, int PCA_dim, int iter, int seed);
	Projection reduceChromatinAccessibility(sp_mat &peaks, int reduced_dim, int method, int iter);
}

#endif
