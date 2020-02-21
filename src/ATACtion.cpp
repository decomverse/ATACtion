#include <ATACtion.h>
#include <RcppArmadillo.h>

using namespace Rcpp;

#define ARMA_USE_CXX11_RNG


// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List reduceChromatinAccessibility(sp_mat expression, int reduced_dim = 30, int method=1, int iters = 10) {
	Projection projection = ATACtion::reduceChromatinAccessibility(expression, reduced_dim, method, iters);

	List res;	
	res["S_r"] = projection.S_r;		
	res["V"] = projection.V;
	res["lambda"] = projection.lambda;
	res["explained_var"] = projection.exp_var;	
	
	return res;
}
