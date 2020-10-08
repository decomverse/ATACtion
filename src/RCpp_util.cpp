#include <RcppArmadillo.h>
#include <RCpp_util.h>
#include <ACTIONet.h>


#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#define ARMA_USE_CXX11_RNG
#define DYNSCHED

// From from Granja, et al. (2019)
// [[Rcpp::export]]
IntegerMatrix tabulate2dCpp(IntegerVector x1, int xmin, int xmax, IntegerVector y1, int ymin, int ymax){
  if(x1.size() != y1.size()){
    stop("width must equal size!");
  }
  IntegerVector x = clone(x1);
  IntegerVector y = clone(y1);
  int n = x.size();
  IntegerVector rx = seq(xmin,xmax);
  IntegerVector ry = seq(ymin,ymax);
  IntegerMatrix mat( ry.size() , rx.size() );
  int xi,yi;
  for(int i = 0; i < n; i++){
    xi = (x[i] - xmin);
    yi = (y[i] - ymin);
    if(yi >= 0 && yi < ry.size()){
      if(xi >= 0 && xi < rx.size()){
        mat( yi , xi ) = mat( yi , xi ) + 1;
      }
    }
  }
  return mat;
}
