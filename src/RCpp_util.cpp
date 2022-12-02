#include <RcppArmadillo.h>
#include <RCpp_util.h>
#include <ATACtion.h>


// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

#define ARMA_USE_CXX11_RNG
#define DYNSCHED

// From from Granja, et al. (2019)
// [[Rcpp::export]]
Rcpp::IntegerMatrix tabulate2dCpp(Rcpp::IntegerVector x1, int xmin, int xmax, Rcpp::IntegerVector y1, int ymin, int ymax){
  if(x1.size() != y1.size()){
    stop("width must equal size!");
  }
  Rcpp::IntegerVector x = clone(x1);
  Rcpp::IntegerVector y = clone(y1);
  int n = x.size();
  Rcpp::IntegerVector rx = seq(xmin,xmax);
  Rcpp::IntegerVector ry = seq(ymin,ymax);
  Rcpp::IntegerMatrix mat( ry.size() , rx.size() );
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
